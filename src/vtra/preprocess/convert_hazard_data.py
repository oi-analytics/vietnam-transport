"""
Purpose
-------

Convert GeoTiff raster hazard datasets to shapefiles based on masking and selecting values from
    - Single-band raster files
    - Multi-band (3-bands) raster files

Input data requirements
-----------------------

1. Correct paths to all hazard datasets
2. Single-band GeoTiff hazard raster files with:
    - values - between 0 and 1000
    - raster grid geometry
    - projection systems: Default assumed = EPSG:32648

3. Multi-band GeoTiff hazard raster files with:
    - 3-bands
    - values - in each band between 0 and 255
    - raster grid geometry
    - projection systems: Default assumed = EPSG:32648

Results
-------

1. Shapefiles whose names show the hazard models and their selected range of values 
    - ID - equal to 1 
    - geometry - Shapely Polygon outline of selected hazard
"""

import os
import subprocess
import json
import sys


from vtra.utils import load_config

import fiona
import fiona.crs
import rasterio
import numpy as np
import pandas as pd

def glofris_data_details(file_name,root_dir):
	"""
    Read names of GLOFRIS files and create attributes 

    Parameters
        - file_name - String name of GeoTff file
        - root_dir - String path to directory of file

    Outputs
        df - Pandas DataFrame written to csv file with columns:
            - file_name - String
            - hazard_type - String
            - year - Integer: 2016 or 2030
            - climate_scenario - String: RCP4.5 or RCP8.5 or none
            - probability - Float: 1/(return period)
            - banded - Boolean: True or False
            - bands - Integer
    """
	for root, dirs, files in os.walk(root_dir):
		for file in files:
			if file.endswith(".tif") or file.endswith(".tiff"):
				fname = file.split('.tif')
				fname = fname[0]
				print (fname)
				if '2030' in fname:
					year = 2030
				else:
					year = 2016

				if 'rcp4p5' in fname:
					sc = 'rcp 4.5'
				elif 'rcp8p5' in fname:
					sc = 'rcp 8.5'
				else:
					sc = 'none'

				f_all.append((fname,'flooding',year,sc,1.0/float(fname[-5:]),'FALSE','none'))

	df = pd.DataFrame(f_all,columns = ['file_name',	'hazard_type',	'year',	'climate_scenario',	'probability','banded',	'bands'])
	df.to_csv(os.path.join(root_dir,'glofris_files.csv'),index = False)



def raster_rewrite(in_raster,out_raster,nodata):
	"""
    Rewrite a raster to reproject and change no data value

    Parameters
        - in_raster - String name of input GeoTff file path
        - out_raster - String name of output GeoTff file path
        - nodata - Float value of data that is treated as no data 

    Outputs
        Reproject and replace raster with nodata = -1  
    """
	with rasterio.open(in_raster) as dataset:
		data_array = dataset.read()
		data_array[np.where(np.isnan(data_array))] = nodata

		with rasterio.open(out_raster, 'w', driver='GTIff',
					height=data_array.shape[1],    # numpy of rows
					width=data_array.shape[2],     # number of columns
					count=dataset.count,                        # number of bands
					dtype=data_array.dtype,  # this must match the dtype of our array
					crs=dataset.crs,
					transform=dataset.transform) as out_data:
			out_data.write(data_array)  # optional second parameter is the band number to write to
			out_data.nodata = -1  # set the raster's nodata value


	os.remove(in_raster)
	os.rename(out_raster,in_raster)


def raster_projections_and_databands(file_path):
	"""
    Extract projection, data bands numbers and valuees from raster

    Parameters
        - file_path - String name of input GeoTff file path

    Outputs
        - counts - Number of bans in raster
        - crs - Projection system of raster
        - data_vals - Numpy array of raster values 
    """
	with rasterio.open(file_path) as dataset:
		counts = dataset.count
		if dataset.crs:
			crs = dataset.crs.to_string()
		else:
			crs = 'invalid/unknown'
		data_array = dataset.read()
		if dataset.count > 1:
			data_list = []
			for i in range(0,dataset.count):
				data_list.append(data_array[i].reshape(dataset.height*dataset.width).tolist())
			data_vals = list(set(list(zip(*data_list))))
		else:
			data_vals = list(set(data_array.reshape(dataset.count*dataset.height*dataset.width).tolist()))
			if all(isinstance(x, int) for x in data_vals) is False:
				data_vals = []

	return counts,crs, data_vals

def convert_geotiff_to_vector_with_threshold(from_threshold,to_threshold, infile, infile_epsg,tmpfile_1, tmpfile_2, outfile):
	"""
    Convert GeoTiff raster file to Shapefile with geometries based on raster threshold ranges  

    Parameters
        - from_threshold - Float value of lower bound of GeoTiff threshold value to be selected
        - to_threshold - Float value of upper bound of GeoTiff threshold value to be selected
        - infile - String name of input GeoTff file path
        - infile_epsg - Integer value of EPSG Projection number of raster
        - tmpfile_1 - Stirng name of tmp file 1
        - tmpfile_2 - Stirng name of tmp file 2
        - outfile - Stirng name of output shapefile

    Outputs
        Shapefile with Polygon geometries of rasters based on raster threshold ranges
    """
	args = [
		"gdal_calc.py",
		'-A', infile,
		'--outfile={}'.format(tmpfile_1),
		'--calc=logical_and(A>={0}, A<{1})'.format(from_threshold,to_threshold),
		'--type=Byte', '--NoDataValue=0',
		'--co=SPARSE_OK=YES',
		'--co=NBITS=1',
		'--quiet',
		'--co=COMPRESS=LZW'
	]
	subprocess.run(args)

	subprocess.run([
		"gdal_edit.py",
		'-a_srs', 'EPSG:{}'.format(infile_epsg),
		tmpfile_1
	])

	subprocess.run([
		"gdal_polygonize.py",
		tmpfile_1,
		'-q',
		'-f', 'ESRI Shapefile',
		tmpfile_2
	])

	subprocess.run([
		"ogr2ogr",
		'-a_srs', 'EPSG:{}'.format(infile_epsg),
		'-t_srs', 'EPSG:4326',
		outfile,
		tmpfile_2
	])

	subprocess.run(["rm", tmpfile_1])
	subprocess.run(["rm", tmpfile_2])
	subprocess.run(["rm", tmpfile_2.replace('shp', 'shx')])
	subprocess.run(["rm", tmpfile_2.replace('shp', 'dbf')])
	subprocess.run(["rm", tmpfile_2.replace('shp', 'prj')])

def convert_geotiff_to_vector_with_multibands(band_colors, infile, infile_epsg,tmpfile_1, tmpfile_2, outfile):
	"""
    Convert multi-band GeoTiff raster file to Shapefile with geometries based on raster band color values  

    Parameters
        - band_colors - Tuple with 3-values each corresponding to the values in raster bands 
        - infile - String name of input GeoTff file path
        - infile_epsg - Integer value of EPSG Projection number of raster
        - tmpfile_1 - Stirng name of tmp file 1
        - tmpfile_2 - Stirng name of tmp file 2
        - outfile - Stirng name of output shapefile

    Outputs
        Shapefile with Polygon geometries of rasters based on raster band values
    """
	args = [
		"gdal_calc.py",
		'-A', infile,
		'--A_band=1',
		'-B', infile,
		'--B_band=2',
		'-C', infile,
		'--C_band=3',
		'--outfile={}'.format(tmpfile_1),
		'--type=Byte', '--NoDataValue=0',
		'--calc=logical_and(A=={0}, B=={1},C=={2})'.format(band_colors[0],band_colors[1],band_colors[2]),
		'--co=SPARSE_OK=YES',
		'--co=NBITS=1',
		'--quiet',
		'--co=COMPRESS=LZW'
	]
	subprocess.run(args)

	subprocess.run([
		"gdal_edit.py",
		'-a_srs', 'EPSG:{}'.format(infile_epsg),
		tmpfile_1
	])


	subprocess.run([
		"gdal_polygonize.py",
		tmpfile_1,
		'-q',
		'-f', 'ESRI Shapefile',
		tmpfile_2
	])

	subprocess.run([
		"ogr2ogr",
		'-a_srs', 'EPSG:{}'.format(infile_epsg),
		'-t_srs', 'EPSG:4326',
		outfile,
		tmpfile_2
	])

	subprocess.run(["rm", tmpfile_1])
	subprocess.run(["rm", tmpfile_2])
	subprocess.run(["rm", tmpfile_2.replace('shp', 'shx')])
	subprocess.run(["rm", tmpfile_2.replace('shp', 'dbf')])
	subprocess.run(["rm", tmpfile_2.replace('shp', 'prj')])

def convert(threshold, infile, tmpfile_1, outfile):
	"""
    Convert GeoTiff raster file to Shapefile with geometries based on raster threshold less that 999  

    Parameters
        - threshold - Float value of lower bound of GeoTiff threshold value to be selected
        - infile - String name of input GeoTff file path
        - tmpfile_1 - Stirng name of tmp file 1
        - outfile - Stirng name of output shapefile

    Outputs
        Shapefile with Polygon geometries of rasters based on raster values above a threshold
    """
	args = [
		"gdal_calc.py",
		'-A', infile,
		'--outfile={}'.format(tmpfile_1),
		'--calc=logical_and(A>={}, A<999)'.format(threshold),
		'--type=Byte', '--NoDataValue=0',
		'--co=SPARSE_OK=YES',
		'--co=NBITS=1',
		'--co=COMPRESS=LZW'
	]
	subprocess.run(args)

	subprocess.run([
		"gdal_polygonize.py",
		tmpfile_1,
		'-q',
		'-f', 'ESRI Shapefile',
		outfile
	])

def main():
	"""
    1. Specify the paths from where to read and write:
        - Input data
        - Hazard data
    
    2. Supply input data and parameters
        - Thresholds of flood hazards
        - Values of bands to be selected
        - Color code of multi-band rasters
        - Specific file names that might require some specific operations
    """
	data_path = load_config()['paths']['data']
	root_dir = os.path.join(data_path,'Hazard_data')

	thresholds = [1,2,3,4,999]
	band_vals_1 = [3,4]
	band_vals_2 = [4,5]
	color_codes_1 = [(255,190,190),(245,0,0),(255,0,0)]
	color_codes_2 = [(255,170,0),(255,128,0)]
	specific_files = ['LSZ_NgheAn_to_PhuYen.tif'] 
	f_all = []
	for root, dirs, files in os.walk(root_dir):
		for file in files:
			if file.endswith(".tif") or file.endswith(".tiff"):
				band_nums, crs, unique_data_values = raster_projections_and_databands(os.path.join(root, file))
				print (file,crs, unique_data_values)
				if 'epsg' in crs:
					crs_split = crs.split(':')
					s_crs = [int(c) for c in crs_split if c.isdigit() is True][0]
				else:
					s_crs = 32648


				if not unique_data_values:
					# threshold based datasets
					for t in range(len(thresholds)-1):
						thr_1 = thresholds[t]
						thr_2 = thresholds[t+1]
						in_file = os.path.join(root,file)
						tmp_1 = os.path.join(root,file.split(".tif")[0] + '_mask.tiff')
						tmp_2 = os.path.join(root,file.split(".tif")[0] + '_mask.shp')
						out_file = os.path.join(root,file.split(".tif")[0] + '_{0}m-{1}m_threshold.shp'.format(thr_1,thr_2))
						convert_geotiff_to_vector_with_threshold(thr_1,thr_2,in_file,s_crs,tmp_1, tmp_2, out_file)
				elif band_nums == 1:
					# code value based dataset
					if file in specific_files:
						code_vals = band_vals_1
					else:
						code_vals = band_vals_2
					for c in code_vals:
						in_file = os.path.join(root,file)
						tmp_1 = os.path.join(root,file.split(".tif")[0] + '_mask.tiff')
						tmp_2 = os.path.join(root,file.split(".tif")[0] + '_mask.shp')
						out_file = os.path.join(root,file.split(".tif")[0] + '_{}_band.shp'.format(c))
						convert_geotiff_to_vector_with_threshold(c,c+1,in_file,s_crs,tmp_1, tmp_2, out_file)
				if band_nums == 3:
					# multi-band color datasets
					# remove nodata values from the bands in the raster
					raster_rewrite(os.path.join(root, file),os.path.join(root, 'test.tif'),0)

					for dv in unique_data_values:
						if dv in color_codes_1:
							thr = 5
							bc = dv
							in_file = os.path.join(root,file)
							tmp_1 = os.path.join(root,file.split(".tif")[0] + '_mask.tiff')
							tmp_2 = os.path.join(root,file.split(".tif")[0] + '_mask.shp')
							out_file = os.path.join(root,file.split(".tif")[0] + '_{}_band.shp'.format(thr))
							convert_geotiff_to_vector_with_multibands(bc,in_file,s_crs,tmp_1, tmp_2, out_file)
						elif dv in  color_codes_2:
							thr = 4
							bc = dv

							in_file = os.path.join(root,file)
							tmp_1 = os.path.join(root,file.split(".tif")[0] + '_mask.tiff')
							tmp_2 = os.path.join(root,file.split(".tif")[0] + '_mask.shp')
							out_file = os.path.join(root,file.split(".tif")[0] + '_{}_band.shp'.format(thr))
							convert_geotiff_to_vector_with_multibands(bc,in_file,s_crs,tmp_1, tmp_2, out_file)





if __name__ == "__main__":
	main()
