# -*- coding: utf-8 -*-
"""
Created on Tue July 17 10:20:55 2018

@authors: Raghav Pant, Tom Russell, elcok
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
	with rasterio.open(file_path) as dataset:
		counts = dataset.count
		if dataset.crs:
			crs = dataset.crs.to_string()
		else:
			crs = 'invalid/unknown'
		# bands = dataset.meta
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

		# print (list(set(resphaped_array.tolist())))
	return counts,crs, data_vals

def convert_geotiff_to_vector_with_threshold(from_threshold,to_threshold, infile, infile_epsg,tmpfile_1, tmpfile_2, outfile):
	"""Threshold raster, convert to polygons, assign crs
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
	"""Threshold raster, convert to polygons, assign crs
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
	"""Threshold raster, convert to polygons
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
	data_path = load_config()['paths']['data']
	root_dir = os.path.join(data_path,'Hazard_data','Glofris')
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
					thresholds = [1,2,3,4,999]
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
					if file == 'LSZ_NgheAn_to_PhuYen.tif':
						code_vals = [3,4]
					else:
						code_vals = [4,5]
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
						if dv in [(255,190,190),(245,0,0),(255,0,0)]:
							thr = 5
							bc = dv
							in_file = os.path.join(root,file)
							tmp_1 = os.path.join(root,file.split(".tif")[0] + '_mask.tiff')
							tmp_2 = os.path.join(root,file.split(".tif")[0] + '_mask.shp')
							out_file = os.path.join(root,file.split(".tif")[0] + '_{}_band.shp'.format(thr))
							convert_geotiff_to_vector_with_multibands(bc,in_file,s_crs,tmp_1, tmp_2, out_file)
						elif dv in  [(255,170,0),(255,128,0)]:
							thr = 4
							bc = dv

							in_file = os.path.join(root,file)
							tmp_1 = os.path.join(root,file.split(".tif")[0] + '_mask.tiff')
							tmp_2 = os.path.join(root,file.split(".tif")[0] + '_mask.shp')
							out_file = os.path.join(root,file.split(".tif")[0] + '_{}_band.shp'.format(thr))
							convert_geotiff_to_vector_with_multibands(bc,in_file,s_crs,tmp_1, tmp_2, out_file)





if __name__ == "__main__":
	main()
