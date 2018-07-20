
# -*- coding: utf-8 -*-
"""
Python script to intersect hazards and network line geometries
Created on Wed Jul 18 2018

@author: Raghav Pant

Create tables of all failure scenarios with the follwing attributes 
Edge_id Hazard_type year climate_scenario hazard_band band_name min_depth max_depth affect_length commune_id commune_name district_id district_name province_id province_name 
"""

import os
import sys
import pandas as pd
import geopandas as gpd
import itertools
from shapely.geometry import Polygon

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..','..'))
from scripts.utils import load_config, line_length

def all_failure_edges(input_shapefile):
	shape_gpd = gpd.read_file(input_shapefile)
	shape_vals = pd.DataFrame(shape_gpd.values)
	del shape_gpd
	fail_edges = shape_vals.groupby(['edge_id'])['length'].sum().reset_index()

	return fail_edges

def spatial_scenario_selection(line_shapefile,polygon_shapefile,name_province,hazard_dictionary,data_dictionary):
	line_gpd = gpd.read_file(line_shapefile)
	poly_gpd = gpd.read_file(polygon_shapefile)

	poly_gpd = poly_gpd[poly_gpd['PRO_NAME_E'] == name_province]


	if len(line_gpd.index) > 0 and len(poly_gpd.index) > 0:
		line_gpd.columns = map(str.lower, line_gpd.columns)
		poly_gpd.columns = map(str.lower, poly_gpd.columns)

		#create spatial index
		poly_sindex = poly_gpd.sindex

		poly_sindex = poly_gpd.sindex
		for l_index, lines in line_gpd.iterrows():
			intersected_polys = poly_gpd.iloc[list(poly_sindex.intersection(lines.geometry.bounds))]
			for p_index, poly in intersected_polys.iterrows():
				if (lines['geometry'].intersects(poly['geometry']) is True) and (poly.geometry.is_valid is True):
					value_dictionary = {'edge_id': lines['edge_id'],'length':lines['length'],
							'province_id':poly['province_i'],'province_name':poly['pro_name_e'],
							'district_id':poly['district_i'],'district_name':poly['dis_name_e'],
							'commune_id':poly['commune_id'],'commune_name':poly['name_eng']}

					data_dictionary.append({**value_dictionary,**hazard_dictionary})

	del line_gpd, poly_gpd
	return data_dictionary


def main():
	data_path = load_config()['paths']['data']
	hazard_intersections_path = load_config()['paths']['output']
	provinces = ['Thanh Hoa','Binh Dinh','Lao Cai']
	bnds = [3,4,5]
	thresholds = [1,2,3,4,999]
	sector = 'Roads'

	commune_shp = os.path.join(data_path,'Vietnam_boundaries','who_boundaries','who_communes.shp')

	hazard_description_file = os.path.join(data_path,'Hazard_data','hazard_data_folder_data_info.xlsx')
	hazard_df = pd.read_excel(hazard_description_file,sheet_name ='file_contents')
	hazard_files = hazard_df['file_name'].values.tolist()

	data_excel = os.path.join(hazard_intersections_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')
	excel_writer = pd.ExcelWriter(data_excel)
	for province in provinces:
		# set all paths for all input files we are going to use
		data_dict = []
		province_name = province.replace(' ','').lower()
		intersection_dir = os.path.join(hazard_intersections_path,'{}_roads_hazard_intersections'.format(province_name))
		for root, dirs, files in os.walk(intersection_dir):
			for file in files:
				if file.endswith(".shp"):
					hazard_dict = {}
					hazard_dict['sector'] = sector
					hazard_shp = os.path.join(root,file)
					hz_file = [h for h in hazard_files if h in file][0]
					hazard_dict['hazard_type'] = hazard_df.loc[hazard_df.file_name == hz_file].hazard_type.values[0]
					hazard_dict['year'] = hazard_df.loc[hazard_df.file_name == hz_file].year.values[0]
					hazard_dict['climate_scenario'] = hazard_df.loc[hazard_df.file_name == hz_file].climate_scenario.values[0]
					hazard_dict['probability'] = hazard_df.loc[hazard_df.file_name == hz_file].probability.values[0]
					band_type = hazard_df.loc[hazard_df.file_name == hz_file].banded.values[0]
					if str(band_type) == 'True':
						hazard_dict['band_num'] = [b for b in bnds if '{}_band'.format(b) in file][0]
						band_names = hazard_df.loc[hazard_df.file_name == hz_file].bands.values[0].split(',')
						hazard_dict['band_name'] = [b for b in band_names if str(hazard_dict['band_num']) in b][0]
						hazard_dict['min_val'] = 0
						hazard_dict['max_val'] = 0
					else:
						hazard_dict['band_num'] = 0
						hazard_dict['band_name'] = 'none'
						hazard_thrs = [(thresholds[t],thresholds[t+1]) for t in range(len(thresholds)-1) if '{0}m-{1}m'.format(thresholds[t],thresholds[t+1]) in file][0]
						hazard_dict['min_val'] = hazard_thrs[0]
						hazard_dict['max_val'] = hazard_thrs[1]


					data_dict = spatial_scenario_selection(hazard_shp,commune_shp,province,hazard_dict,data_dict)
					print ('Done with',file)
		
		data_df = pd.DataFrame(data_dict)
		data_df_cols = data_df.columns.values.tolist()
		selected_cols = [cols for cols in data_df_cols if cols != 'length']
		data_df = data_df.groupby(selected_cols)['length'].sum().reset_index()
		data_df.to_excel(excel_writer,province_name,index = False)
		excel_writer.save()
		del data_df



if __name__ == "__main__":
	main()



