
# -*- coding: utf-8 -*-
"""
Python script to intersect hazards and network line geometries
Created on Wed Jul 18 2018

@author: Raghav Pant
"""

import os
import sys
import pandas as pd
import geopandas as gpd
import itertools
from shapely.geometry import Polygon

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..','..'))
from scripts.utils import load_config, line_length

def networkedge_hazard_intersection(edge_shapefile,hazard_shapefile,output_shapefile):
	line_gpd = gpd.read_file(edge_shapefile)
	poly_gpd = gpd.read_file(hazard_shapefile)

	if len(line_gpd.index) > 0 and len(poly_gpd.index) > 0:
		line_gpd.columns = map(str.lower, line_gpd.columns)
		poly_gpd.columns = map(str.lower, poly_gpd.columns)

		line_bounding_box = line_gpd.total_bounds
		line_bounding_box_coord = list(itertools.product([line_bounding_box[0],line_bounding_box[2]], [line_bounding_box[1],line_bounding_box[3]]))
		line_bounding_box_geom = Polygon(line_bounding_box_coord)
		line_bounding_box_gpd = gpd.GeoDataFrame(pd.DataFrame([[1],[line_bounding_box_geom]]).T,crs='epsg:4326')
		line_bounding_box_gpd.columns = ['ID','geometry']

		poly_bounding_box = poly_gpd.total_bounds
		poly_bounding_box_coord = list(itertools.product([poly_bounding_box[0],poly_bounding_box[2]], [poly_bounding_box[1],poly_bounding_box[3]]))
		poly_bounding_box_geom = Polygon(poly_bounding_box_coord)
		poly_bounding_box_gpd = gpd.GeoDataFrame(pd.DataFrame([[1],[poly_bounding_box_geom]]).T,crs='epsg:4326')
		poly_bounding_box_gpd.columns = ['ID','geometry']

		poly_sindex = poly_bounding_box_gpd.sindex

		selected_polys = poly_bounding_box_gpd.iloc[list(poly_sindex.intersection(line_bounding_box_gpd['geometry'].iloc[0].bounds))]
		if len(selected_polys.index) > 0:
			data = []
			#create spatial index
			poly_sindex = poly_gpd.sindex
			for l_index, lines in line_gpd.iterrows():
				intersected_polys = poly_gpd.iloc[list(poly_sindex.intersection(lines.geometry.bounds))]
				for p_index, poly in intersected_polys.iterrows():
					if (lines['geometry'].intersects(poly['geometry']) is True) and (poly.geometry.is_valid is True):
						# data.append({'edge_id': lines['edge_id'],'width':lines['width'],'length':1000.0*line_length(lines['geometry'].intersection(poly['geometry'])),'geometry':lines['geometry'].intersection(poly['geometry'])})
						data.append({'edge_id': lines['edge_id'],'length':1000.0*line_length(lines['geometry'].intersection(poly['geometry'])),'geometry':lines['geometry'].intersection(poly['geometry'])})
			if data: 
				# intersections_data = gpd.GeoDataFrame(data,columns=['edge_id','width','length','geometry'],crs='epsg:4326')
				intersections_data = gpd.GeoDataFrame(data,columns=['edge_id','length','geometry'],crs='epsg:4326')
				# intersections_data['area'] = intersections_data['width']*intersections_data['length']
				intersections_data.to_file(output_shapefile)

				del intersections_data

	del line_gpd, poly_gpd


def main():
	data_path = load_config()['paths']['data']
	output_path = load_config()['paths']['output']

	provinces = ['Thanh Hoa','Binh Dinh','Lao Cai']

	for province in provinces:
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()
		
		road_shp = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))
		hazard_dir = os.path.join(data_path,'Hazard_data','Glofris')
		for root, dirs, files in os.walk(hazard_dir):
			for file in files:
				if file.endswith(".shp"):
					hazard_shp = os.path.join(root,file)
					out_shp_name = 'vietbando_{}_edges_'.format(province_name) + file
					out_shp = os.path.join(output_path,'{}_roads_hazard_intersections'.format(province_name),out_shp_name)

					networkedge_hazard_intersection(road_shp,hazard_shp,out_shp)

					print ('Done with vietbando_{0}_edges.shp and {1}'.format(province_name,file))

if __name__ == "__main__":
	main()



