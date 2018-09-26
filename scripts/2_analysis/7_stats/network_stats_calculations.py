"""Summarise hazard data

Get OD data and process it
Author: Raghav Pant
Date: April 20, 2018
"""
import geopandas as gpd
import pandas as pd
import os
import igraph as ig
import numpy as np
import sys
import subprocess
from shapely.geometry import Point
from shapely.geometry import Polygon
from scipy.spatial import Voronoi
import itertools
import operator
import ast
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches



from vtra.utils import *
from vtra.transport_network_creation import *

mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 13.
mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 14.
mpl.rcParams['xtick.labelsize'] = 13.
mpl.rcParams['ytick.labelsize'] = 13.

'''
Create the database connection
'''
def main():
	data_path,calc_path,output_path,figure_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output'],load_config()['paths']['figures']

	'''
	Get the modal shares
	'''
	modes_cols = ['road','rail','air','inland','coastal']

	for m in range(len(modes_cols)):
		mode_data_path = os.path.join(data_path,'Results','Flow_shapefiles','weighted_edges_flows_national_{0}.shp'.format(modes_cols[m]))
		edges_df = gpd.read_file(mode_data_path)

		if modes_cols[m] in ('air','inland','coastal'):
			edges_df = edges_df[edges_df['max_tons'] > 0]

		length = sum(edges_df['length'].values.tolist())
		print ('Total length of {0} network = {1}'.format(modes_cols[m],length))

		if modes_cols[m] == 'road':
			paved_df = edges_df[edges_df['road_cond'] == 'paved']
			paved_length = sum(paved_df['length'].values.tolist())

			unpaved_df = edges_df[edges_df['road_cond'] == 'unpaved']
			unpaved_length = sum(unpaved_df['length'].values.tolist())

			print ('Road paved length = ',paved_length)
			print ('Road unpaved length = ',unpaved_length)
			print ('Road paved length percentage = ',100.0*paved_length/length)
			print ('Road unpaved length percentage = ',100.0*unpaved_length/length)

			for road_class in range(1,7):
				road_class_df = edges_df[edges_df['road_class'] == road_class]
				class_length = sum(road_class_df['length'].values.tolist())

				print ('Road class {0} length = {1}'.format(road_class,class_length))
				print ('Road class {0} length percentage = {1}'.format(road_class,100.0*class_length/length))


			road_width_list = []
			unique_rd_widths = sorted(list(set(edges_df['width'].values.tolist())))
			for u_rdw in unique_rd_widths:
				unique_rd_widths_df = edges_df[edges_df['width'] == u_rdw]
				u_rdw_count = len(unique_rd_widths_df.index)
				u_rdw_length = sum(unique_rd_widths_df['length'].values.tolist())
				road_width_list.append((u_rdw,u_rdw_count,u_rdw_length,100.0*u_rdw_length/length))


			road_width_list_df = pd.DataFrame(road_width_list,columns = ['Road width','Link Count','Length (km)','Length (%)'])
			road_width_list_df.to_csv(os.path.join(data_path,'Results','network_stats','national_road_widths.csv'),index = False)

			# plt.figure(figsize=(8,4))
			# ax = edges_df['width'].hist(bins=10)

			# plt.xlabel('Road width (meters)', fontweight='bold')
			# plt.ylabel('Count', fontweight='bold')

			# plt.tight_layout()
			# output_file = os.path.join(figure_path, 'national_road_width-histogram.png')
			# plt.savefig(output_file,dpi=500)

			# plt.close()

		if modes_cols[m] == 'rail':
			route_codes = list(set(edges_df['railwaylin'].values.tolist()))
			rail_route_list = []
			for u_r in route_codes:
				unique_r_df = edges_df[edges_df['railwaylin'] == u_r]
				u_r_length = sum(unique_r_df['length'].values.tolist())
				rail_route_list.append((u_r,u_r_length,100.0*u_r_length/length))


			rail_route_list_df = pd.DataFrame(rail_route_list,columns = ['Rail route code','Length (km)','Length (%)'])
			rail_route_list_df.to_csv(os.path.join(data_path,'Results','network_stats','national_rail_routes.csv'),index = False)


	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']

	province_stats = []
	level_names = ['National','Provincial','Local','Other']
	for prn in range(len(province_list)):
		province = province_list[prn]
		province_name = province.replace(' ','').lower()

		edges_in = os.path.join(data_path,'Results','Flow_shapefiles','weighted_edges_commune_center_flows_{0}_5_tons.shp'.format(province_name))
		edges_df = gpd.read_file(edges_in)
		length = sum(edges_df['length'].values.tolist())
		for level in range(0,4):
			level_paved_df = edges_df[(edges_df['level'] == level) & (edges_df['road_cond'] == 'paved')]
			level_unpaved_df = edges_df[(edges_df['level'] == level) & (edges_df['road_cond'] == 'unpaved')]
			level_paved_length = sum(level_paved_df['length'].values.tolist())
			level_unpaved_length = sum(level_unpaved_df['length'].values.tolist())
			province_stats.append((province,level_names[level],level_paved_length,level_unpaved_length,1.0*level_paved_length/length,1.0*level_unpaved_length/length))

	province_stats_df = pd.DataFrame(province_stats,columns = ['Province name','Road level','Paved (km)','Unpaved (km)', 'Paved(%)', 'Unpaved (%)'])
	province_stats_df.to_csv(os.path.join(data_path,'Results','network_stats','province_roads_conditions.csv'),index = False)


if __name__ == '__main__':
	main()
