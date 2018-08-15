# -*- coding: utf-8 -*-
"""
Python script to assign commodity flows on the road network in Tanzania
Created on Wed March 06 2018

@author: Raghav Pant
"""

import pandas as pd
import os
import psycopg2
import networkx as nx
import csv
import itertools
import operator
import ast
from sqlalchemy import create_engine
import numpy as np
import igraph as ig 
import copy
from collections import Counter
import sys
import math
import copy 
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config
from scripts.transport_network_creation import province_shapefile_to_network, add_igraph_generalised_costs_province_roads, province_shapefile_to_dataframe

def assign_failed_edge_attributes(x,edge_properties):
	'''
	'''
	edge_index = [e for e in range(len(edge_properties)) if edge_properties[e][0] == x.edge_id][0]
	return edge_properties[edge_index][1], edge_properties[edge_index][2]


def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	# provinces to consider 
	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']
	province_terrian = ['mountain','flat','flat']

	# shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles')
	# flow_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	# excl_wrtr = pd.ExcelWriter(flow_output_excel)

	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')
	rd_prop_file = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')

	'''
	Path OD flow disruptions
	'''
	# for prn in range(len(province_list)):
	for prn in range(0,1):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()

		all_edge_fail_scenarios = pd.read_excel(fail_scenarios_data,sheet_name = province_name)
		all_edges = list(set(all_edge_fail_scenarios['edge_id'].values.tolist()))
		all_edge_fail_scenarios['road_cond'] = 'unknown'
		all_edge_fail_scenarios['asset_type'] = 'unknown'
		all_edge_fail_scenarios['width'] = 0
		# print ('done with file reading')

		'''
		First do single edge failures
		'''
		edges_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))
		edges = province_shapefile_to_dataframe(edges_in,province_terrian[prn],rd_prop_file)
		edge_attr = list(zip(edges['edge_id'].values.tolist(),edges['road_cond'].values.tolist(),edges['asset_type'].values.tolist(),edges['width'].values.tolist()))
		# print (edge_attr)
		
		edge_attr = [e for e in edge_attr if e[0] in all_edges]
		for e in edge_attr:
			all_edge_fail_scenarios.loc[all_edge_fail_scenarios['edge_id'] == e[0], 'road_cond'] = e[1]
			all_edge_fail_scenarios.loc[all_edge_fail_scenarios['edge_id'] == e[0], 'asset_type'] = e[2]
			all_edge_fail_scenarios.loc[all_edge_fail_scenarios['edge_id'] == e[0], 'width'] = e[3]

		


		# all_edge_fail_scenarios['attributes'] = all_edge_fail_scenarios.apply(lambda x: assign_failed_edge_attributes(x,edge_attr),axis = 1)
		# all_edge_fail_scenarios['asset_type','width'] = all_edge_fail_scenarios['attributes'].apply(pd.Series)
		# all_edge_fail_scenarios.drop('attributes',axis=1,inplace=True)


		single_edge_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_{0}.csv'.format(province_name))
		edge_impacts = pd.read_csv(single_edge_path)
		edge_impacts_attr = list(zip(edge_impacts['edge_id'].values.tolist(),edge_impacts['min_econ_loss'].values.tolist(),edge_impacts['max_econ_loss'].values.tolist())) 

		all_edge_fail_scenarios['min_econ_loss'] = 0
		all_edge_fail_scenarios['max_econ_loss'] = 0

		edge_impacts_attr = [e for e in edge_impacts_attr if e[0] in all_edges]
		for e in edge_impacts_attr:
			all_edge_fail_scenarios.loc[all_edge_fail_scenarios['edge_id'] == e[0], 'min_econ_loss'] = e[1]
			all_edge_fail_scenarios.loc[all_edge_fail_scenarios['edge_id'] == e[0], 'max_econ_loss'] = e[2]

		# all_edge_fail_scenarios['attributes'] = all_edge_fail_scenarios.apply(lambda x: assign_failed_edge_attributes(x,edge_impacts_attr),axis = 1)
		# all_edge_fail_scenarios['min_econ_loss','max_econ_loss'] = all_edge_fail_scenarios['attributes'].apply(pd.Series)
		# all_edge_fail_scenarios.drop('attributes',axis=1,inplace=True)

		df_path = os.path.join(output_path,'hazard_scenarios','roads_hazard_intersections_{}.csv'.format(province_name))
		all_edge_fail_scenarios.to_csv(df_path,index = False)

		
			

if __name__ == "__main__":
	main()