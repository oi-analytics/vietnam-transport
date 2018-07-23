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

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config
from scripts.transport_network_creation import province_shapefile_to_network, add_igraph_time_costs_province_roads

def igraph_scenario_edge_failures(network_graph,edge_failure_set,flow_dataframe):
	edge_fail_dictionary = []
	edge_path_index = []
	for edge in edge_failure_set: 
		edge_index = [x for x in network_graph.es if x['edge_id'] == edge]
		if edge_index:
			edge_index = edge_index[0]
		
		# print (edge_index.index)
		# print (edge_index)
		# fr_nd = edge_index['f_node']
		# t_nd = edge_index['t_node']
		# fr_id = [x for x in network_graph.vs if x['node'] == fr_nd][0]
		# t_id = [x for x in network_graph.vs if x['node'] == t_nd][0]

			network_graph.delete_edges(edge_index.index)

			edge_path_index += flow_dataframe.loc[flow_dataframe['edge_path'].str.contains(edge)].index.tolist()


	edge_path_index = list(set(edge_path_index))
	# print (edge_path_index)
	if edge_path_index:
		for e in edge_path_index:
			od_pair = ast.literal_eval(flow_dataframe.iloc[e]['od_nodes'])

			origin_node = [x for x in network_graph.vs if x['name'] == od_pair[0]]
			destination_node = [x for x in network_graph.vs if x['name'] == od_pair[-1]]

			if not origin_node or not destination_node:
				'''
				no alternative path exists
				'''
				edge_fail_dictionary.append({'edge_id':edge,'old_cost':flow_dataframe.iloc[e]['cost'],
									'old_distance':flow_dataframe.iloc[e]['distance'],'old_time':flow_dataframe.iloc[e]['time'],
									'new_cost':0,'new_distance':0,'new_time':0})
			else:
				new_route = network_graph.get_shortest_paths(origin_node[0],to = destination_node[0], weights = 'min_cost', mode = 'OUT', output='epath')[0]
				if not new_route:
					'''
					no alternative path exists
					'''
					# path_index = path_index_list[e]
					edge_fail_dictionary.append({'edge_id':edge,'old_cost':flow_dataframe.iloc[e]['cost'],
										'old_distance':flow_dataframe.iloc[e]['distance'],'old_time':flow_dataframe.iloc[e]['time'],
										'new_cost':0,'new_distance':0,'new_time':0})
				else:
					new_cost = 0
					new_distance = 0
					new_travel_time = 0
					for n in new_route:
						# new_edge_flow_dict += Counter({str(network_graph.es[n]['edge']):commodity_values[0]})
						new_cost += network_graph.es[n]['min_cost']
						new_distance += network_graph.es[n]['length']
						new_travel_time += network_graph.es[n]['min_time']
					
					edge_fail_dictionary.append({'edge_id':edge,'old_cost':flow_dataframe.iloc[e]['cost'],
										'old_distance':flow_dataframe.iloc[e]['distance'],'old_time':flow_dataframe.iloc[e]['time'],
										'new_cost':new_cost,'new_distance':new_distance,'new_time':new_travel_time})


	return edge_fail_dictionary

def main():
	data_path = load_config()['paths']['data']
	hazard_intersections_path = load_config()['paths']['output']
	provinces = ['Lao Cai','Binh Dinh','Thanh Hoa']
	bnds = [3,4,5]
	thresholds = [1,2,3,4,999]
	sector = 'Roads'

	flow_paths_data = os.path.join(hazard_intersections_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	fail_scenarios_data = os.path.join(hazard_intersections_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')

	'''
	Path OD flow disruptions
	'''
	for province in provinces:
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()

		flow_df = pd.read_excel(flow_paths_data,sheet_name = province_name)
		# print (flow_df.index.values.tolist())
		
		# pth_key_list = flow_df['path_index'].values.tolist()

		# npth_list = flow_df['od_nodes'].values.tolist()
		# npth_list = [ast.literal_eval(npaths) for npaths in npth_list]

		# epth_list = flow_df['edge_path'].values.tolist()
		# epth_list = [ast.literal_eval(epaths) for epaths in epth_list]

		# rev_list = flow_df['net_rev'].values.tolist()

		fail_df = pd.read_excel(fail_scenarios_data,sheet_name = province_name)
		'''
		Select individual edge first
		columns of failure excel are
		band_name	
		band_num	
		climate_scenario	
		commune_id	
		commune_name	
		district_id	
		district_name	
		edge_id	
		hazard_type	
		max_val	
		min_val	
		probability	
		province_id	
		province_name	
		sector	
		year	
		length 
		'''

		'''
		First do single edge failures
		'''
		edges_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))
		G = province_shapefile_to_network(edges_in)
		# print (G.es)
		G = add_igraph_time_costs_province_roads(G,0.019)

		single_ef_list = list(set(fail_df['edge_id'].values.tolist()))
		ef_list = []
		for edge in single_ef_list:
			ef_dict = igraph_scenario_edge_failures(G,[edge],flow_df)
			if ef_dict:
				ef_list += ef_dict
			
			print ('Done with region edge',edge)

		df = pd.DataFrame(ef_list)
		df_path = os.path.join(hazard_intersections_path,'failure_results','all_edge_failure_{}.csv'.format(province_name))
		df.to_csv(df_path,index = False)

if __name__ == "__main__":
	main()