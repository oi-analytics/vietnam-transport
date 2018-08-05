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
from scripts.transport_network_creation import province_shapefile_to_network, add_igraph_generalised_costs_province_roads

def igraph_scenario_edge_failures(network_in,edge_failure_set,flow_dataframe,vehicle_wt,path_criteria,rev_criteria,tons_criteria,cost_criteria,distance_criteria,time_criteria):
	network_graph = copy.deepcopy(network_in)
	edge_fail_dictionary = []
	edge_path_index = []
	for edge in edge_failure_set: 
		edge_index = [x for x in network_graph.es if x['edge_id'] == edge]
		if edge_index:
			edge_index = edge_index[0]
			network_graph.delete_edges(edge_index.index)
			edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains(edge)].index.tolist()


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
				edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
									'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
									'econ_loss':flow_dataframe.iloc[e][rev_criteria],'new_distance':0,'new_time':0,'no_acess': True})
			else:
				tons = flow_dataframe.iloc[e][tons_criteria]
				vh_nums = math.ceil(1.0*tons/vehicle_wt)
				network_graph = add_igraph_generalised_costs_province_roads(network_graph,vh_nums,tons)
				new_route = network_graph.get_shortest_paths(origin_node[0],destination_node[0], weights = cost_criteria, output='epath')[0]
				if not new_route:
					'''
					no alternative path exists
					'''
					# path_index = path_index_list[e]
					edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
									'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
									'econ_loss':flow_dataframe.iloc[e][rev_criteria],'new_distance':0,'new_time':0,'no_acess': True})
				else:
					new_dist = sum([network_graph.es[n]['length'] for n in new_route])
					new_time = sum([network_graph.es[n][time_criteria] for n in new_route])
					new_travel_cost = sum([network_graph.es[n][cost_criteria] for n in new_route])
					
					edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
										'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
										'econ_loss':new_travel_cost - flow_dataframe.iloc[e][cost_criteria],'new_distance':new_dist,'new_time':new_time,'no_acess': False})


	return edge_fail_dictionary

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	truck_unit_wt = 20.0
	# provinces to consider 
	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']
	province_terrian = ['mountain','flat','flat']

	types = ['min','max']
	path_types = ['min_edge_path','max_edge_path']
	rev_types = ['min_netrev','max_netrev']
	tons_types = ['min_croptons','max_croptons']
	dist_types = ['min_distance','max_distance']
	time_types = ['min_time','max_time']
	cost_types = ['min_gcost','max_gcost']

	# shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles')
	# flow_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	# excl_wrtr = pd.ExcelWriter(flow_output_excel)

	flow_paths_data = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')

	rd_prop_file = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')

	'''
	Path OD flow disruptions
	'''
	# for prn in range(len(province_list)):
	for prn in range(0,1):
		province_ods_df = []
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()

		cols = ['od_nodes','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
				'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']
		flow_df = pd.read_excel(flow_paths_data,sheet_name = province_name)

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
		G = province_shapefile_to_network(edges_in,province_terrian[prn],rd_prop_file)

		single_ef_list = list(set(fail_df['edge_id'].values.tolist()))
		edge_fail_ranges = []
		for t in range(len(types)):
			ef_list = []
			for edge in single_ef_list:
				ef_dict = igraph_scenario_edge_failures(G,[edge],flow_df,truck_unit_wt,path_types[t],rev_types[t],tons_types[t],cost_types[t],dist_types[t],time_types[t])
				if ef_dict:
					ef_list += ef_dict
			
				print ('Done with province {0} edge {1} type {2}'.format(province_name,edge,types[t]))

			df = pd.DataFrame(ef_list)
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_{0}_{1}.csv'.format(province_name,types[t]))
			df.to_csv(df_path,index = False)

			egde_impact = df[['edge_id','econ_value','tons','econ_loss']]
			egde_impact = egde_impact.groupby(['edge_id'])['econ_value','tons','econ_loss'].sum().reset_index()
			egde_impact.rename(columns={'econ_value': '{}_econ_value'.format(types[t]), 'tons': '{}_tons'.format(types[t]),'econ_loss':'{}_econ_loss'.format(types[t])}, inplace=True)
			edge_fail_ranges.append(egde_impact)

		egde_impact = edge_fail_ranges[0]
		egde_impact = pd.merge(egde_impact,edge_fail_ranges[1],how='left', on=['edge_id'])
		df_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_{0}.csv'.format(province_name))
		egde_impact.to_csv(df_path,index = False)

		
			

if __name__ == "__main__":
	main()