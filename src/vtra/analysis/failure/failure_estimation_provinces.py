# -*- coding: utf-8 -*-
"""
Python script to assign commodity flows on the road network
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

from vtra.utils import *
from vtra.transport_network_creation import *

def igraph_scenario_edge_failures_changing_tonnages(network_df_in,edge_failure_set,flow_dataframe,vehicle_wt,path_criteria,rev_criteria,tons_criteria,cost_criteria,distance_criteria,time_criteria):
	network_graph_df = copy.deepcopy(network_df_in)
	edge_fail_dictionary = []
	edge_path_index = []
	for edge in edge_failure_set:
		network_graph_df = network_graph_df[network_graph_df.edge_id != edge]
		edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains(edge)].index.tolist()


	edge_path_index = list(set(edge_path_index))
	# print (edge_path_index)
	if edge_path_index:
		network_graph = ig.Graph.TupleList(network_graph_df.itertuples(index=False), edge_attrs=list(network_graph_df.columns)[2:])
		# only keep connected network
		network_graph = network_graph.clusters().giant()

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
									'econ_loss':flow_dataframe.iloc[e][rev_criteria],'new_distance':0,'new_time':0,'no_access': True})
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
									'econ_loss':flow_dataframe.iloc[e][rev_criteria],'new_distance':0,'new_time':0,'no_access': True})
				else:
					new_dist = 0
					new_time = 0
					new_travel_cost = 0
					for n in new_route:
						new_dist += network_graph.es[n]['length']
						new_time += network_graph.es[n][time_criteria]
						new_travel_cost += network_graph.es[n][cost_criteria]

					edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
										'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
										'econ_loss':new_travel_cost - flow_dataframe.iloc[e][cost_criteria],'new_distance':new_dist,'new_time':new_time,'no_access': False})


	return edge_fail_dictionary

def network_od_path_estimations(graph,source,target,cost_criteria,time_criteria):
	# compute min cost paths and values
	paths = graph.get_shortest_paths(source,target,weights=cost_criteria,output="epath")

	edge_path_list = []
	path_dist_list = []
	path_time_list = []
	path_gcost_list = []

	for path in paths:
		edge_path = []
		path_dist = 0
		path_time = 0
		path_gcost = 0
		if path:
			for n in path:
				edge_path.append(graph.es[n]['edge_id'])
				path_dist += graph.es[n]['length']
				path_time += graph.es[n][time_criteria]
				path_gcost += graph.es[n][cost_criteria]

		edge_path_list.append(edge_path)
		path_dist_list.append(path_dist)
		path_time_list.append(path_time)
		path_gcost_list.append(path_gcost)


	return edge_path_list, path_dist_list,path_time_list,path_gcost_list

def igraph_scenario_edge_failures(network_df_in,edge_failure_set,flow_dataframe,vehicle_wt,path_criteria,rev_criteria,tons_criteria,cost_criteria,distance_criteria,time_criteria,vehicle_criteria):
	network_graph_df = copy.deepcopy(network_df_in)
	edge_fail_dictionary = []
	edge_path_index = []
	for edge in edge_failure_set:
		network_graph_df = network_graph_df[network_graph_df.edge_id != edge]
		edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains("'{}'".format(edge))].index.tolist()


	edge_path_index = list(set(edge_path_index))
	# print (edge_path_index)
	if edge_path_index:
		network_graph = ig.Graph.TupleList(network_graph_df.itertuples(index=False), edge_attrs=list(network_graph_df.columns)[2:])
		# only keep connected network
		network_graph = network_graph.clusters().giant()
		network_graph = add_igraph_generalised_costs_roads(network_graph,1,vehicle_wt)
		nodes_name = np.asarray([x['name'] for x in network_graph.vs])

		select_flows = flow_dataframe[flow_dataframe.index.isin(edge_path_index)]
		no_access = select_flows[(~select_flows['origin'].isin(nodes_name)) | (~select_flows['destination'].isin(nodes_name))]
		if len(no_access.index) > 0:
			for iter_,value in no_access.iterrows():
				edge_fail_dictionary.append({'edge_id':edge,'econ_value':value[rev_criteria],'tons':value[tons_criteria],
									'vehicle_nums':value[vehicle_criteria],
									'old_distance':value[distance_criteria],'old_time':value[time_criteria],'old_cost':value[cost_criteria],
									'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})


		po_access = select_flows[(select_flows['origin'].isin(nodes_name)) & (select_flows['destination'].isin(nodes_name))]
		if len(po_access.index) > 0:
			po_access = po_access.set_index('origin')
			origins = list(set(po_access.index.values.tolist()))
			# print (origins)
			# print (po_access)
			for origin in origins:
				# destinations = po_access.loc[origin,'destination']
				# print (destinations,type(destinations))
				# if isinstance(destinations,str):
				# if len(po_access.loc[origin].shape) == 1:
				# 	destinations = po_access.loc[origin,'destination']
				# 	croptons = po_access.loc[origin,tons_criteria]
				# 	rev = po_access.loc[origin,rev_criteria]
				# 	veh_nums = po_access.loc[origin,vehicle_criteria]
				# 	dist = po_access.loc[origin,distance_criteria]
				# 	time = po_access.loc[origin,time_criteria]
				# 	gcost = po_access.loc[origin,cost_criteria]

				# 	# compute min cost paths and values
				# 	paths = network_graph.get_shortest_paths(origin,destinations,weights=cost_criteria,output="epath")[0]
				# 	if len(paths) > 0:
				# 		new_dist = 0
				# 		new_time = 0
				# 		new_gcost = 0
				# 		for n in paths:
				# 			new_dist += network_graph.es[n]['length']
				# 			new_time += network_graph.es[n][time_criteria]
				# 			new_gcost += network_graph.es[n][cost_criteria]

				# 		edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev,'tons':croptons,'vehicle_nums':veh_nums,
				# 						'old_distance':dist,'old_time':time,'old_cost':gcost,
				# 						'new_distance':new_dist,'new_time':new_time,'new_cost':new_gcost,'no_access': 0})
				# 	else:
				# 		edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev,'tons':croptons,'vehicle_nums':veh_nums,
				# 					'old_distance':dist,'old_time':time,'old_cost':gcost,
				# 					'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})
				# else:
				# 	destinations = po_access.loc[origin,'destination'].values.tolist()
				# 	croptons = po_access.loc[origin,tons_criteria].values.tolist()
				# 	rev = po_access.loc[origin,rev_criteria].values.tolist()
				# 	veh_nums = po_access.loc[origin,vehicle_criteria].values.tolist()
				# 	dist = po_access.loc[origin,distance_criteria].values.tolist()
				# 	time = po_access.loc[origin,time_criteria].values.tolist()
				# 	gcost = po_access.loc[origin,cost_criteria].values.tolist()

				# 	# compute min cost paths and values
				# 	paths = network_graph.get_shortest_paths(origin,destinations,weights=cost_criteria,output="epath")
				# 	for p in range(len(paths)):
				# 		if len(paths[p]) > 0:
				# 			new_dist = 0
				# 			new_time = 0
				# 			new_gcost = 0
				# 			for n in paths[p]:
				# 				new_dist += network_graph.es[n]['length']
				# 				new_time += network_graph.es[n][time_criteria]
				# 				new_gcost += network_graph.es[n][cost_criteria]

				# 			edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev[p],'tons':croptons[p],'vehicle_nums':veh_nums[p],
				# 							'old_distance':dist[p],'old_time':time[p],'old_cost':gcost[p],
				# 							'new_distance':new_dist,'new_time':new_time,'new_cost':new_gcost,'no_access': 0})
				# 		else:
				# 			edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev[p],'tons':croptons[p],'vehicle_nums':veh_nums[p],
				# 						'old_distance':dist[p],'old_time':time[p],'old_cost':gcost[p],
				# 						'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})

				destinations = po_access.loc[[origin],'destination'].values.tolist()
				croptons = po_access.loc[[origin],tons_criteria].values.tolist()
				rev = po_access.loc[[origin],rev_criteria].values.tolist()
				veh_nums = po_access.loc[[origin],vehicle_criteria].values.tolist()
				dist = po_access.loc[[origin],distance_criteria].values.tolist()
				time = po_access.loc[[origin],time_criteria].values.tolist()
				gcost = po_access.loc[[origin],cost_criteria].values.tolist()

				# compute min cost paths and values
				paths = network_graph.get_shortest_paths(origin,destinations,weights=cost_criteria,output="epath")
				for p in range(len(paths)):
					if len(paths[p]) > 0:
						new_dist = 0
						new_time = 0
						new_gcost = 0
						for n in paths[p]:
							new_dist += network_graph.es[n]['length']
							new_time += network_graph.es[n][time_criteria]
							new_gcost += network_graph.es[n][cost_criteria]
						edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev[p],'tons':croptons[p],'vehicle_nums':veh_nums[p],
										'old_distance':dist[p],'old_time':time[p],'old_cost':gcost[p],
										'new_distance':new_dist,'new_time':new_time,'new_cost':new_gcost,'no_access': 0})
					else:
						edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev[p],'tons':croptons[p],'vehicle_nums':veh_nums[p],
									'old_distance':dist[p],'old_time':time[p],'old_cost':gcost[p],
									'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})

	return edge_fail_dictionary

def network_failure_assembly(edge_failure_dataframe,region_name,vehicle_wt,gdf_edges,save_edges = True,output_path =''):
	"""
	Assign net revenue to roads assets in Vietnam

	Inputs are:
	start_points - GeoDataFrame of start points for shortest path analysis.
	end_points - GeoDataFrame of potential end points for shorest path analysis.
	G - iGraph network of the province.
	save_edges -

	Outputs are:
	Shapefile with all edges and the total net reveneu transferred along each edge
	GeoDataFrame of total net revenue transferred along each edge
	"""

	failure_columns = edge_failure_dataframe.columns.values.tolist()
	failure_columns = [f for f in failure_columns if f != 'edge_id']
	# print (failure_columns)
	for fc in failure_columns:
		gdf_edges[fc] = 0

	for iter_, row in edge_failure_dataframe.iterrows():
		# print (row[1:])
		gdf_edges.loc[gdf_edges['edge_id'] == row['edge_id'],failure_columns] = row[failure_columns].values

	if save_edges == True:
		# gdf_edges.to_file(os.path.join(output_path,'weighted_edges_district_center_flows_{0}_{1}_tons.shp'.format(region_name,int(vehicle_wt))))
		gdf_edges.to_file(os.path.join(output_path,'weighted_edges_commune_center_failures_{0}_{1}_tons.shp'.format(region_name,int(vehicle_wt))))

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	# truck_unit_wt = [5.0,20.0]
	truck_unit_wt = [5.0]
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
	vechicle_types = ['min_vehicle_nums','max_vehicle_nums']

	shp_output_path = os.path.join(output_path,'failure_shapefiles')
	# flow_paths_data = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	flow_paths_data = os.path.join(output_path,'flow_mapping_paths','province_roads_commune_center_flow_paths.xlsx')
	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')

	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	cols = ['origin','destination','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']

	'''
	Path OD flow disruptions
	'''
	for prn in range(len(province_list)):
	# for prn in range(0,1):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()

		fail_df = pd.read_excel(fail_scenarios_data,sheet_name = province_name)
		single_ef_list = list(set(fail_df['edge_id'].values.tolist()))
		print ('scenarios in province {0} are {1}'.format(province_name,len(single_ef_list)))
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
		# G = province_shapefile_to_network(edges_in,province_terrian[prn],rd_prop_file)
		G_df = province_shapefile_to_dataframe(edges_in,province_terrian[prn],rd_prop_file)
		for tr_wt in truck_unit_wt:
			# flow_df = pd.read_excel(flow_paths_data,sheet_name = province_name + '_{}_tons'.format(int(tr_wt)))
			edge_fail_ranges = []
			for t in range(len(types)):
			# for t in range(1,2):
				ef_list = []
				for edge in single_ef_list:
					ef_dict = igraph_scenario_edge_failures(G_df,[edge],flow_df,tr_wt,path_types[t],rev_types[t],tons_types[t],cost_types[t],dist_types[t],time_types[t],vechicle_types[t])
					if ef_dict:
						ef_list += ef_dict

					print ('Done with province {0} edge {1} type {2}'.format(province_name,edge,types[t]))

				df = pd.DataFrame(ef_list)
				df_path = os.path.join(output_path,'failure_results','single_edge_failures_commune_access_all_path_impacts_{0}_{1}_{2}_tons.csv'.format(province_name,types[t],int(tr_wt)))
				# df.to_csv(df_path,index = False)
				# df = pd.read_csv(df_path)

				df['econ_loss'] = df['no_access']*df['econ_value'] + (1 - df['no_access'])*df['vehicle_nums']*(df['new_cost'] - df['old_cost'])
				df['tons_loss'] = df['no_access']*df['tons']
				edge_impact = df[['edge_id','tons_loss','econ_loss']]
				edge_impact = edge_impact.groupby(['edge_id'])['tons_loss','econ_loss'].sum().reset_index()
				edge_impact.rename(columns={'tons_loss': '{}_tons_loss'.format(types[t]),'econ_loss':'{}_econ_loss'.format(types[t])}, inplace=True)
				edge_fail_ranges.append(edge_impact)

			edge_impact = edge_fail_ranges[0]
			edge_impact = pd.merge(edge_impact,edge_fail_ranges[1],how='left', on=['edge_id']).fillna(0)
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_commune_access_totals_{0}_{1}_tons.csv'.format(province_name,int(tr_wt)))
			edge_impact.to_csv(df_path,index = False)

			# edge_impact = pd.read_csv(df_path).fillna(0)
			network_failure_assembly(edge_impact,province_name,tr_wt,G_df,save_edges = True,output_path =shp_output_path)



if __name__ == "__main__":
	main()
