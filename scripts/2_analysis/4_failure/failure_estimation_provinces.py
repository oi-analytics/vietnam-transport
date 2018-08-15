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

# def igraph_scenario_edge_failures(network_df_in,edge_failure_set,flow_dataframe,vehicle_wt,path_criteria,rev_criteria,tons_criteria,cost_criteria,distance_criteria,time_criteria):
# 	network_graph_df = copy.deepcopy(network_df_in)
# 	edge_fail_dictionary = []
# 	edge_path_index = []
# 	for edge in edge_failure_set:
# 		network_graph_df = network_graph_df[network_graph_df.edge_id != edge]
# 		edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains(edge)].index.tolist() 


# 	edge_path_index = list(set(edge_path_index))
# 	# print (edge_path_index)
# 	if edge_path_index:
# 		network_graph = ig.Graph.TupleList(network_graph_df.itertuples(index=False), edge_attrs=list(network_graph_df.columns)[2:])
# 		# only keep connected network
# 		network_graph = network_graph.clusters().giant()

# 		for e in edge_path_index:
# 			od_pair = ast.literal_eval(flow_dataframe.iloc[e]['od_nodes'])

# 			origin_node = [x for x in network_graph.vs if x['name'] == od_pair[0]]
# 			destination_node = [x for x in network_graph.vs if x['name'] == od_pair[-1]]

# 			if not origin_node or not destination_node:
# 				'''
# 				no alternative path exists
# 				'''
# 				edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
# 									'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
# 									'econ_loss':flow_dataframe.iloc[e][rev_criteria],'new_distance':0,'new_time':0,'no_access': True})
# 			else:
# 				tons = flow_dataframe.iloc[e][tons_criteria]
# 				vh_nums = math.ceil(1.0*tons/vehicle_wt)
# 				network_graph = add_igraph_generalised_costs_province_roads(network_graph,vh_nums,tons)
# 				new_route = network_graph.get_shortest_paths(origin_node[0],destination_node[0], weights = cost_criteria, output='epath')[0]
# 				if not new_route:
# 					'''
# 					no alternative path exists
# 					'''
# 					# path_index = path_index_list[e]
# 					edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
# 									'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
# 									'econ_loss':flow_dataframe.iloc[e][rev_criteria],'new_distance':0,'new_time':0,'no_access': True})
# 				else:
# 					new_dist = 0
# 					new_time = 0
# 					new_travel_cost = 0
# 					for n in new_route:
# 						new_dist += network_graph.es[n]['length']
# 						new_time += network_graph.es[n][time_criteria]
# 						new_travel_cost += network_graph.es[n][cost_criteria]
					
# 					edge_fail_dictionary.append({'edge_id':edge,'econ_value':flow_dataframe.iloc[e][rev_criteria],'tons':flow_dataframe.iloc[e][tons_criteria],
# 										'old_distance':flow_dataframe.iloc[e][distance_criteria],'old_time':flow_dataframe.iloc[e][time_criteria],
# 										'econ_loss':new_travel_cost - flow_dataframe.iloc[e][cost_criteria],'new_distance':new_dist,'new_time':new_time,'no_access': False})


# 	return edge_fail_dictionary

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
		network_graph = add_igraph_generalised_costs_province_roads(network_graph,1,vehicle_wt)
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
				if len(po_access.loc[sc].index) == 1:
					destinations = po_access.loc[origin,'destination']
					croptons = po_access.loc[origin,tons_criteria]
					rev = po_access.loc[origin,rev_criteria]
					veh_nums = po_access.loc[origin,vehicle_criteria]
					dist = po_access.loc[origin,distance_criteria]
					time = po_access.loc[origin,time_criteria]
					gcost = po_access.loc[origin,cost_criteria]

					# compute min cost paths and values
					paths = network_graph.get_shortest_paths(origin,destinations,weights=cost_criteria,output="epath")[0]
					if len(paths) > 0:
						new_dist = 0
						new_time = 0
						new_gcost = 0
						for n in paths:
							new_dist += network_graph.es[n]['length']
							new_time += network_graph.es[n][time_criteria]
							new_gcost += network_graph.es[n][cost_criteria]

						edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev,'tons':croptons,'vehicle_nums':veh_nums,
										'old_distance':dist,'old_time':time,'old_cost':gcost,
										'new_distance':new_dist,'new_time':new_time,'new_cost':new_gcost,'no_access': 0})
					else:
						edge_fail_dictionary.append({'edge_id':edge,'econ_value':rev,'tons':croptons,'vehicle_nums':veh_nums,
									'old_distance':dist,'old_time':time,'old_cost':gcost,
									'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})
				else:
					destinations = po_access.loc[origin,'destination'].values.tolist()
					croptons = po_access.loc[origin,tons_criteria].values.tolist()
					rev = po_access.loc[origin,rev_criteria].values.tolist()
					veh_nums = po_access.loc[origin,vehicle_criteria].values.tolist()
					dist = po_access.loc[origin,distance_criteria].values.tolist()
					time = po_access.loc[origin,time_criteria].values.tolist()
					gcost = po_access.loc[origin,cost_criteria].values.tolist()
				
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

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	truck_unit_wt = [5.0,20.0]
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

	flow_paths_data = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')

	rd_prop_file = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')

	cols = ['origin','destination','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']

	'''
	Path OD flow disruptions
	'''
	for prn in range(len(province_list)):
	# for prn in range(1,3):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()

		fail_df = pd.read_excel(fail_scenarios_data,sheet_name = province_name)
		single_ef_list = list(set(fail_df['edge_id'].values.tolist()))
		print (len(single_ef_list))
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
			flow_df = pd.read_excel(flow_paths_data,sheet_name = province_name + '_{}_tons'.format(int(tr_wt)))
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
				df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_{0}_{1}_{2}_tons.csv'.format(province_name,types[t],int(tr_wt)))
				df.to_csv(df_path,index = False)

				df['econ_loss'] = df['no_access']*df['econ_value'] + (1 - df['no_access'])*df['vehicle_nums']*(df['new_cost'] - df['old_cost'])
				egde_impact = df[['edge_id','econ_value','tons','econ_loss']]
				egde_impact = egde_impact.groupby(['edge_id'])['econ_value','tons','econ_loss'].sum().reset_index()
				egde_impact.rename(columns={'econ_value': '{}_econ_value'.format(types[t]), 'tons': '{}_tons'.format(types[t]),'econ_loss':'{}_econ_loss'.format(types[t])}, inplace=True)
				edge_fail_ranges.append(egde_impact)

			egde_impact = edge_fail_ranges[0]
			egde_impact = pd.merge(egde_impact,edge_fail_ranges[1],how='left', on=['edge_id'])
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_{0}_{1}_tons.csv'.format(province_name,int(tr_wt)))
			egde_impact.to_csv(df_path,index = False)
		
			

if __name__ == "__main__":
	main()