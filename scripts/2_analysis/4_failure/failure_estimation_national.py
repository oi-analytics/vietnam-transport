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
from scripts.utils import *
from scripts.transport_network_creation import *

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

def igraph_scenario_edge_failures(network_df_in,edge_failure_set,flow_dataframe,vehicle_wt,min_factor,max_factor,path_criteria,cost_criteria,time_criteria):
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
		network_graph = add_igraph_generalised_costs_network(network_graph,1,vehicle_wt,min_factor,max_factor)
		nodes_name = np.asarray([x['name'] for x in network_graph.vs])

		select_flows = flow_dataframe[flow_dataframe.index.isin(edge_path_index)]
		no_access = select_flows[(~select_flows['origin'].isin(nodes_name)) | (~select_flows['destination'].isin(nodes_name))]
		if len(no_access.index) > 0:
			for iter_,value in no_access.iterrows():
				edge_dict = {'edge_id':edge,
							'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1}
				
				edge_fail_dictionary.append({'edge_id':edge,'origin':value['origin'],'destination':value['destination'],
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
				if len(po_access.loc[origin].shape) == 1:
					destinations = po_access.loc[origin,'destination']

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

						edge_fail_dictionary.append({'edge_id':edge,'origin':origin,'destination':destinations,
											'new_distance':new_dist,'new_time':new_time,'new_cost':new_gcost,'no_access': 0})
					else:
						edge_fail_dictionary.append({'edge_id':edge,'origin':origin,'destination':destinations,
											'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})
						
				else:
					destinations = po_access.loc[origin,'destination'].values.tolist()				
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

							edge_fail_dictionary.append({'edge_id':edge,'origin':origin,'destination':destinations[p],
											'new_distance':new_dist,'new_time':new_time,'new_cost':new_gcost,'no_access': 0})
						else:
							edge_fail_dictionary.append({'edge_id':edge,'origin':origin,'destination':destinations[p],
											'new_distance':0,'new_time':0,'new_cost':0,'no_access': 1})

	return edge_fail_dictionary

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	types = ['min','max']
	path_types = ['min_edge_path','max_edge_path']
	tons_types = ['min_tons','max_tons']
	dist_types = ['min_distance','max_distance']
	time_types = ['min_time','max_time']
	cost_types = ['min_gcost','max_gcost']
	vehicle_types = ['min_vehicle_nums','max_vehicle_nums']
	rice_type = ['min_rice','max_rice']
	ind_crop_cols =['sugar','wood','steel','constructi','cement','fertilizer','coal','petroluem','manufactur',
				'fishery','meat','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp'] 

	flow_paths_data = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_paths.xlsx')
	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','national_scale_hazard_intersections.xlsx')

	cols = ['origin','destination','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']

	
	'''
	Get the modal shares
	'''
	modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Waterways','waterways'),('Waterways','waterways')]
	modes_file_paths = [('Roads','national_roads'),('Railways','national_rail')]
	modes = ['road','rail','inland','coastal']
	veh_wt = [20,800,800,1200]
	usage_factors = [(0,0),(0,0),(0.2,0.25),(0.2,0.25)]
	speeds = [(0,0),(40,60),(9,20),(9,20)]
	mode_cols = ['road','rail','inland','coastal']
	new_mode_cols = ['o','d','road','rail','inland','coastal']

	md_prop_file = os.path.join(data_path,'mode_properties','mode_costs.xlsx')
	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	for m in range(len(modes_file_paths)):
		mode_data_path = os.path.join(data_path,modes_file_paths[m][0],modes_file_paths[m][1])
		for file in os.listdir(mode_data_path):
			try:
				if file.endswith(".shp") and 'edges' in file.lower().strip():
					edges_in = os.path.join(mode_data_path, file)
			except:
				return ('Network nodes and edge files necessary')
			
		if modes[m] == 'road': 
			G_df =  national_road_shapefile_to_dataframe(edges_in,rd_prop_file)
		else:
			G_df =  network_shapefile_to_dataframe(edges_in,md_prop_file,modes[m],speeds[m][0],speeds[m][1])

		# G = add_igraph_generalised_costs_network(G,1,veh_wt[m],usage_factors[m][0],usage_factors[m][1])
		# nodes_name = np.asarray([x['name'] for x in G.vs])
		# nodes_index = np.asarray([x.index for x in G.vs])
		# node_dict = dict(zip(nodes_name,nodes_index))


		fail_df = pd.read_excel(fail_scenarios_data,sheet_name = modes[m])
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
		# G = province_shapefile_to_network(edges_in,province_terrian[prn],rd_prop_file)
		# G_df = province_shapefile_to_dataframe(edges_in,province_terrian[prn],rd_prop_file)
		flow_df = pd.read_excel(flow_paths_data,sheet_name = modes[m])
		edge_fail_ranges = []
		for t in range(len(types)):
		# for t in range(1,2):
			ef_list = []
			for edge in single_ef_list:
				ef_dict = igraph_scenario_edge_failures(G_df,[edge],flow_df,veh_wt[m],usage_factors[m][0],usage_factors[m][1],path_types[t],cost_types[t],time_types[t])
				if ef_dict:
					ef_list += ef_dict
			
				print ('Done with mode {0} edge {1} type {2}'.format(modes[m],edge,types[t]))

			df = pd.DataFrame(ef_list)
			df.to_csv(os.path.join(output_path,'failure_results','single_edge_failures_all_paths_national_{0}_{1}.csv'.format(modes[m],types[t])),index = False)
			
			select_cols = ['origin','destination','o_region','d_region',dist_types[t],time_types[t],cost_types[t],vehicle_types[t]] + ind_crop_cols + [rice_type[t],tons_types[t]]
			flow_df_select = flow_df[select_cols]
			flow_df_select = pd.merge(flow_df_select,df,on = ['origin','destination'],how = 'left').fillna(0)
			flow_df_select = flow_df_select[(flow_df_select[tons_types[t]] > 0) & (flow_df_select['edge_id'] != 0)]
			
			flow_df_select['dist_diff'] = flow_df_select['new_distance'] - flow_df_select[dist_types[t]]
			flow_df_select['time_diff'] = flow_df_select['new_time'] - flow_df_select[time_types[t]]
			flow_df_select['transport_loss'] = (1 - flow_df_select['no_access'])*flow_df_select[vehicle_types[t]]*(flow_df_select['new_cost'] - flow_df_select[cost_types[t]])
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_national_{0}_{1}.csv'.format(modes[m],types[t]))
			flow_df_select.to_csv(df_path,index = False)

			select_cols = ['edge_id','o_region','d_region','no_access'] + ind_crop_cols + [rice_type[t],tons_types[t]]
			egde_impact = flow_df_select[select_cols]
			egde_impact = egde_impact[egde_impact['no_access'] == 1]
			egde_impact = egde_impact.groupby(['edge_id', 'o_region','d_region'])[ind_crop_cols + [rice_type[t],tons_types[t]]].sum().reset_index()
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_national_{0}_{1}.csv'.format(modes[m],types[t]))
			egde_impact.to_csv(df_path,index = False)
			# edge_fail_ranges.append(egde_impact)
			egde_impact = flow_df_select[select_cols+['transport_loss']]
			egde_impact = egde_impact[egde_impact['no_access'] == 0]
			egde_impact = egde_impact.groupby(['edge_id', 'o_region','d_region'])['transport_loss'].sum().reset_index()
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_transport_loss_national_{0}_{1}.csv'.format(modes[m],types[t]))
			egde_impact.to_csv(df_path,index = False)

			

if __name__ == "__main__":
	main()