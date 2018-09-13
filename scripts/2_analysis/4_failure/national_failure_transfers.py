# -*- coding: utf-8 -*-
"""
Python script to assign commodity flows on the road network in Tanzania
Created on Wed March 06 2018

@author: Raghav Pant
DONT THINK THIS SCRIPT IS NEEDED
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

def swap_min_max(x,min_col,max_col):
	'''
	'''
	if x[min_col] > x[max_col]:
		return x[max_col],x[min_col]
	else:
		return x[min_col],x[max_col]

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

def network_failure_assembly(edge_failure_dataframe,mode_fail,transport_mode,gdf_edges,save_edges = True,output_path =''):
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
	# min_ind_cols = []
	# max_ind_cols = []
	# ch_min_ind_cols = []
	# ch_max_ind_cols = []
	# for ind in industry_columns:
	# 	min_ind_cols.append('min_{}'.format(ind))
	# 	max_ind_cols.append('max_{}'.format(ind))
	# 	if ind in ('rice','tons'):
	# 		ch_min_ind_cols.append('min_{}'.format(ind))
	# 		ch_max_ind_cols.append('max_{}'.format(ind))
	# 	else:
	# 		ch_min_ind_cols.append(ind)
	# 		ch_max_ind_cols.append(ind)

	# print (len(ch_min_ind_cols))
	# print (len(ch_max_ind_cols))
	# print (len(min_ind_cols))
	# print (len(max_ind_cols))
	for t in range(len(transport_mode)):
		mode_df = gdf_edges[t]
		mode_transfer_df = edge_failure_dataframe[edge_failure_dataframe['edge_id'].str.contains(transport_mode[t])]
		# print (mode_transfer_df)

		failure_columns = mode_transfer_df.columns.values.tolist()
		failure_columns = [f for f in failure_columns if f != 'edge_id']

		for fc in failure_columns:
			mode_df[fc] = 0

		print (mode_df)
		for iter_, row in mode_transfer_df.iterrows():
			# print (row[1:])
			mode_df.loc[mode_df['edge_id'] == row['edge_id'],failure_columns] = row[failure_columns].values
	

		industry_columns = list(set([f.split('min_')[1] for f in failure_columns if 'min' in f]))

		for ind in industry_columns:
			mode_df['swap'] = mode_df.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)),axis = 1)
			mode_df[['min_{}'.format(ind),'max_{}'.format(ind)]] = mode_df['swap'].apply(pd.Series)
			mode_df.drop('swap',axis=1,inplace=True)

		if save_edges == True:
			mode_df.to_file(os.path.join(output_path,'weighted_edges_failures_national_{0}_transfer_from_{1}_10_shift.shp'.format(transport_mode[t],mode_fail)))

		del mode_df, mode_transfer_df

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	types = ['min','max']
	veh_wt = [20,800,800,1200]
	usage_factors = [(0,0),(0,0),(0.2,0.25),(0.2,0.25),(0,0)]
	speeds = [(0,0),(40,60),(9,20),(9,20),(0,0)]
	multi_md_len = 3.0

	shp_output_path = os.path.join(output_path,'failure_shapefiles')

	modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Waterways','waterways')]
	modes = ['road','rail','inland']
	md_prop_file = os.path.join(data_path,'mode_properties','mode_costs.xlsx')
	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')
	G_multi_df = []
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
			G_df = network_shapefile_to_dataframe(edges_in,md_prop_file,modes[m],speeds[m][0],speeds[m][1])

		G_multi_df.append(G_df)

	'''
	Get the modal shares
	'''
	# fail_mode = ['road','rail']
	fail_mode = ['road']

	for m in range(len(fail_mode)):
		'''
		First do single edge failures
		'''
		edge_fail_ranges = []
		for t in range(len(types)):
			mode_dict = {'road':{},'rail':{},'waterways':{}}
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_national_{0}_{1}_multi_modal_options_10_shift.csv'.format(fail_mode[m],types[t]))
			flow_df_select = pd.read_csv(df_path).fillna(0)
			flow_df_select = flow_df_select[['origin','destination','new_path','{}_tons'.format(types[t])]]
			print ('Old number of paths',len(flow_df_select.index))
			flow_df_select = flow_df_select.drop_duplicates(subset=['origin','destination','new_path'], keep='first')
			# flow_df_select = flow_df_select.groupby(['new_path'])['{}_tons'.format(types[t])].sum().reset_index()
			print ('New Number of paths',len(flow_df_select.index))
			for iter_,row in flow_df_select.iterrows():
				new_path = ast.literal_eval(row['new_path'])
				if new_path:
					for p in new_path:
						if 'road' in p:
							if p in (mode_dict['road'].keys()):
								mode_dict['road'][p].append((row['origin'],row['destination'],row['{}_tons'.format(types[t])]))
							else:
								mode_dict['road'][p] = [(row['origin'],row['destination'],row['{}_tons'.format(types[t])])]
						elif 'rail' in p:
							if p in (mode_dict['rail'].keys()):
								mode_dict['rail'][p].append((row['origin'],row['destination'],row['{}_tons'.format(types[t])]))
							else:
								mode_dict['rail'][p] = [(row['origin'],row['destination'],row['{}_tons'.format(types[t])])]
						elif 'water' in p:
							if p in (mode_dict['waterways'].keys()):
								mode_dict['waterways'][p].append((row['origin'],row['destination'],row['{}_tons'.format(types[t])]))
							else:
								mode_dict['waterways'][p] = [(row['origin'],row['destination'],row['{}_tons'.format(types[t])])]
			del flow_df_select
			edge_impact = [] 
			for key,values in mode_dict.items():
				values_tup_list = []
				for a,b in values.items():
					c = list(set(b))
					v = sum([z for (x,y,z) in c])
					values_tup_list.append((a,v))
				# values_df = pd.DataFrame(values_tup_list,columns = ['edge_id','{}_tons'.format(types[t])])
				# values_df = values_df.drop_duplicates(subset=['origin','destination','edge_id'], keep='first')
				# values_df = values_df[['edge_id','{}_tons'.format(types[t])]]
				# values_df = values_df.groupby(['edge_id'])['{}_tons'.format(types[t])].sum().reset_index()
				edge_impact.append(pd.DataFrame(values_tup_list,columns = ['edge_id','{}_tons'.format(types[t])])) 

			edge_impact = pd.concat(edge_impact, axis=0, sort = 'False', ignore_index=True)

			edge_fail_ranges.append(edge_impact)
			del edge_impact

		edge_impact = edge_fail_ranges[0]
		edge_impact = pd.merge(edge_impact,edge_fail_ranges[1],how='left', on=['edge_id']).fillna(0)
		df_path = os.path.join(output_path,'failure_results','single_edge_failures_transfers_national_{0}_10_shift.csv'.format(fail_mode[m]))
		# edge_impact = pd.read_csv(df_path)
		edge_impact.to_csv(df_path,index = False)

		network_failure_assembly(edge_impact,fail_mode[m],['road','rail','water'],G_multi_df,save_edges = True,output_path =shp_output_path)


			

if __name__ == "__main__":
	main()