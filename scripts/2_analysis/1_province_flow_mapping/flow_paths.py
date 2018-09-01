# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 15:41:39 2018

@author: elcok
"""

import geopandas as gpd
import pandas as pd
import os
import igraph as ig
import numpy as np
import sys
import subprocess
from shapely.geometry import Point
import itertools
import operator
import ast
import math

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

def network_od_paths_check(points_dataframe,node_dict,graph,vehicle_wt):
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
	save_paths = []
	cr_tons = ['min_croptons','max_croptons']
	g_costs = ['min_gcost','max_gcost']
	for iter_,row in points_dataframe.iterrows():
		try:
			od_pair = ast.literal_eval(row['od_nodes'])
			pos0_i = graph.vs[node_dict[od_pair[0]]]
			pos1_i = graph.vs[node_dict[od_pair[1]]]
			# print (od_pair,pos0_i,pos1_i)
			od_paths = []
			if pos0_i != pos1_i:
				for t in range(len(cr_tons)):
					tons = row[cr_tons[t]]
					vh_nums = math.ceil(1.0*tons/vehicle_wt)
					graph = add_igraph_generalised_costs_roads(graph,vh_nums,tons)
					path = graph.get_shortest_paths(pos0_i,pos1_i,weights=g_costs[t],output="epath")
				
					# get the path edges, path length 
					get_path = [graph.es[n]['edge_id'] for n in path][0]
					if get_path not in od_paths:
						od_paths.append(get_path)

				if len(od_paths) > 1:
					print ('different paths',od_pair)
		except:
			print(iter_)

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

def network_od_paths_assembly(points_dataframe,node_dict,graph,vehicle_wt,region_name,gdf_edges,save_edges = True,output_path ='',excel_writer =''):
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
	# save_paths = []
	# for iter_,row in points_dataframe.iterrows():
	# 	od_pair = ast.literal_eval(row['od_nodes'])
	# 	if od_pair[0] in node_dict.keys() and od_pair[1] in node_dict.keys():
	# 		pos0_i = graph.vs[node_dict[od_pair[0]]]
	# 		pos1_i = graph.vs[node_dict[od_pair[1]]]
	# 		# print (od_pair,pos0_i,pos1_i)
	# 		if pos0_i != pos1_i:
	# 			tons = row['max_croptons']
	# 			if tons > 0:
	# 				vh_nums = math.ceil(1.0*tons/vehicle_wt)
	# 			else:
	# 				vh_nums = 1.0

	# 			get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,pos0_i,pos1_i,'min_gcost','min_time')
	# 			get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,pos0_i,pos1_i,'max_gcost','max_time')
				
	# 			if get_min_path == get_max_path:
	# 				minmax_path_match = 'True'
	# 			else:
	# 				minmax_path_match = 'False'
		

	# 			save_paths.append((od_pair,get_min_path,get_max_path,row['min_netrev'],row['max_netrev'],
	# 								row['min_croptons'],row['max_croptons'],get_min_dist,get_max_dist,get_min_time,
	# 								get_max_time,vh_nums*get_min_gcost,vh_nums*get_max_gcost,minmax_path_match))				
	save_paths = []
	points_dataframe = points_dataframe.set_index('origin')
	# print (points_dataframe)
	origins = list(set(points_dataframe.index.values.tolist()))
	# print (origins)
	for origin in origins:
		try:
			destinations = points_dataframe.loc[origin,'destination'].values.tolist()

			min_croptons = points_dataframe.loc[origin,'min_croptons'].values.tolist()
			max_croptons = points_dataframe.loc[origin,'max_croptons'].values.tolist()

			min_rev = points_dataframe.loc[origin,'min_netrev'].values.tolist()
			max_rev = points_dataframe.loc[origin,'max_netrev'].values.tolist()

			min_veh_nums = points_dataframe.loc[origin,'min_vehicle_nums'].values.tolist()
			max_veh_nums = points_dataframe.loc[origin,'max_vehicle_nums'].values.tolist()

			get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
			get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,origin,destinations,'max_gcost','max_time')

			# get_min_gcost = np.multiply(np.array(min_veh_nums),np.array(get_min_gcost))
			# get_max_gcost = np.multiply(np.array(max_veh_nums),np.array(get_max_gcost))

			save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,min_rev,max_rev,min_croptons,max_croptons,
								get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost,min_veh_nums,max_veh_nums))
			print ("done with {0} in province {1}".format(origin,region_name))
		except:
			print(origin)


	# print (len(save_paths))

	# save_paths = []
	# for iter_,row in points_dataframe.iterrows():
	# 	try:
	# 		od_pair = ast.literal_eval(row['od_nodes'])
	# 		pos0_i = graph.vs[node_dict[od_pair[0]]]
	# 		pos1_i = graph.vs[node_dict[od_pair[1]]]
	# 		# print (od_pair,pos0_i,pos1_i)
	# 		if pos0_i != pos1_i:
	# 			tons = row['max_croptons']
	# 			if tons > 0:
	# 				vh_nums = math.ceil(1.0*tons/vehicle_wt)
	# 			else:
	# 				vh_nums = 1.0

	# 			get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,pos0_i,pos1_i,'min_gcost','min_time')
	# 			get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,pos0_i,pos1_i,'max_gcost','max_time')
				
	# 			if get_min_path == get_max_path:
	# 				minmax_path_match = 'True'
	# 			else:
	# 				minmax_path_match = 'False'
		

	# 			save_paths.append((od_pair,get_min_path,get_max_path,row['min_netrev'],row['max_netrev'],
	# 								row['min_croptons'],row['max_croptons'],get_min_dist,get_max_dist,get_min_time,
	# 								get_max_time,vh_nums*get_min_gcost,vh_nums*get_max_gcost,minmax_path_match))				

	# 	except:
	# 		print(iter_)

	# cols = ['od_nodes','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
	# 		'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost','minmax_path_match']
	cols = ['origin','destination','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost','min_vehicle_nums','max_vehicle_nums']
	save_paths_df = pd.DataFrame(save_paths,columns = cols)
	save_paths_df.to_excel(excel_writer,region_name + '_{}_tons'.format(int(vehicle_wt)),index = False)
	excel_writer.save()
	del save_paths_df

	# all_edges = [x['edge_id'] for x in graph.es]
	# all_edges_geom = [x['geometry'] for x in graph.es]
	
	# gdf_edges = gpd.GeoDataFrame(pd.DataFrame([all_edges,all_edges_geom]).T,crs='epsg:4326')
	# gdf_edges.columns = ['edge_id','geometry']
	
	gdf_edges['min_netrev'] = 0
	gdf_edges['max_netrev'] = 0
	gdf_edges['min_tons'] = 0
	gdf_edges['max_tons'] = 0

	for path in save_paths:
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[2]),'min_netrev'] += path[4]
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[3]),'max_netrev'] += path[5]
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[2]),'min_tons'] += path[6]
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[3]),'max_tons'] += path[7]
	
	gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_netrev','max_netrev'),axis = 1)
	gdf_edges[['min_netrev','max_netrev']] = gdf_edges['swap'].apply(pd.Series)
	gdf_edges.drop('swap',axis=1,inplace=True)

	gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_tons','max_tons'),axis = 1)
	gdf_edges[['min_tons','max_tons']] = gdf_edges['swap'].apply(pd.Series)
	gdf_edges.drop('swap',axis=1,inplace=True)

	if save_edges == True:
		# gdf_edges.to_file(os.path.join(output_path,'weighted_edges_district_center_flows_{0}_{1}_tons.shp'.format(region_name,int(vehicle_wt))))
		gdf_edges.to_file(os.path.join(output_path,'weighted_edges_commune_center_flows_{0}_{1}_tons.shp'.format(region_name,int(vehicle_wt))))


if __name__ == '__main__':

	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	truck_unit_wt = [5.0,20.0]
	truck_unit_wt = [5.0]
	# provinces to consider 
	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']
	province_terrian = ['mountain','flat','flat']

	shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles')
	# od_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_ods.xlsx')
	# flow_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	# excl_wrtr = pd.ExcelWriter(flow_output_excel)

	od_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_commune_center_flow_ods.xlsx')
	flow_output_excel = os.path.join(output_path,'flow_mapping_paths','province_roads_commune_center_flow_paths.xlsx')
	excl_wrtr = pd.ExcelWriter(flow_output_excel)

	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	

	for prn in range(len(province_list)):
	# for prn in range(0,1):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()
		
		edges_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))

		G = province_shapefile_to_network(edges_in,province_terrian[prn],rd_prop_file)
		gdf_edges = province_shapefile_to_dataframe(edges_in,province_terrian[prn],rd_prop_file)
		nodes_name = np.asarray([x['name'] for x in G.vs])
		nodes_index = np.asarray([x.index for x in G.vs])
		node_dict = dict(zip(nodes_name,nodes_index))

		all_ods = pd.read_excel(od_output_excel,sheet_name = province_name)
		all_ods = all_ods[['origin','destination','min_croptons','max_croptons','min_netrev','max_netrev']]
		for tr_wt in truck_unit_wt:
			all_ods['min_vehicle_nums'] = np.maximum(1,np.ceil(all_ods['min_croptons']/tr_wt))
			all_ods['max_vehicle_nums'] = np.maximum(1,np.ceil(all_ods['max_croptons']/tr_wt))
			# print (all_ods)
			G = add_igraph_generalised_costs_roads(G,1,tr_wt)
			network_od_paths_assembly(all_ods,node_dict,G,tr_wt,province_name,gdf_edges,save_edges = True,output_path =shp_output_path,excel_writer =excl_wrtr)
	
