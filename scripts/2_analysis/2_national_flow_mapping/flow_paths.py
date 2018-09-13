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

def network_od_path_estimations_changing_tonnages(graph,source,target,tonnage,vehicle_wt,utilization_factors,cost_criteria,time_criteria):
	# compute min cost paths and values

	graph = add_igraph_generalised_costs_network(graph,np.ceil(tonnage/vehicle_wt),tonnage,utilization_factors[0],utilization_factors[1])
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


def network_od_paths_assembly_changing_tonnages(points_dataframe,node_dict,graph,vehicle_wt,utilization_factors,transport_mode,industry_columns,gdf_edges,save_edges = True,output_path ='',excel_writer =''):
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
	for iter_,row in points_dataframe.iterrows():
		# origin = row['origin']
		# destinations = [row['destination']]
		# tons = row['min_tons']
		# get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations_changing_tonnages(graph,origin,destinations,tons,vehicle_wt,utilization_factors,'min_gcost','min_time')
		# tons = row['max_tons']
		# get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations_changing_tonnages(graph,origin,destinations,tons,vehicle_wt,utilization_factors,'max_gcost','max_time')
		# save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,
		# 					get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost))
		# print ("done with {0} in network {1}".format(origin,transport_mode))
		try: 
			origin = row['origin']
			destinations = [row['destination']]
			tons = 0.9*row['min_tons']
			get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations_changing_tonnages(graph,origin,destinations,tons,vehicle_wt,utilization_factors,'min_gcost','min_time')
			tons = 0.9*row['max_tons']
			get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations_changing_tonnages(graph,origin,destinations,tons,vehicle_wt,utilization_factors,'max_gcost','max_time')
			save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,
								get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost))
			print ("done with {0} in network {1}".format(origin,transport_mode))
		except:
			print(origin)	

	cols = ['origin','destination','min_edge_path','max_edge_path',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']
	save_paths_df = pd.DataFrame(save_paths,columns = cols)

	points_dataframe = points_dataframe.reset_index()
	save_paths_df = pd.merge(save_paths_df,points_dataframe,how='left', on=['origin','destination']).fillna(0)
	
	save_paths_df = save_paths_df[(save_paths_df['max_tons'] > 0) & (save_paths_df['origin'] != 0)]
	if transport_mode != 'air':
		save_paths_df['min_vehicle_nums'] = np.maximum(1,np.ceil(save_paths_df['min_tons']/vehicle_wt))
		save_paths_df['max_vehicle_nums'] = np.maximum(1,np.ceil(save_paths_df['max_tons']/vehicle_wt))

	save_paths_df.to_excel(excel_writer,transport_mode,index = False)
	excel_writer.save()
	del save_paths
	del save_paths_df
	
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


	# for i in range(len(min_ind_cols)):
	# 	gdf_edges[min_ind_cols[i]] = 0
	# 	gdf_edges[max_ind_cols[i]] = 0

	# for iter_,path in save_paths_df.iterrows():
	# 	min_path = path['min_edge_path']
	# 	max_path = path['max_edge_path']
		
	# 	gdf_edges.loc[gdf_edges['edge_id'].isin(min_path),min_ind_cols] += path[ch_min_ind_cols].values
	# 	gdf_edges.loc[gdf_edges['edge_id'].isin(max_path),max_ind_cols] += path[ch_max_ind_cols].values
	
	# # gdf_edges[min_ind_cols] = gdf_edges['min_vals'].apply(pd.Series)
	# # gdf_edges[max_ind_cols] = gdf_edges['max_vals'].apply(pd.Series)
	# # gdf_edges.drop('min_vals',axis=1,inplace=True)
	# # gdf_edges.drop('max_vals',axis=1,inplace=True)


	# for ind in industry_columns:
	# 	gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)),axis = 1)
	# 	gdf_edges[['min_{}'.format(ind),'max_{}'.format(ind)]] = gdf_edges['swap'].apply(pd.Series)
	# 	gdf_edges.drop('swap',axis=1,inplace=True)

	# if save_edges == True:
	# 	gdf_edges.to_file(os.path.join(output_path,'weighted_edges_flows_national_{0}_2.shp'.format(transport_mode)))

	# del gdf_edges, save_paths_df

def network_od_paths_assembly(points_dataframe,node_dict,graph,vehicle_wt,transport_mode,industry_columns,gdf_edges,save_edges = True,output_path ='',excel_writer =''):
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
	points_dataframe = points_dataframe.set_index('origin')
	# print (points_dataframe)
	origins = list(set(points_dataframe.index.values.tolist()))
	# print (origins)
	for origin in origins:
		# origin_df =  points_dataframe.loc[origin]
		# # print (list(origin_df['destination']),origin_df.shape)
		# if len(origin_df.shape) == 1:
		# 	# print (origin_df['destination'])
		# 	destinations = [origin_df['destination']]

		# 	min_tons = [origin_df['min_tons']]
		# 	max_tons = [origin_df['max_tons']]

		# else:
		# 	# print (origin_df)
		# 	destinations = origin_df.loc[origin,'destination'].values.tolist()

		# 	min_tons = origin_df.loc[origin,'min_tons'].values.tolist()
		# 	max_tons = origin_df.loc[origin,'max_tons'].values.tolist()


		# get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
		# get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,origin,destinations,'max_gcost','max_time')

		# # get_min_gcost = np.multiply(np.array(min_veh_nums),np.array(get_min_gcost))
		# # get_max_gcost = np.multiply(np.array(max_veh_nums),np.array(get_max_gcost))

		# save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,min_tons,max_tons,
		# 						get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost))
		# print ("done with {0} in network {1}".format(origin,transport_mode))
		try:
			origin_df =  points_dataframe.loc[origin]
			# print (list(origin_df['destination']),origin_df.shape)
			if len(origin_df.shape) == 1:
				# print (origin_df['destination'])
				destinations = [origin_df['destination']]

				# min_tons = [origin_df['min_tons']]
				# max_tons = [origin_df['max_tons']]

				# min_veh_nums = [origin_df['min_vehicle_nums']]
				# max_veh_nums = [origin_df['max_vehicle_nums']]
				# get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
				# print (get_min_path)
			else:
				# print (origin_df)
				destinations = origin_df['destination'].values.tolist()

				# min_tons = origin_df['min_tons'].values.tolist()
				# max_tons = origin_df['max_tons'].values.tolist()

				# min_veh_nums = origin_df['min_vehicle_nums'].values.tolist()
				# max_veh_nums = origin_df['max_vehicle_nums'].values.tolist()

			get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
			get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,origin,destinations,'max_gcost','max_time')

			# get_min_gcost = np.multiply(np.array(min_veh_nums),np.array(get_min_gcost))
			# get_max_gcost = np.multiply(np.array(max_veh_nums),np.array(get_max_gcost))

			save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,
								get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost))
			print ("done with {0} in network {1}".format(origin,transport_mode))
		except:
			print(origin)


	cols = ['origin','destination','min_edge_path','max_edge_path',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']
	save_paths_df = pd.DataFrame(save_paths,columns = cols)

	points_dataframe = points_dataframe.reset_index()
	save_paths_df = pd.merge(save_paths_df,points_dataframe,how='left', on=['origin','destination']).fillna(0)
	
	save_paths_df = save_paths_df[(save_paths_df['max_tons'] > 0) & (save_paths_df['origin'] != 0)]
	if transport_mode != 'air':
		save_paths_df['min_vehicle_nums'] = np.maximum(1,np.ceil(save_paths_df['min_tons']/vehicle_wt))
		save_paths_df['max_vehicle_nums'] = np.maximum(1,np.ceil(save_paths_df['max_tons']/vehicle_wt))

	save_paths_df.to_excel(excel_writer,transport_mode,index = False)
	excel_writer.save()
	del save_paths
	
	min_ind_cols = []
	max_ind_cols = []
	ch_min_ind_cols = []
	ch_max_ind_cols = []
	for ind in industry_columns:
		min_ind_cols.append('min_{}'.format(ind))
		max_ind_cols.append('max_{}'.format(ind))
		if ind in ('rice','tons'):
			ch_min_ind_cols.append('min_{}'.format(ind))
			ch_max_ind_cols.append('max_{}'.format(ind))
		else:
			ch_min_ind_cols.append(ind)
			ch_max_ind_cols.append(ind)

	print (len(ch_min_ind_cols))
	print (len(ch_max_ind_cols))
	print (len(min_ind_cols))
	print (len(max_ind_cols))


	for i in range(len(min_ind_cols)):
		gdf_edges[min_ind_cols[i]] = 0
		gdf_edges[max_ind_cols[i]] = 0

	for iter_,path in save_paths_df.iterrows():
		min_path = path['min_edge_path']
		max_path = path['max_edge_path']
		
		gdf_edges.loc[gdf_edges['edge_id'].isin(min_path),min_ind_cols] += path[ch_min_ind_cols].values
		gdf_edges.loc[gdf_edges['edge_id'].isin(max_path),max_ind_cols] += path[ch_max_ind_cols].values
	
	# gdf_edges[min_ind_cols] = gdf_edges['min_vals'].apply(pd.Series)
	# gdf_edges[max_ind_cols] = gdf_edges['max_vals'].apply(pd.Series)
	# gdf_edges.drop('min_vals',axis=1,inplace=True)
	# gdf_edges.drop('max_vals',axis=1,inplace=True)


	for ind in industry_columns:
		gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)),axis = 1)
		gdf_edges[['min_{}'.format(ind),'max_{}'.format(ind)]] = gdf_edges['swap'].apply(pd.Series)
		gdf_edges.drop('swap',axis=1,inplace=True)

	if save_edges == True:
		gdf_edges.to_file(os.path.join(output_path,'weighted_edges_flows_national_{0}.shp'.format(transport_mode)))

	del gdf_edges, save_paths_df

'''
Create the database connection
'''
def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	'''
	Get the modal shares
	'''
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Airports','airnetwork'),('Waterways','waterways')]
	modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Airports','airnetwork'),('Waterways','waterways'),('Waterways','waterways')]
	modes_file_paths = [('Roads','national_roads')]
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail')]
	modes = ['road','rail','air','inland','coastal']
	veh_wt = [20,800,0,800,1200]
	usage_factors = [(0,0),(0,0),(0,0),(0.2,0.25),(0.2,0.25)]
	speeds = [(0,0),(40,60),(700,800),(9,20),(9,20)]
	mode_cols = ['road','rail','air','inland','coastal']
	new_mode_cols = ['o','d','road','rail','air','inland','coastal']

	shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles')
	od_output_excel = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_ods.xlsx')

	flow_output_excel = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_paths_roads_90_percent.xlsx')
	excl_wrtr = pd.ExcelWriter(flow_output_excel)

	md_prop_file = os.path.join(data_path,'mode_properties','mode_costs.xlsx')
	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	ind_cols = ['cement','coal','constructi','fertilizer','fishery','manufactur','acof','cash','cass','maiz','pepp','rcof',
				'rubb','swpo','teas','meat','rice','petroluem','steel','sugar','wood','tons']

	for m in range(len(modes_file_paths)):
		mode_data_path = os.path.join(data_path,modes_file_paths[m][0],modes_file_paths[m][1])
		for file in os.listdir(mode_data_path):
			try:
				if file.endswith(".shp") and 'edges' in file.lower().strip():
					edges_in = os.path.join(mode_data_path, file)
			except:
				return ('Network nodes and edge files necessary')
			
		if modes[m] == 'road': 
			G =  national_road_shapefile_to_network(edges_in,rd_prop_file)
			gdf_edges = national_road_shapefile_to_dataframe(edges_in,rd_prop_file)
		else:
			G =  network_shapefile_to_network(edges_in,md_prop_file,modes[m],speeds[m][0],speeds[m][1])
			gdf_edges = network_shapefile_to_dataframe(edges_in,md_prop_file,modes[m],speeds[m][0],speeds[m][1])

		# G = add_igraph_generalised_costs_network(G,1,veh_wt[m],usage_factors[m][0],usage_factors[m][1])
		nodes_name = np.asarray([x['name'] for x in G.vs])
		nodes_index = np.asarray([x.index for x in G.vs])
		node_dict = dict(zip(nodes_name,nodes_index))

		# print ([x for x in G.vs])
		# print (len(nodes_name))
		# od_output_csv = os.path.join(output_path,'flow_mapping_paths','national_scale_{}_ods.csv'.format(modes[m]))
		# all_ods = pd.read_csv(od_output_csv)
		all_ods = pd.read_excel(od_output_excel,sheet_name = modes[m])
		all_ods = all_ods[all_ods['max_tons'] > 0.5]
		print (len(all_ods.index))
		# if mode[m] != 'air' 
		# 	all_ods['min_vehicle_nums'] = np.maximum(1,np.ceil(all_ods['min_tons']/vehicle_wt))
		# 	all_ods['max_vehicle_nums'] = np.maximum(1,np.ceil(all_ods['max_tons']/vehicle_wt))
			# all_ods = all_ods[['origin','destination','min_croptons','max_croptons','min_netrev','max_netrev']]
		# network_od_paths_assembly(all_ods,node_dict,G,veh_wt[m],modes[m],ind_cols,gdf_edges,save_edges = True,output_path =shp_output_path,excel_writer =excl_wrtr)
		network_od_paths_assembly_changing_tonnages(all_ods,node_dict,G,veh_wt[m],usage_factors[m],modes[m],ind_cols,gdf_edges,save_edges = True,output_path =shp_output_path,excel_writer =excl_wrtr)
if __name__ == '__main__':
	main()