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
from scripts.transport_network_creation import national_shapefile_to_network, add_igraph_generalised_costs_roads

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

def network_od_paths_assembly(points_dataframe,node_dict,graph,vehicle_wt,transport_mode,save_edges = True,output_path ='',excel_writer =''):
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

		# 	min_veh_nums = [origin_df['min_vehicle_nums']]
		# 	max_veh_nums = [origin_df['max_vehicle_nums']]
		# 	get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
		# 	print (get_min_path)
		# else:
		# 	# print (origin_df)
		# 	destinations = origin_df.loc[origin,'destination'].values.tolist()

		# 	min_tons = origin_df.loc[origin,'min_tons'].values.tolist()
		# 	max_tons = origin_df.loc[origin,'max_tons'].values.tolist()

		# 	min_veh_nums = origin_df.loc[origin,'min_vehicle_nums'].values.tolist()
		# 	max_veh_nums = origin_df.loc[origin,'max_vehicle_nums'].values.tolist()

		# get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
		# get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,origin,destinations,'max_gcost','max_time')

		# # get_min_gcost = np.multiply(np.array(min_veh_nums),np.array(get_min_gcost))
		# # get_max_gcost = np.multiply(np.array(max_veh_nums),np.array(get_max_gcost))

		# save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,min_tons,max_tons,
		# 						get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost,min_veh_nums,max_veh_nums))
		# print ("done with {0} in network {1}".format(origin,transport_mode))
		try:
			origin_df =  points_dataframe.loc[origin]
			# print (list(origin_df['destination']),origin_df.shape)
			if len(origin_df.shape) == 1:
				# print (origin_df['destination'])
				destinations = [origin_df['destination']]

				min_tons = [origin_df['min_tons']]
				max_tons = [origin_df['max_tons']]

				min_veh_nums = [origin_df['min_vehicle_nums']]
				max_veh_nums = [origin_df['max_vehicle_nums']]
				# get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
				# print (get_min_path)
			else:
				# print (origin_df)
				destinations = origin_df['destination'].values.tolist()

				min_tons = origin_df['min_tons'].values.tolist()
				max_tons = origin_df['max_tons'].values.tolist()

				min_veh_nums = origin_df['min_vehicle_nums'].values.tolist()
				max_veh_nums = origin_df['max_vehicle_nums'].values.tolist()

			get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(graph,origin,destinations,'min_gcost','min_time')
			get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(graph,origin,destinations,'max_gcost','max_time')

			# get_min_gcost = np.multiply(np.array(min_veh_nums),np.array(get_min_gcost))
			# get_max_gcost = np.multiply(np.array(max_veh_nums),np.array(get_max_gcost))

			save_paths += list(zip([origin]*len(destinations),destinations,get_min_path,get_max_path,min_tons,max_tons,
								get_min_dist,get_max_dist,get_min_time,get_max_time,get_min_gcost,get_max_gcost,min_veh_nums,max_veh_nums))
			# print ("done with {0} in network {1}".format(origin,transport_mode))
		except:
			print(origin)


	cols = ['origin','destination','min_edge_path','max_edge_path','min_tons','max_tons',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost','min_vehicle_nums','max_vehicle_nums']
	save_paths_df = pd.DataFrame(save_paths,columns = cols)
	save_paths_df.to_excel(excel_writer,transport_mode + '_{}_tons'.format(int(vehicle_wt)),index = False)
	excel_writer.save()
	del save_paths_df

	all_edges = [x['edge_id'] for x in graph.es]
	all_veh = [x['vehicle_co'] for x in graph.es]
	all_edges_geom = [x['geometry'] for x in graph.es]
	
	gdf_edges = gpd.GeoDataFrame(pd.DataFrame([all_edges,all_veh,all_edges_geom]).T,crs='epsg:4326')
	gdf_edges.columns = ['edge_id','vehicle_co','geometry']
	gdf_edges['vehicle_co'] = pd.to_numeric(gdf_edges['vehicle_co'])
	
	gdf_edges['min_tons'] = 0
	gdf_edges['max_tons'] = 0
	gdf_edges['min_veh'] = 0
	gdf_edges['max_veh'] = 0

	for path in save_paths:
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[2]),'min_tons'] += path[4]
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[3]),'max_tons'] += path[5]
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[2]),'min_veh'] += path[-2]
		gdf_edges.loc[gdf_edges['edge_id'].isin(path[3]),'max_veh'] += path[-1]
	

	gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_tons','max_tons'),axis = 1)
	gdf_edges[['min_tons','max_tons']] = gdf_edges['swap'].apply(pd.Series)
	gdf_edges.drop('swap',axis=1,inplace=True)

	gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_veh','max_veh'),axis = 1)
	gdf_edges[['min_veh','max_veh']] = gdf_edges['swap'].apply(pd.Series)
	gdf_edges.drop('swap',axis=1,inplace=True)

	if save_edges == True:
		gdf_edges.to_file(os.path.join(output_path,'weighted_edges_flows_national_{0}_{1}_tons.shp'.format(transport_mode,int(vehicle_wt))))

	del gdf_edges

'''
Create the database connection
'''
def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	exchange_rate = 1.05*(1000000/21000)
	vehicle_wt = 20.0

	# population_points_in = os.path.join(data_path,'Points_of_interest','population_points.shp')
	# commune_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','commune_level_stats.shp')

	# crop_data_path = os.path.join(data_path,'Agriculture_crops','crop_data')
	# rice_month_file = os.path.join(data_path,'rice_atlas_vietnam','rice_production.shp')
	# crop_month_fields = ['P_Jan','P_Feb','P_Mar','P_Apr','P_May','P_Jun','P_Jul','P_Aug','P_Sep','P_Oct','P_Nov','P_Dec']
	# crop_names = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']

	'''
	Get the modal shares
	'''
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Airports','airnetwork'),('Waterways','waterways')]
	modes_file_paths = [('Roads','national_roads')]
	modes = ['road','rail','air','water']
	mode_cols = ['road','rail','air','inland','coastal']
	new_mode_cols = ['o','d','road','rail','air','water']
	shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles')
	od_output_excel = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_ods.xlsx')

	flow_output_excel = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_paths.xlsx')
	excl_wrtr = pd.ExcelWriter(flow_output_excel)

	rd_prop_file = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')

	for m in range(len(modes_file_paths)):
		mode_data_path = os.path.join(data_path,modes_file_paths[m][0],modes_file_paths[m][1])
		for file in os.listdir(mode_data_path):
			try:
				if file.endswith(".shp") and 'edges' in file.lower().strip():
					edges_in = os.path.join(mode_data_path, file)
				if file.endswith(".shp") and 'nodes' in file.lower().strip():
					nodes_in = os.path.join(mode_data_path, file)
			except:
				print ('Network nodes and edge files necessary')
			
		if modes[m] == 'road': 
			G =  national_shapefile_to_network(edges_in,rd_prop_file)
			G = add_igraph_generalised_costs_roads(G,1,vehicle_wt)
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
		all_ods['min_vehicle_nums'] = np.maximum(1,np.ceil(all_ods['min_tons']/vehicle_wt))
		all_ods['max_vehicle_nums'] = np.maximum(1,np.ceil(all_ods['max_tons']/vehicle_wt))
		# all_ods = all_ods[['origin','destination','min_croptons','max_croptons','min_netrev','max_netrev']]
		network_od_paths_assembly(all_ods,node_dict,G,vehicle_wt,modes[m],save_edges = True,output_path =shp_output_path,excel_writer =excl_wrtr)
			

if __name__ == '__main__':
	main()