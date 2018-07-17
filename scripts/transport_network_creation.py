# -*- coding: utf-8 -*-
"""
Python script to create transport networks in Vietnam  
Created on Wed June 27 2018

@author: Raghav Pant
"""

import pandas as pd
import os
import psycopg2
import networkx as nx
import csv
import igraph as ig 
import numpy as np
import geopandas as gpd
from scripts.utils import line_length

def assign_assumed_width_to_province_roads(x):
	'''
	Assign widths to roads assets in Vietnam
	The widths are assigned based on our understanding of: 
	1. The reported width in the data which is not reliable
	2. A design specification based understanding of the assumed width based on ranges of values
	
	Inputs are:
	asset_width - Numeric value for width of asset
	
	Outputs are:
	modified_width - assigned width of the road asset based on design specifications
	'''
	if 0 <= x.WIDTH < 4.25:
		return 3.5
	elif 4.25 <= x.WIDTH < 6.0:
		return 5.0
	elif 6.0 <= x.WIDTH < 8.0:
		return 7.0
	elif 8.0 <= x.WIDTH < 11.5:
		return 9.0
	elif 11.5 <= x.WIDTH < 17.5:
		return 14.0
	elif 17.5 <= x.WIDTH < 24.5:
		return 21.0
	elif 24.5 <= x.WIDTH < 100:
		return 9.0
	else:
		return x.WIDTH

def assign_asset_type_to_province_roads(x):
	'''
	Assign asset types to roads assets in Vietnam
	The types are assigned based on our understanding of: 
	1. The reported asset code in the data
	
	Inputs are:
	asset code - Numeric value for code of asset
	
	Outputs are:
	asset type - Which is either of (Bridge,Dam,Culvert,Tunnel,Spillway,Road)
	'''
	if x.CODE in (12,25):
		return 'Bridge'
	elif x.CODE == (23):
		return 'Dam'
	elif x.CODE == (24):
		return 'Culvert'
	elif x.CODE == (26):
		return 'Tunnel'
	elif x.CODE == (27):
		return 'Spillway'
	else:
		return 'Road'


def assign_minmax_travel_speeds_province_roads_apply(x):
	'''
	Assign travel speeds to roads assets in Vietnam
	The speeds are assigned based on our understanding of: 
	1. The types of assets
	2. The levels of classification of assets: 0-National,1-Provinical,2-Local,3-Other
	3. The terrain where the assets are located: Flat or Mountain or No information
		
	Inputs are:
	asset_code - Numeric code for type of asset
	asset_level - Numeric code for level of asset
	asset_terrain - String value of the terrain of asset
		
	Outputs are:
	speed_min - Minimum assigned speed in km/hr
	speed_max - Maximum assigned speed in km/hr
	'''
	asset_code = x.CODE
	asset_level = x.LEVEL
	asset_terrain='flat'

	if (not asset_terrain) or (asset_terrain == 'flat'):
		if asset_code == 17: # This is an expressway
			return 100,120
		elif asset_code in (15,4): # This is a residential road or a mountain pass
			return 40,60
		elif asset_level == 0: # This is any other national network asset
			return 80,100
		elif asset_level == 1:# This is any other provincial network asset
			return 60,80
		elif asset_level == 2: # This is any other local network asset
			return 40,60
		else:			# Anything else not included above
			return 20,40

	else:
		if asset_level < 3:
			return 40, 60
		else:
			return 20,40


def province_shapefile_to_network(edges_in):
	"""
	input parameters:
		edges_in : string of path to edges file/network file. 
		
	output:
		SG: connected graph of the shapefile
	
	"""
	
	edges = gpd.read_file(edges_in)
	
	# assign minimum and maximum speed to network
	edges['SPEED'] = edges.apply(assign_minmax_travel_speeds_province_roads_apply,axis=1)
	edges[['MIN_SPEED', 'MAX_SPEED']] = edges['SPEED'].apply(pd.Series)
	edges.drop('SPEED',axis=1,inplace=True)

	# correct the widths of the road assets
	edges['MOD_WIDTH'] = edges.apply(assign_assumed_width_to_province_roads,axis=1)
	# edges.drop('WIDTH',axis=1,inplace=True)

	# create an asset type column
	edges['ASSET_TYPE'] = edges.apply(assign_asset_type_to_province_roads,axis=1)

	# get the right linelength
	edges['LENGTH'] = edges.geometry.apply(line_length)

	# assign costs to the edges
	# edges['MAX_COST'] = cost_param*edges['LENGTH']/edges['MIN_SPEED']
	# edges['MIN_COST'] = cost_param*edges['LENGTH']/edges['MAX_SPEED']
	# make sure that From and To node are the first two columns of the dataframe
	# to make sure the conversion from dataframe to igraph network goes smooth
	edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2],axis=1)
	
	# create network from edge file
	G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

	# only keep connected network
	return G.clusters().giant()

def add_igraph_time_costs_province_roads(G,cost_param):
	G.es['MAX_COST'] = list(cost_param*(np.array(G.es['LENGTH'])/np.array(G.es['MAX_SPEED'])))
	G.es['MIN_COST'] = list(cost_param*(np.array(G.es['LENGTH'])/np.array(G.es['MIN_SPEED'])))

	return G

'''
Functions we are not using at present for provincial analysis. Will clean them later
'''

def assign_assumed_width_to_province_roads_from_file(asset_width,width_range_list):
	'''
	Assign widths to roads assets in Vietnam
	The widths are assigned based on our understanding of: 
	1. The reported width in the data which is not reliable
	2. A design specification based understanding of the assumed width based on ranges of values
	
	Inputs are:
	asset_width - Numeric value for width of asset
	width_range_list - List of tuples containing (from_width,to_width,assumed_width)
	
	Outputs are:
	assumed_width - assigned width of the raod asset based on design specifications
	'''
	
	assumed_width = asset_width
	for width_vals in width_range_list:
		if width_vals[0] <= assumed_width <= width_vals[1]:
			asset_width = width_vals[2]
			break

	return assumed_width

def assign_minmax_travel_speeds_province_roads(asset_code,asset_level,asset_terrain):
	'''
	Assign travel speeds to roads assets in Vietnam
	The speeds are assigned based on our understanding of: 
	1. The types of assets
	2. The levels of classification of assets: 0-National,1-Provincial,2-Local,3-Other
	3. The terrain where the assets are located: Flat or Mountain or No information
		
	Inputs are:
	asset_code - Numeric code for type of asset
	asset_level - Numeric code for level of asset
	asset_terrain - String value of the terrain of asset
		
	Outputs are:
	speed_min - Minimum assigned speed in km/hr
	speed_max - Maximum assigned speed in km/hr
	'''

	if (not asset_terrain) or (asset_terrain == 'flat'):
		if asset_code == 17: # This is an expressway
			return 100,120
		elif asset_code in (15,4): # This is a residential road or a mountain pass
			return 40,60
		elif asset_level == 0: # This is any other national network asset
			return 80,100
		elif asset_level == 1:# This is any other provincial network asset
			return 60,80
		elif asset_level == 2: # This is any other local network asset
			return 40,60
		else:			# Anything else not included above
			return 20,40

	else:
		if asset_level < 3:
			return 40, 60
		else:
			return 20,40

def shapefile_to_network(edges_in,path_width_table):
	"""
	input parameters:
		edges_in : string of path to edges file/network file. 
		
	output:
		SG: connected graph of the shapefile
	"""
	
	edges = gpd.read_file(edges_in)
	
	# assign minimum and maximum speed to network
	edges['SPEED'] = edges.apply(assign_minmax_travel_speeds_roads_apply,axis=1)
	edges[['MIN_SPEED', 'MAX_SPEED']] = edges['SPEED'].apply(pd.Series)
	edges.drop('SPEED',axis=1,inplace=True)
	
	# get the right linelength
	edges['LENGTH'] = edges.geometry.apply(line_length)

	# get the width of edges
	width_range_list = [tuple(x) for x in pd.read_excel(path_width_table,sheet_name ='widths').values]
	
	edges['WIDTH'] = edges.WIDTH.apply(lambda x: assign_assumed_width_to_roads(x,width_range_list))
	
	# make sure that From and To node are the first two columns of the dataframe
	# to make sure the conversion from dataframe to igraph network goes smooth
	edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2],axis=1)
	
	# create network from edge file
	G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

	# only keep connected network
	return G.clusters().giant()

def assign_travel_times_and_variability(speed_attributes,variability_attributes,mode_type,variability_location,mode_attribute,distance):
	travel_time = 0
	variability_time = 0

	st = [s[2] for s in speed_attributes if s[0].lower().strip() == mode_type and s[1].lower().strip() == mode_attribute]
	if len(st) > 0:
		travel_time = distance/sum(st)
		vt = [v[3] for v in variability_attributes if v[0] == mode_type and v[1] == variability_location]
		# print (vt)
		variability_time = (1.0*sum(vt)/100)*travel_time

	return travel_time, variability_time

def assign_network_dictionary(network_dictionary,edge,from_node,to_node,distance,speed):
	
	network_dictionary['edge'].append(edge)
	network_dictionary['from_node'].append(from_node)
	network_dictionary['to_node'].append(to_node)
	network_dictionary['distance'].append(distance)
	network_dictionary['speed'].append(speed)
	network_dictionary['travel_cost'].append(0.019*distance/speed)

	return network_dictionary

def create_network_dictionary(network_dictionary,layer_name, edge_id, from_node_id, to_node_id,edge_speed_attribute,edge_geom_attribute,cursor,connection):
	'''
	input parameters:
	layer_name: sql layer name to extraxt information from
	edge_id: unique ID of edge in sql layer
	from_node_id: unqiue source node ID corresponding to edge ID in sql layer
	to_node_id: unique target node ID corresponding to edge ID in sql layer
	edge_geom_attribute: for calculating the length of the edge
	cursor: the postGIS cursor 
	connection: the postGIS connection

	output:
	network_dict = {'edgeid':[],'from_node':[],'to_node':[],'distance':[]}

	'''
	sql_query = '''select {0}, {1}, {2}, {3}, st_length({4}::geography)/1000 from {5}
				'''.format(edge_id,from_node_id,to_node_id,edge_speed_attribute,edge_geom_attribute,layer_name)

	cursor.execute(sql_query)
	read_layer = cursor.fetchall()
	for row in read_layer:
		if str(row[0]).isdigit():
			e_id = int(row[0])
		else:
			e_id = row[0]
		
		if str(row[1]).isdigit():
			fn_id = int(row[1])
		else:
			fn_id = row[1]
		
		if str(row[2]).isdigit():
			tn_id = int(row[2])
		else:
			tn_id = row[2]

		sp = float(row[3])
		lgth = float(row[4])

			
		network_dictionary = assign_network_dictionary(network_dictionary,e_id,fn_id,tn_id,lgth,sp)

	return network_dictionary

def create_networkx_topology(network_dictionary,tonnage,teus):
	# all_net_dict = {'edge':[],'from_node':[],'to_node':[],'waiting_cost':[],'travel_cost':[],'transport_price_ton':[],
	# 			'transport_price_teu':[],'variability_cost':[]}

	G = nx.Graph()
	edges = network_dictionary['edge']
	from_nodes = network_dictionary['from_node']
	to_nodes = network_dictionary['to_node']
	wait_costs = network_dictionary['waiting_cost']
	travel_costs = network_dictionary['travel_cost']
	price_tons = network_dictionary['transport_price_ton']
	price_teus = network_dictionary['transport_price_teu']
	variable_costs = network_dictionary['variability_cost'] 

	for e in range(len(edges)):
		generalised_cost = wait_costs[e] + travel_costs[e] + tonnage*price_tons[e] + teus*price_teus[e] + variable_costs[e]
		G.add_edge(from_nodes[e],to_nodes[e], edge = edges[e], cost = generalised_cost)

	return G

def create_igraph_topology(network_dictionary):
	# all_net_dict = {'edge':[],'from_node':[],'to_node':[],'waiting_cost':[],'travel_cost':[],'transport_price_ton':[],
	# 			'transport_price_teu':[],'variability_cost':[]}

	G = ig.Graph()
	edges = network_dictionary['edge']
	from_nodes = network_dictionary['from_node']
	to_nodes = network_dictionary['to_node']

	unique_nodes = list(set(from_nodes + to_nodes))
	igraph_ids = list(range(len(unique_nodes)))
	G.add_vertices(igraph_ids)

	edge_list = []
	for e in range(len(edges)):
		fn = from_nodes[e]
		tn = to_nodes[e]
		ig_fn = igraph_ids[unique_nodes.index(fn)]
		ig_tn = igraph_ids[unique_nodes.index(tn)]

		edge_list.append((ig_fn,ig_tn))

	G.add_edges(edge_list)

	G.vs['node'] = unique_nodes

	G.es['from_node'] = from_nodes
	G.es['to_node'] = to_nodes
	G.es['edge'] = edges
	G.es['distance'] = network_dictionary['distance']
	G.es['travel_cost'] = network_dictionary['travel_cost'] 

	return G

def add_igraph_costs(G,tonnage,teus):
	# all_net_dict = {'edge':[],'from_node':[],'to_node':[],'waiting_cost':[],'travel_cost':[],'transport_price_ton':[],
	# 			'transport_price_teu':[],'variability_cost':[]}

	if teus > 0:
		generalised_cost = np.array(G.es['waiting_cost']) + np.array(G.es['travel_cost']) + teus*np.array(G.es['price_teus']) + np.array(G.es['variable_costs'])
	else:
		generalised_cost = np.array(G.es['waiting_cost'])+ np.array(G.es['travel_cost']) + tonnage*np.array(G.es['price_tons']) + np.array(G.es['variable_costs'])

	
	G.es['cost'] = list(generalised_cost)

	return G

def get_networkx_edges(G, node_path,data_cond):
	node_tup = zip(node_path[:-1],node_path[1:])
	t_list = []
	for tup in node_tup:
		t = [d for (u,v,d) in G.edges(data = data_cond) if (u,v) == tup or (v,u) == tup]
		if len(t) > 0:
			t_list.append(t[0])
	return t_list

def get_igraph_edges(G,edge_path,data_cond):
	t_list = [G.es[n][data_cond] for n in edge_path]
	return (t_list)

