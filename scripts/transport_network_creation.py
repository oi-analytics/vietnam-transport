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

def assign_assumed_width_to_roads(asset_width,width_range_list):
	'''
	==========================================================================================
	Assign widths to roads assets in Vietnam
	The widths are assigned based on our understanding of: 
	1. The reported width in the data which is not reliable
	2. A design specification based understanding of the assumed width based on ranges of values
	
	Inputs are:
	asset_width - Numeric value for width of asset
	width_range_list - List of tuples containing (from_width,to_width,assumed_width)
	
	Outputs are:
	assumed_width - assigned width of the raod asset based on design specifications
	============================================================================================ 
	'''
	assumed_width = asset_width
	for width_vals in width_range_list:
		if width_range_list[0] <= assumed_width <= width_range_list[1]:
			asset_width = width_range_list[2]
			break

	return assumed_width


def assign_minmax_travel_speeds_roads(asset_code,asset_level,asset_terrain):
    '''
    ====================================================================================
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
    ==================================================================================== 
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

def assign_minmax_travel_speeds_roads_apply(x):
    '''
    ====================================================================================
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
    ==================================================================================== 
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

def shapefile_to_network(edges_in):
    """
    input parameters:
        edges_in : string of path to edges file/network file. 
        
    output:
        SG: connected graph of the shapefile
    
    """
    
    edges = gpd.read_file(edges_in)
    
    # assign minimum and maximum speed to network
    edges['speed'] = edges.apply(assign_minmax_travel_speeds_roads_apply,axis=1)
    edges[['min_speed', 'max_speed']] = edges['speed'].apply(pd.Series)
    edges.drop('speed',axis=1,inplace=True)


    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2],axis=1)
    
    # create network from edge file
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G.clusters().giant()


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
