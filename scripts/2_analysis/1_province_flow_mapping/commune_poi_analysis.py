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

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from scripts.utils import load_config,extract_value_from_gdf,get_nearest_node,gdf_clip,count_points_in_polygon
from scripts.transport_network_creation import shapefile_to_network

def netrev_edges(region_name,start_points,end_points,graph,save_edges = True,output_path=''):
    """
   ====================================================================================
    Assign net revenue to roads assets in Vietnam
    	
    Inputs are:
    start_points - GeoDataFrame of start points for shortest path analysis.
    end_points - GeoDataFrame of potential end points for shorest path analysis.
    G - iGraph network of the province.
    save_edges - 
    	
    Outputs are:
    Shapefile with all edges and the total net reveneu transferred along each edge
    GeoDataFrame of total net revenue transferred along each edge
    ==================================================================================== 
    """
    save_paths = []
    for iter_,place in start_points.iterrows():
        try:
            closest_center = end_points.loc[end_points['OBJECTID'] 
            == place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]
           
            pos0_i = graph.vs[node_dict[place['NEAREST_G_NODE']]]
            pos1_i = graph.vs[node_dict[closest_center]]
    
            path = graph.get_shortest_paths(pos0_i,pos1_i,weights='LENGTH',output="epath")
    
            get_path = [graph.es[n]['EDGE_ID'] for n in path][0]
            save_paths.append((place['netrev'],get_path))
        except:
            print(iter_)
                    
    all_edges = [x['EDGE_ID'] for x in graph.es]
    all_edges_geom = [x['geometry'] for x in graph.es]
    
    gdf_edges = gpd.GeoDataFrame(pd.DataFrame([all_edges,all_edges_geom]).T,crs='epsg:4326')
    gdf_edges.columns = ['EDGE_ID','geometry']
    
    gdf_edges['netrev'] = 0
    for path in save_paths:
        gdf_edges.loc[gdf_edges['EDGE_ID'].isin(path[1]),'netrev'] += path[0]
    
    if save_edges == True:
        gdf_edges.to_file(os.path.join(output_path,'weighted_edges_{}.shp'.format(region_name)))
    return gdf_edges
    

if __name__ == '__main__':
    
    data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']
    
# =============================================================================
#     #province to consider 
# =============================================================================
#    province = 'Thanh Hoa'
    province = 'Binh Dinh'
#    province = 'Lao Cai'
    
# =============================================================================
#     # set all paths for all input files we are going to use
# =============================================================================
    province_name = province.replace(' ','').lower()
    
    edges_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))
    nodes_in = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_nodes.shp'.format(province_name))
    population_points_in = os.path.join(data_path,'Points_of_interest','population_points.shp')
    commune_center_in = os.path.join(data_path,'Points_of_interest','commune_committees_points.shp')

    province_path = os.path.join(data_path,'Vietnam_boundaries','who_boundaries','who_provinces.shp')
    commune_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','commune_level_stats.shp')

    # load provinces and get geometry of the right province
    provinces = gpd.read_file(province_path)
    province_geom = provinces.loc[provinces.NAME_ENG == province].geometry.values[0]
        
    #clip all to province
    prov_pop = gdf_clip(population_points_in,province_geom)
    prov_commune_center = gdf_clip(commune_center_in,province_geom)
    prov_communes = gdf_clip(commune_path,province_geom)

# =============================================================================
#     # load nodes and edges
# =============================================================================
    nodes = gpd.read_file(nodes_in)
    sindex_nodes = nodes.sindex
    
# =============================================================================
#     # get revenue values for each village
# =============================================================================
    
    # first create sindex of all villages to count number of villages in commune
    prov_pop_sindex = prov_pop.sindex
    
    # create new column in prov_communes with amount of villages
    prov_communes['n_villages'] = prov_communes.geometry.apply(lambda x: count_points_in_polygon(x,prov_pop_sindex)) 
    prov_communes['netrev_village'] = (prov_communes['netrevenue']*prov_communes['nfirm'])/prov_communes['n_villages'] 

    commune_sindex = prov_communes.sindex
    # give each village a net revenue based on average per village in commune
    prov_pop['netrev'] = prov_pop.geometry.apply(lambda x: extract_value_from_gdf(x,commune_sindex,prov_communes,'netrev_village'))
    
    # and use average if commune has no stats
    prov_pop.loc[prov_pop['netrev'] == 0,'netrev'] = prov_pop['netrev'].mean()
   
# =============================================================================
#     # get nearest node in network for all start and end points
# =============================================================================
    prov_pop['NEAREST_G_NODE'] = prov_pop.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes,'NODE_ID'))
    prov_commune_center['NEAREST_G_NODE'] = prov_commune_center.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes,'NODE_ID'))

# =============================================================================
#     # prepare for shortest path routing, we'll use the spatial index of the centers
#     # to find the nearest center for each population point
# =============================================================================
    sindex_commune_center = prov_commune_center.sindex
    prov_pop['NEAREST_C_CENTER'] = prov_pop.geometry.apply(lambda x: get_nearest_node(x,sindex_commune_center,prov_commune_center,'OBJECTID'))

# =============================================================================
#     # load network
# =============================================================================
    G = shapefile_to_network(edges_in)

    nodes_name = np.asarray([x['name'] for x in G.vs])
    nodes_index = np.asarray([x.index for x in G.vs])
    node_dict = dict(zip(nodes_name,nodes_index))

# =============================================================================
#     # get updated edges
# =============================================================================
    edges_updated = netrev_edges(province_name,prov_pop,prov_commune_center,G,save_edges = True,output_path=output_path)
