# -*- coding: utf-8 -*-
"""
Created on Sat Jul 14 15:41:39 2018

@author: elcok
"""

import geopandas as gpd
import pandas as pd
import os
import json
import igraph as ig
import numpy as np

from geopy.distance import vincenty
from boltons.iterutils import pairwise

def load_config():
    
   # Define current directory and data directory
    config_path = os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', '..','..', 'config.json')
    )
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    data_path = config['paths']['data']
    calc_path = config['paths']['calc']
    output_path = config['paths']['output']
    
    return data_path,calc_path,output_path

def line_length(line, ellipsoid='WGS-84'):
    """Length of a line in meters, given in geographic coordinates.

    Adapted from https://gis.stackexchange.com/questions/4022/looking-for-a-pythonic-way-to-calculate-the-length-of-a-wkt-linestring#answer-115285

    Args:
        line: a shapely LineString object with WGS-84 coordinates.
        
        ellipsoid: string name of an ellipsoid that `geopy` understands (see http://geopy.readthedocs.io/en/latest/#module-geopy.distance).

    Returns:
        Length of line in meters.
    """
    if line.geometryType() == 'MultiLineString':
        return sum(line_length(segment) for segment in line)

    return sum(
        vincenty(a, b, ellipsoid=ellipsoid).kilometers
        for a, b in pairwise(line.coords)
    )

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

    edges['LENGTH'] = edges.geometry.apply(line_length)

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2],axis=1)
    
    # create network from edge file
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G #.clusters().giant()

def netrev_edges(province,start_points,end_points,G,save_edges = True,output_path=''):
    save_paths = []
    for iter_,place in start_points.iterrows():
        try:
            closest_center = end_points.loc[end_points['OBJECTID'] 
            == place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]
           
            pos0_i = G.vs[node_dict[place['NEAREST_G_NODE']]]
            pos1_i = G.vs[node_dict[closest_center]]
    
            path = G.get_shortest_paths(pos0_i,pos1_i,weights='LENGTH',output="epath")
    
            get_path = [G.es[n]['EDGE_ID'] for n in path][0]
            save_paths.append((place['netrev'],get_path))
        except:
            print(iter_)
                    
    all_edges = [x['EDGE_ID'] for x in G.es]
    all_edges_geom = [x['geometry'] for x in G.es]
    
    gdf_edges = gpd.GeoDataFrame(pd.DataFrame([all_edges,all_edges_geom]).T,crs='epsg:4326')
    gdf_edges.columns = ['EDGE_ID','geometry']
    
    gdf_edges['netrev'] = 0
    for path in save_paths:
        gdf_edges.loc[gdf_edges['EDGE_ID'].isin(path[1]),'netrev'] += path[0]
    
    if save_edges == True:
        gdf_edges.to_file(os.path.join(output_path,'weighted_edges_{}.shp'.format(province)))
    return gdf_edges
     
    
def province_clip(shape_in,province_geom):
    gdf = gpd.read_file(shape_in)
    return gdf.loc[gdf['geometry'].apply(lambda x: x.within(province_geom))].reset_index(drop=True)
 
def get_nearest_node(x,sindex_nodes,nodes):
    return nodes.loc[list(sindex_nodes.nearest(x.bounds[:2]))]['NODE_ID'].values[0]

def get_nearest_center(x,sindex_nodes,nodes):
    return nodes.loc[list(sindex_nodes.nearest(x.bounds[:2]))]['OBJECTID'].values[0]
    

def count_points_in_polygon(x,prov_pop_sindex):
    return len(list(prov_pop_sindex.intersection(x.bounds)))

def get_netrev(x,commune_sindex,prov_communes):
        return prov_communes.loc[list(commune_sindex.intersection(x.bounds[:2]))]['netrev_village'].values[0]



if __name__ == '__main__':
    
    data_path,calc_path,output_path = load_config()
    
# =============================================================================
#     #province to consider 
# =============================================================================
    province = 'Thanh Hoa'
    
# =============================================================================
#     # set all paths for all input files we are going to use
# =============================================================================
    edges_in = os.path.join(data_path,'Roads','thanhhoa_roads','vietbando_thanhhoa_edges.shp')
    nodes_in = os.path.join(data_path,'Roads','thanhhoa_roads','vietbando_thanhhoa_nodes.shp')
    population_points_in = os.path.join(data_path,'Points_of_interest','population_points.shp')
    commune_center_in = os.path.join(data_path,'Points_of_interest','commune_committees_points.shp')

    province_path = os.path.join(data_path,'Vietnam_boundaries','who_boundaries','who_provinces.shp')
    commune_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','commune_level_stats.shp')

    # load provinces and get geometry of the right province
    provinces = gpd.read_file(province_path)
    province_geom = provinces.loc[provinces.NAME_ENG == province].geometry.values[0]
        
    #clip all to province
    prov_pop = province_clip(population_points_in,province_geom)
    prov_commune_center = province_clip(commune_center_in,province_geom)
    prov_communes = province_clip(commune_path,province_geom)

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
    prov_pop['netrev'] = prov_pop.geometry.apply(lambda x: get_netrev(x,commune_sindex,prov_communes))
    
    # and use average if commune has no stats
    prov_pop.loc[prov_pop['netrev'] == 0,'netrev'] = prov_pop['netrev'].mean()
   
# =============================================================================
#     # get nearest node in network for all start and end points
# =============================================================================
    prov_pop['NEAREST_G_NODE'] = prov_pop.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes))
    prov_commune_center['NEAREST_G_NODE'] = prov_commune_center.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes))

# =============================================================================
#     # prepare for shortest path routing, we'll use the spatial index of the centers
#     # to find the nearest center for each population point
# =============================================================================
    sindex_commune_center = prov_commune_center.sindex
    prov_pop['NEAREST_C_CENTER'] = prov_pop.geometry.apply(lambda x: get_nearest_center(x,sindex_commune_center,prov_commune_center))

# =============================================================================
#     # load network
# =============================================================================
    G = shapefile_to_network(edges_in)

    nodes_name = np.asarray([x['name'] for x in G.vs])
    nodes_index = np.asarray([x.index for x in G.vs])
    node_dict = dict(zip(nodes_name,nodes_index))

# =============================================================================
#     # loop through population points
# =============================================================================
    edges_updated = netrev_edges(province,prov_pop,prov_commune_center,G,save_edges = True,output_path='')
