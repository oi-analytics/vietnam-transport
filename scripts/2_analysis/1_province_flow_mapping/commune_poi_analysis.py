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
def get_nearest_node(x,sindex_nodes,nodes):
    return nodes.loc[list(sindex_nodes.nearest(x.bounds[:2]))]['NODE_ID'].values[0]

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

    

    
def province_clip(shape_in,province_geom):
    
    gdf = gpd.read_file(shape_in)
    return gdf.loc[gdf['geometry'].apply(lambda x: x.within(province_geom))]
    
if __name__ == '__main__':
    
    data_path,calc_path,output_path = load_config()
    
    #province to consider 
    province = 'Thanh Hoa'
    
    # set all paths for all input files we are going to use
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

    # load nodes
    nodes = gpd.read_file(nodes_in)
    sindex_nodes = nodes.sindex
    
    # get nearest node in network for all start and end points
    prov_pop['G_NODE'] = prov_pop.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes))
    prov_commune_center['G_NODE'] = prov_commune_center.geometry.apply(lambda x: get_nearest_node(x,sindex_nodes,nodes))

    # prepare for shortest path routing, we'll use the spatial index of the centers
    # to find the nearest center for each population point
    sindex_commune_center = prov_commune_center.sindex
 
    # load network
    G = shapefile_to_network(edges_in)
    nodes = np.asarray([x['name'] for x in G.vs])
 
    
    for iter_,place in prov_pop.iterrows():
        
        if iter_ > 18819:
            break
    
#        closest_center = 
        
        
        pos0_i = [x for x in G.vs if x['name'] == place['G_NODE']][0]
        
        
#        pos1_i = sg.vs[closest_node(destinations[route[1]],nodes)]
#    path = sg.get_shortest_paths(pos0_i,pos1_i,weights='t_time',output="epath")

    