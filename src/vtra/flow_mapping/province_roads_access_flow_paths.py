"""Mapping the commune access OD node level matrix values to road network paths in Provinces
For all roads in the Provinces:
    ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']

The code estimates 2 values - A MIN and a MAX value of flows between each selected OD node pair
    Based on MIN-MAX generalised costs estimates

Input data requirements
-----------------------
1. Correct paths to all files and correct input parameters
2. MINIMUM MANDATORY DATA SPECIFICATION
    A. Excel file with mode sheets containing network graph structure and attributes
        Should contain following column names and attributes:
            edge_id - String/Integer/Float Edge ID
            from_node - String/Integer/Float node ID that should be present in node_id column
            to_node - String/Integer/Float node ID that should be present in node_id column
            length - Float length of edge in km
            min_time - Float minimum time of travel in hours on edge 
            max_time - Float maximum time of travel in hours on edge  
            min_time_cost - Float minimum cost of time in USD on edge   
            max_time_cost - Float maximum cost of time in USD on edge 
            min_tariff_cost - Float minimum tariff cost in USD on edge   
            max_tariff_cost - Float maximum tariff cost in USD on edge
             
    B. Edge shapefiles for all national-scale networks
        Should contain following column names and attributes:
            edge_id - String/Integer/Float Edge ID
            geometry - Shapely LineString geometry of edges

    C. Excel file with mode sheets containing node-level OD values
        Should contain following column names and attributes:
            origin - String node ID of Origin
            destination - String node ID of Destination
            min_netrev - Float values of miniimum daily OD Net Revenue in USD
            max_netrev - Float values of maximum daily OD Net Revenue in USD
            min_tons -  Float values of minimum daily OD in tons
            max_tons - Float values of maximum daily OD in tons

Results
-------
1. Excel sheets with results of flow mapping based on MIN-MAX generalised costs estimates:
    origin - String node ID of Origin
    destination - String node ID of Destination
    min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
    max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
    min_netrev - Float values of estimated daily Net Revenue for paths with minimum generalised cost flows
    max_netrev - Float values of estimated daily Net Revenue for paths with maximum generalised cost flows
    min_croptons - Float values of estimated daily crop tonnage for paths with minimum generalised cost flows
    max_croptons - Float values of estimated daily crop tonnage for paths with maximum generalised cost flows
    min_distance - Float values of estimated distance for paths with minimum generalised cost flows
    max_distance - Float values of estimated distance for paths with maximum generalised cost flows
    min_time - Float values of estimated time for paths with minimum generalised cost flows
    max_time - Float values of estimated time for paths with maximum generalised cost flows
    min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
    max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
    min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
    max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows

2. Shapefiles with all flows on edges mapping based on MIN-MAX generalised costs estimates
    edge_id - String/Integer/Float Edge ID
    geometry - Shapely LineString geomtry of edges
    min_netrev - Float values of estimated daily Net Revenue in USD on edges
    max_netrev - Float values of estimated daily Net Revenue in USD on edges
    min_tons - Float values of estimated daily crops in tons on edges
    max_tons - Float values of estimated daily crops in tons on edges  
"""

import ast
import itertools
import math
import operator
import os
import subprocess
import sys

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd
from shapely import wkt
from shapely.geometry import Point
from vtra.transport_network_creation import *
from vtra.utils import *

def network_od_path_estimations(graph, source, target, cost_criteria, time_criteria):
    """
    Estimate the paths, distances, times, and costs for given OD pair

    Parameters
    ---------
    graph - igraph network structure 
    source - String/Float/Integer name of Origin node ID
    source - String/Float/Integer name of Destination node ID
    cost_criteria - String name of generalised cost criteria to be used: min_gcost or max_gcost
    time_criteria - String name of time criteria to be used: min_time or max_time    

    Outputs
    -------
    edge_path_list - List of lists of Strings/Floats/Integers of edge ID's in routes
    path_dist_list - List of float values of estimated distances of routes
    path_time_list - List of float values of estimated times of routes
    path_gcost_list - List of float values of estimated generalised costs of routes

    """
    paths = graph.get_shortest_paths(source, target, weights=cost_criteria, output="epath")

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

    return edge_path_list, path_dist_list, path_time_list, path_gcost_list


def network_od_paths_assembly(points_dataframe, graph, vehicle_wt, region_name, excel_writer=''):
    """
    Assemble estimates of OD paths, distances, times, costs and tonnages on networks

    Parameters
    ---------
    points_dataframe - Pandas DataFrame of OD nodes and their tonnages
    graph - igraph network structure 
    vehicle_wt - Float unit weight of vehicle
    region_name - String name of Province
    excel_writer - Name of the excel writer to save Pandas dataframe to Excel file     

    Outputs
    -------
    save_paths_df - Pandas DataFrame 
        With columns:
            origin - String node ID of Origin
            destination - String node ID of Destination
            min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
            max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
            min_netrev - Float values of estimated netrevenue for paths with minimum generalised cost flows
            max_netrev - Float values of estimated netrevenue for paths with maximum generalised cost flows
            min_croptons - Float values of estimated crop tons for paths with minimum generalised cost flows
            max_croptons - Float values of estimated crop tons for paths with maximum generalised cost flows
            min_distance - Float values of estimated distance for paths with minimum generalised cost flows
            max_distance - Float values of estimated distance for paths with maximum generalised cost flows
            min_time - Float values of estimated time for paths with minimum generalised cost flows
            max_time - Float values of estimated time for paths with maximum generalised cost flows
            min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
            max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
            min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
            max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows

    """
    save_paths = []
    points_dataframe = points_dataframe.set_index('origin')
    origins = list(set(points_dataframe.index.values.tolist()))
    for origin in origins:
        try:
            destinations = points_dataframe.loc[[origin], 'destination'].values.tolist()

            min_croptons = points_dataframe.loc[[origin], 'min_croptons'].values.tolist()
            max_croptons = points_dataframe.loc[[origin], 'max_croptons'].values.tolist()

            min_rev = points_dataframe.loc[[origin], 'min_netrev'].values.tolist()
            max_rev = points_dataframe.loc[[origin], 'max_netrev'].values.tolist()

            min_veh_nums = points_dataframe.loc[[origin], 'min_vehicle_nums'].values.tolist()
            max_veh_nums = points_dataframe.loc[[origin], 'max_vehicle_nums'].values.tolist()

            get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(
                graph, origin, destinations, 'min_gcost', 'min_time')
            get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(
                graph, origin, destinations, 'max_gcost', 'max_time')


            save_paths += list(zip([origin]*len(destinations), destinations, get_min_path, get_max_path, min_rev, max_rev, min_croptons, max_croptons,
                                   get_min_dist, get_max_dist, get_min_time, get_max_time, get_min_gcost, get_max_gcost, min_veh_nums, max_veh_nums))
            print("done with {0} in province {1}".format(origin, region_name))
        except:
            print('* no path between {}-{}'.format(origin,destinations))

    cols = ['origin', 'destination', 'min_edge_path', 'max_edge_path', 'min_netrev', 'max_netrev', 'min_croptons', 'max_croptons',
            'min_distance', 'max_distance', 'min_time', 'max_time', 'min_gcost', 'max_gcost', 'min_vehicle_nums', 'max_vehicle_nums']
    save_paths_df = pd.DataFrame(save_paths, columns=cols)
    save_paths_df.to_excel(excel_writer, region_name +
                           '_{}_tons'.format(int(vehicle_wt)), index=False)
    excel_writer.save()
    del save_paths_df

    return save_paths

def write_province_flow_paths_to_network_shapefile(save_paths,gdf_edges, save_edges=True, shape_output_path=''):
    """
    Write results to Shapefiles

    Parameters
    ---------
    save_paths - List of lists of OD flow paths and their min-max tonnages and revenues
    region_name - String name of province
    gdf_edges - GeoDataFrame of network edge set
    save_Edges - Boolean condition to tell code to save created edge shapefile
    shape_output_path - Path where the output shapefile will be stored 

    Outputs
    -------
    gdf_edges - Shapefile 
        With minimum and maximum tonnage and net revenue flows for each edges in network
    """

    gdf_edges['min_netrev'] = 0
    gdf_edges['max_netrev'] = 0
    gdf_edges['min_tons'] = 0
    gdf_edges['max_tons'] = 0

    for path in save_paths:
        gdf_edges.loc[gdf_edges['edge_id'].isin(path[2]), 'min_netrev'] += path[4]
        gdf_edges.loc[gdf_edges['edge_id'].isin(path[3]), 'max_netrev'] += path[5]
        gdf_edges.loc[gdf_edges['edge_id'].isin(path[2]), 'min_tons'] += path[6]
        gdf_edges.loc[gdf_edges['edge_id'].isin(path[3]), 'max_tons'] += path[7]

    gdf_edges['swap'] = gdf_edges.apply(
        lambda x: swap_min_max(x, 'min_netrev', 'max_netrev'), axis=1)
    gdf_edges[['min_netrev', 'max_netrev']] = gdf_edges['swap'].apply(pd.Series)
    gdf_edges.drop('swap', axis=1, inplace=True)

    gdf_edges['swap'] = gdf_edges.apply(
        lambda x: swap_min_max(x, 'min_tons', 'max_tons'), axis=1)
    gdf_edges[['min_tons', 'max_tons']] = gdf_edges['swap'].apply(pd.Series)
    gdf_edges.drop('swap', axis=1, inplace=True)

    if save_edges == True:
        gdf_edges.to_file(shape_output_path,encoding='utf-8')


if __name__ == '__main__':
    """
    Specify the paths from where you want to read and write:
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    1. Names of the three Provinces: ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
        List of string types 
    2. Assumed terrains of the provinces: flat or mountanious
        List of string types
    3. Assumed unit weights of trucks
        List of float types

    Give the paths to the input data files:
    1. Network edges EXcel File
    2. OD flows Excel file
    3. Road properties Excel file
    
    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    province_list = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
    province_terrian = ['mountain', 'flat', 'flat']
    truck_unit_wt = [5.0]

    """Give the paths to the input data files
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    network_data_excel = os.path.join(data_path,'post_processed_networks','province_roads_edges.xlsx')
    od_output_excel = os.path.join(
        output_path, 'flow_ods','province_roads_commune_center_flow_ods.xlsx')
    rd_prop_file = os.path.join(data_path, 'mode_properties', 'road_properties.xlsx')
    
    """Specify the output files and paths to be created 
    """
    output_dir = os.path.join(output_path, 'flow_mapping_shapefiles')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    output_dir = os.path.join(output_path, 'flow_mapping_paths')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    flow_output_excel = os.path.join(
        output_dir, 'province_roads_commune_center_access_flow_paths.xlsx')
    excl_wrtr = pd.ExcelWriter(flow_output_excel)

    """Start the OD flow mapping process    
    """
    for prn in range(len(province_list)):
        province = province_list[prn]
        province_name = province.replace(' ', '').lower()

        """Load igraph network and GeoDataFrame
        """
        print ('* Loading {} igraph network and GeoDataFrame'.format(province))
        edges_in = pd.read_excel(network_data_excel,sheet_name = province_name,encoding='utf-8')
        G = ig.Graph.TupleList(edges_in.itertuples(index=False), edge_attrs=list(edges_in.columns)[2:])
        del edges_in
        gdf_edges = gpd.read_file(os.path.join(network_data_path,'{}_roads_edges.shp'.format(province_name)),encoding='utf-8')
        gdf_edges = gdf_edges[['edge_id','geometry']]
        
        """Load OD nodes pairs and tonnages
        """
        print ('* Loading {} OD nodes pairs and tonnages'.format(province))
        all_ods = pd.read_excel(od_output_excel, sheet_name=province_name)
        all_ods = all_ods[['origin', 'destination', 'min_croptons',
                           'max_croptons', 'min_netrev', 'max_netrev']]
        
        for tr_wt in truck_unit_wt:
            """Calculate OD paths
            """
            print ('* Calculating {} OD paths'.format(province))
            all_ods['min_vehicle_nums'] = np.maximum(1, np.ceil(all_ods['min_croptons']/tr_wt))
            all_ods['max_vehicle_nums'] = np.maximum(1, np.ceil(all_ods['max_croptons']/tr_wt))
            G = add_igraph_generalised_costs_roads(G, 1, tr_wt)
            all_paths = network_od_paths_assembly(all_ods, G, tr_wt,province_name, excel_writer=excl_wrtr)
            
            """Create network shapefiles with flows
            """
            print ('* Creating {} network shapefiles with flows'.format(province))
            shp_output_path = os.path.join(
                output_path,'flow_mapping_shapefiles',
                'weighted_edges_commune_center_access_flows_{0}_{1}_tons.shp'.format(province_name,int(tr_wt)))

            write_province_flow_paths_to_network_shapefile(
                all_paths, gdf_edges, save_edges=True, shape_output_path=shp_output_path)


