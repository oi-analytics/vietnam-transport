"""
Purpose
-------

Mapping the commune access OD node level matrix values to road network paths in Provinces

For all roads in the Provinces: ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']

The code estimates 2 values - A MIN and a MAX value of flows between each selected OD node pair
    - Based on MIN-MAX generalised costs estimates

Input data requirements
-----------------------
1. Correct paths to all files and correct input parameters
2. Excel file with mode sheets containing network graph structure and attributes
    - edge_id - String Edge ID
    - from_node - String node ID that should be present in node_id column
    - to_node - String node ID that should be present in node_id column
    - length - Float length of edge in km
    - min_time - Float minimum time of travel in hours on edge 
    - max_time - Float maximum time of travel in hours on edge  
    - min_time_cost - Float minimum cost of time in USD on edge   
    - max_time_cost - Float maximum cost of time in USD on edge 
    - min_tariff_cost - Float minimum tariff cost in USD on edge   
    - max_tariff_cost - Float maximum tariff cost in USD on edge
         
3. Edge shapefiles for all national-scale networks with attributes:
    - edge_id - String Edge ID
    - geometry - Shapely LineString geometry of edges

4. Excel file with mode sheets containing node-level OD values with attributes:
    - origin - String node ID of Origin
    - destination - String node ID of Destination
    - min_netrev - Float values of miniimum daily OD Net Revenue in USD
    - max_netrev - Float values of maximum daily OD Net Revenue in USD
    - min_tons -  Float values of minimum daily OD in tons
    - max_tons - Float values of maximum daily OD in tons

Results
-------
1. Excel sheets with results of flow mapping based on MIN-MAX generalised costs estimates:
    - origin - String node ID of Origin
    - destination - String node ID of Destination
    - min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
    - max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
    - min_netrev - Float values of estimated daily Net Revenue for paths with minimum generalised cost flows
    - max_netrev - Float values of estimated daily Net Revenue for paths with maximum generalised cost flows
    - min_croptons - Float values of estimated daily crop tonnage for paths with minimum generalised cost flows
    - max_croptons - Float values of estimated daily crop tonnage for paths with maximum generalised cost flows
    - min_distance - Float values of estimated distance for paths with minimum generalised cost flows
    - max_distance - Float values of estimated distance for paths with maximum generalised cost flows
    - min_time - Float values of estimated time for paths with minimum generalised cost flows
    - max_time - Float values of estimated time for paths with maximum generalised cost flows
    - min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
    - max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
    - min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
    - max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows

2. Shapefiles with all flows on edges mapping based on MIN-MAX generalised costs estimates:
    - edge_id - String/Integer/Float Edge ID
    - geometry - Shapely LineString geomtry of edges
    - min_netrev - Float values of estimated daily Net Revenue in USD on edges
    - max_netrev - Float values of estimated daily Net Revenue in USD on edges
    - min_tons - Float values of estimated daily crops in tons on edges
    - max_tons - Float values of estimated daily crops in tons on edges

References
----------
1. Pant, R., Koks, E.E., Russell, T., Schoenmakers, R. & Hall, J.W. (2018).
   Analysis and development of model for addressing climate change/disaster risks in multi-modal transport networks in Vietnam.
   Final Report, Oxford Infrastructure Analytics Ltd., Oxford, UK.
2. All input data folders and files referred to in the code below.
"""

import ast
import itertools
import math
import operator
import os
import subprocess
import sys
import copy

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd
from shapely import wkt
from shapely.geometry import Point
from vtra.transport_flow_and_failure_functions import *
from vtra.utils import *


def network_od_paths_assembly_provincial(points_dataframe, graph, 
    vehicle_wt, region_name,excel_writer=''):
    """
    Assemble estimates of OD paths, distances, times, costs and tonnages on networks

    Parameters
    ---------
    - points_dataframe - Pandas DataFrame of OD nodes and their tonnages
    - graph - igraph network structure 
    - vehicle_wt - Float unit weight of vehicle
    - region_name - String name of Province
    - excel_writer - Name of the excel writer to save Pandas dataframe to Excel file     

    Outputs
    -------
    save_paths_df - Pandas DataFrame with columns:
        - origin - String node ID of Origin
        - destination - String node ID of Destination
        - min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
        - max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
        - min_netrev - Float values of estimated netrevenue for paths with minimum generalised cost flows
        - max_netrev - Float values of estimated netrevenue for paths with maximum generalised cost flows
        - min_croptons - Float values of estimated crop tons for paths with minimum generalised cost flows
        - max_croptons - Float values of estimated crop tons for paths with maximum generalised cost flows
        - min_distance - Float values of estimated distance for paths with minimum generalised cost flows
        - max_distance - Float values of estimated distance for paths with maximum generalised cost flows
        - min_time - Float values of estimated time for paths with minimum generalised cost flows
        - max_time - Float values of estimated time for paths with maximum generalised cost flows
        - min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
        - max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
        - min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
        - max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows

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
                graph, origin, destinations,vehicle_wt,vehicle_wt, 'min_gcost', 'min_time')
            get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(
                graph, origin, destinations,vehicle_wt,vehicle_wt, 'max_gcost', 'max_time')


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
    del save_paths

    return save_paths_df

if __name__ == '__main__':
    """
    1. Specify the paths from where you want to read and write:
        - Input data
        - Intermediate calcuations data
        - Output results
    
    2. Supply input data and parameters
        - Names of the three Provinces: List of string types 
        - Assumed terrains of the provinces: List of string types
        - Assumed unit weights of trucks: List of float types

    3. Give the paths to the input data files:
        - Network edges EXcel File
        - OD flows Excel file
        - Road properties Excel file
    
    4. Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    province_list = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
    province_terrian = ['mountain', 'flat', 'flat']
    truck_unit_wt = [5.0]
    percentage = [100.0]

    """Give the paths to the input data files
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    network_data_excel = os.path.join(data_path,'post_processed_networks','province_roads_edges.xlsx')
    od_output_excel = os.path.join(
        output_path, 'flow_ods','province_roads_commune_center_flow_ods.xlsx')
    rd_prop_file = os.path.join(data_path, 'mode_properties', 'road_properties.xlsx')
    
    """Specify the output files and paths to be created 
    """
    flow_shp_dir = os.path.join(output_path, 'flow_mapping_shapefiles')
    if os.path.exists(flow_shp_dir) == False:
        os.mkdir(flow_shp_dir)

    flow_csv_dir = os.path.join(output_path, 'flow_mapping_combined')
    if os.path.exists(flow_csv_dir) == False:
        os.mkdir(flow_csv_dir)

    flow_paths_dir = os.path.join(output_path, 'flow_mapping_paths')
    if os.path.exists(flow_paths_dir) == False:
        os.mkdir(flow_paths_dir)

    for perct in percentage:
        flow_output_excel = os.path.join(
            flow_paths_dir, 'province_roads_commune_center_access_flow_paths_{}_percent.xlsx'.format(int(perct)))
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
            ods = pd.read_excel(od_output_excel, sheet_name=province_name)
            ods = ods[['origin', 'destination', 'min_croptons',
                               'max_croptons', 'min_netrev', 'max_netrev']]
            
            all_ods = copy.deepcopy(ods)
            all_ods_tons_cols = [col for col in all_ods.columns.values.tolist() if col not in ['origin','destination']] 
            all_ods[all_ods_tons_cols] = 0.01*perct*all_ods[all_ods_tons_cols]

            for tr_wt in truck_unit_wt:
                """Calculate OD paths
                """
                print ('* Calculating {} OD paths'.format(province))
                all_ods['min_vehicle_nums'] = np.maximum(1, np.ceil(all_ods['min_croptons']/tr_wt))
                all_ods['max_vehicle_nums'] = np.maximum(1, np.ceil(all_ods['max_croptons']/tr_wt))
                # G = add_igraph_generalised_costs_roads(G, 1, tr_wt)
                all_paths = network_od_paths_assembly_provincial(all_ods, G, tr_wt,province_name, excel_writer=excl_wrtr)
                
                """Create network shapefiles with flows
                """
                print ('* Creating {} network csv and shapefiles with flows'.format(province))
                shp_output_path = os.path.join(
                    flow_shp_dir,
                    'weighted_edges_commune_center_access_flows_{0}_{1}_tons_{2}_percent.shp'.format(province_name,int(tr_wt),int(perct)))
                csv_output_path = os.path.join(
                    flow_csv_dir,
                    'weighted_edges_commune_center_access_flows_{0}_{1}_tons_{2}_percent.csv'.format(province_name,int(tr_wt),int(perct)))

                write_flow_paths_to_network_files(
                    all_paths,['netrev','croptons'],['netrev','croptons'], gdf_edges, 
                    save_csv=True, save_shapes=False, shape_output_path=shp_output_path,csv_output_path=csv_output_path)

            del all_ods

