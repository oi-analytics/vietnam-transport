"""Mapping the OD node level matrix values to network paths
For all transport modes at national scale:
    ['road', 'rail', 'air', 'inland', 'coastal']

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
            min_tons -  Float values of minimum daily OD in tons
            max_tons - Float values of maximum daily OD in tons
        Should also contain names of the industry columns specified in the inputs
            
Results
-------
1. Excel sheets with results of flow mapping based on MIN-MAX generalised costs estimates:
        origin - String node ID of Origin
        destination - String node ID of Destination
        o_region - String name of Province of Origin node ID
        d_region - String name of Province of Destination node ID
        min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
        max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
        min_distance - Float values of estimated distance for paths with minimum generalised cost flows
        max_distance - Float values of estimated distance for paths with maximum generalised cost flows
        min_time - Float values of estimated time for paths with minimum generalised cost flows
        max_time - Float values of estimated time for paths with maximum generalised cost flows
        min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
        max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
        min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
        max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows
        industry_columns - All daily tonnages of industry columns given in the OD matrix data
2. Shapefiles
        edge_id - String/Integer/Float Edge ID
        geometry - Shapely LineString geomtry of edges
        min_{industry} - Float values of estimated minimum daily industries/commodities/total volumes in tons on edges 
        max_{industry} - Float values of estimated maximum daily industries/commodities/total volumes in tons on edges
            For all industry/commodities in the list: 
            ['cement', 'coal', 'constructi', 'fertilizer', 'fishery', 
            'manufactur', 'acof', 'cash', 'cass', 'maiz', 'pepp', 'rcof',
            'rubb', 'swpo', 'teas', 'meat', 'rice', 'petroluem', 'steel', 
            'sugar', 'wood', 'tons'] 

"""
import os
import subprocess
import sys
import copy

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd
from shapely import wkt
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from vtra.transport_flow_and_failure_functions import *
from vtra.utils import *


def network_od_paths_assembly_national(points_dataframe, 
    graph, vehicle_wt, transport_mode, 
    excel_writer=''):
    """
    Assembles estimates of OD paths, distances, times, costs and tonnages on networks 

    Parameters
    ---------
    points_dataframe - Pandas DataFrame of OD nodes and their tonnages
    graph - igraph network structure 
    vehicle_wt - Float unit weight of vehicle
    transport_mode - String name of modes
    excel_writer - Name of the excel writer to save Pandas dataframe to Excel file     

    Outputs
    -------
    save_paths_df - Pandas DataFrame 
        With columns:
            origin - String node ID of Origin
            destination - String node ID of Destination
            min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
            max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
            min_distance - Float values of estimated distance for paths with minimum generalised cost flows
            max_distance - Float values of estimated distance for paths with maximum generalised cost flows
            min_time - Float values of estimated time for paths with minimum generalised cost flows
            max_time - Float values of estimated time for paths with maximum generalised cost flows
            min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
            max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
            min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
            max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows
            industry_columns - All tonnages of industry columns given in the OD matrix data

    """
    save_paths = []
    for iter_, row in points_dataframe.iterrows():
        try:
            origin = row['origin']
            destinations = [row['destination']]
            tons = row['min_tons']
            get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations(
                graph, origin, destinations, tons, vehicle_wt, 'min_gcost', 'min_time')
            tons = row['max_tons']
            get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations(
                graph, origin, destinations, tons, vehicle_wt, 'max_gcost', 'max_time')
            save_paths += list(zip([origin]*len(destinations), destinations, get_min_path, get_max_path,
                                   get_min_dist, get_max_dist, get_min_time, get_max_time, get_min_gcost, get_max_gcost))
            print("done with {} in network {}".format(origin, transport_mode))
        except:
            print('* no path between {}-{}'.format(origin,destinations))

    cols = ['origin', 'destination', 'min_edge_path', 'max_edge_path',
            'min_distance', 'max_distance', 'min_time', 'max_time', 'min_gcost', 'max_gcost']
    save_paths_df = pd.DataFrame(save_paths, columns=cols)

    points_dataframe = points_dataframe.reset_index()
    save_paths_df = pd.merge(save_paths_df, points_dataframe, how='left', on=[
                             'origin', 'destination']).fillna(0)

    save_paths_df = save_paths_df[(save_paths_df['max_tons'] > 0)
                                  & (save_paths_df['origin'] != 0)]
    if transport_mode != 'air':
        save_paths_df['min_vehicle_nums'] = np.maximum(
            1, np.ceil(save_paths_df['min_tons']/vehicle_wt))
        save_paths_df['max_vehicle_nums'] = np.maximum(
            1, np.ceil(save_paths_df['max_tons']/vehicle_wt))

    save_paths_df.to_excel(excel_writer, transport_mode, index=False)
    excel_writer.save()
    del save_paths
    
    return save_paths_df

def main():
    """
    Specify the paths from where you want to read and write:
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    1. Names of modes
        List of strings
    2. Unit weight of vehicle assumed for each mode
        List of float types 
    3. Names of all industry sector and crops in VITRANSS2 and IFPRI datasets
        List of string types
    4. Names of commodity/industry columns for which min-max tonnage column names already exist
        List of string types
    5. Percentage of OD flow we want to send along path
        FLoat type

    Give the paths to the input data files:
    1. Network edges Excel file
    2. OD flows Excel file
    3. Costs of modes Excel file 
    4. Road properties Excel file
    
    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    modes = ['road', 'rail', 'air', 'inland', 'coastal']
    veh_wt = [20, 800, 0, 800, 1200]
    ind_cols = ['cement', 'coal', 'constructi', 'fertilizer', 'fishery', 'manufactur', 'acof', 'cash', 'cass', 'maiz', 'pepp', 'rcof',
                'rubb', 'swpo', 'teas', 'meat', 'rice', 'petroluem', 'steel', 'sugar', 'wood', 'tons']
    min_max_exist = ['rice','tons']
    percentage = [10,90,100]
    
    """Give the paths to the input data files
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    network_data_excel = os.path.join(data_path,'post_processed_networks','national_edges.xlsx')
    od_output_excel = os.path.join(
        output_path, 'flow_ods', 'national_scale_flow_ods.xlsx')
    md_prop_file = os.path.join(data_path, 'mode_properties', 'mode_costs.xlsx')
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
            flow_paths_dir, 
            'national_scale_flow_paths_{}_percent.xlsx'.format(int(perct)))
    
        excl_wrtr = pd.ExcelWriter(flow_output_excel)

        """Start the OD flow mapping process    
        """
        for m in range(len(modes)):
            """Load mode igraph network and GeoDataFrame
            """
            print ('* Loading {} igraph network and GeoDataFrame'.format(modes[m]))
            edges_in = pd.read_excel(network_data_excel,sheet_name = modes[m],encoding='utf-8')
            G = ig.Graph.TupleList(edges_in.itertuples(index=False), edge_attrs=list(edges_in.columns)[2:])
            del edges_in
            gdf_edges = gpd.read_file(os.path.join(network_data_path,'{}_edges.shp'.format(modes[m])),encoding='utf-8')
            gdf_edges = gdf_edges[['edge_id','geometry']]

            """Load mode OD nodes pairs and tonnages
            """
            print ('* Loading {} OD nodes pairs and tonnages'.format(modes[m]))
            ods = pd.read_excel(od_output_excel, sheet_name=modes[m])
            ods = ods[ods['max_tons'] > 0.5]
            
            all_ods = copy.deepcopy(ods)
            all_ods_tons_cols = [col for col in all_ods.columns.values.tolist() if col not in ['origin','o_region','destination','d_region']] 
            all_ods[all_ods_tons_cols] = 0.01*perct*all_ods[all_ods_tons_cols]
            """Calculate mode OD paths
            """
            print ('* Calculating {} OD paths'.format(modes[m]))
            all_paths = network_od_paths_assembly_national(
                all_ods, G, veh_wt[m], modes[m],excel_writer=excl_wrtr)

            del all_ods
            """Create network shapefiles with flows
            """
            print ('* Creating {} network shapefiles with flows'.format(modes[m]))
            
            shp_output_path = os.path.join(flow_shp_dir,'weighted_flows_national_{}_{}_percent.shp'.format(modes[m],int(perct)))
            csv_output_path = os.path.join(flow_csv_dir,'weighted_flows_national_{}_{}_percent.csv'.format(modes[m],int(perct)))


            write_flow_paths_to_network_files(all_paths,
                ind_cols,min_max_exist,gdf_edges, 
                save_csv=True, save_shapes=False, shape_output_path=shp_output_path,csv_output_path=csv_output_path)


if __name__ == '__main__':
    main()
