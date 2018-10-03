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

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd
from shapely import wkt
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from vtra.transport_network_creation import *
from vtra.utils import *

def network_od_path_estimations_changing_tonnages(graph,
    source, target, tonnage, vehicle_wt, utilization_factors, cost_criteria, time_criteria):
    """
    Estimate the paths, distances, times, and costs for given OD pair

    Parameters
    ---------
    graph - igraph network structure 
    source - String/Float/Integer name of Origin node ID
    source - String/Float/Integer name of Destination node ID
    tonnage - Float value of tonnage 
    vehicle_wt - Float unit weight of vehicle
    utilization_factors - Tuple of float types for uncertainity in cost function
    cost_criteria - String name of generalised cost criteria to be used: min_gcost or max_gcost
    time_criteria - String name of time criteria to be used: min_time or max_time    

    Outputs
    -------
    edge_path_list - List of lists of Strings/Floats/Integers of edge ID's in routes
    path_dist_list - List of float values of estimated distances of routes
    path_time_list - List of float values of estimated times of routes
    path_gcost_list - List of float values of estimated generalised costs of routes

    """

    graph = add_igraph_generalised_costs_network(graph, np.ceil(
        tonnage/vehicle_wt), tonnage, utilization_factors[0], utilization_factors[1])
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


def network_od_paths_assembly_changing_tonnages(points_dataframe, 
    graph, vehicle_wt, utilization_factors, transport_mode, 
    percentage = 100.0, excel_writer=''):
    """
    Assembles estimates of OD paths, distances, times, costs and tonnages on networks 

    Parameters
    ---------
    points_dataframe - Pandas DataFrame of OD nodes and their tonnages
    graph - igraph network structure 
    vehicle_wt - Float unit weight of vehicle
    utilization_factors - Tuple of float types for uncertainity in cost function
    transport_mode - String name of modes
    percentage - Float value of the percentage of OD flow we want to allocate on paths
        Default = 100.0 
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
            tons = 1.0*percentage/100.0*row['min_tons']
            get_min_path, get_min_dist, get_min_time, get_min_gcost = network_od_path_estimations_changing_tonnages(
                graph, origin, destinations, tons, vehicle_wt, utilization_factors, 'min_gcost', 'min_time')
            tons = 1.0*percentage/100.0*row['max_tons']
            get_max_path, get_max_dist, get_max_time, get_max_gcost = network_od_path_estimations_changing_tonnages(
                graph, origin, destinations, tons, vehicle_wt, utilization_factors, 'max_gcost', 'max_time')
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

def write_national_flow_paths_to_network_shapefile(save_paths_df,transport_mode,
    industry_columns,min_max_exist,gdf_edges, save_edges=True, shape_output_path=''):
    """
    Write results to Shapefiles

    Parameters
    ---------
    save_paths_df - Pandas DataFrame of OD flow paths and their tonnages
    transport_mode - String name of modes
    industry_columns - List of string names of all OD commodities/industries indentified
    min_max_exist - List of string names of commodity/industry columns for which min-max tonnage column names already exist
    gdf_edges - GeoDataFrame of network edge set
    save_Edges - Boolean condition to tell code to save created edge shapefile
    shape_output_path - Path where the output shapefile will be stored 

    Outputs
    -------
    gdf_edges - Shapefile 
        With minimum and maximum tonnage flows of all commodities/industries for each edge of network
    """

    min_ind_cols = []
    max_ind_cols = []
    ch_min_ind_cols = []
    ch_max_ind_cols = []
    for ind in industry_columns:
        min_ind_cols.append('min_{}'.format(ind))
        max_ind_cols.append('max_{}'.format(ind))
        if ind in min_max_exist:
            ch_min_ind_cols.append('min_{}'.format(ind))
            ch_max_ind_cols.append('max_{}'.format(ind))
        else:
            ch_min_ind_cols.append(ind)
            ch_max_ind_cols.append(ind)

    for i in range(len(min_ind_cols)):
        gdf_edges[min_ind_cols[i]] = 0
        gdf_edges[max_ind_cols[i]] = 0

    for iter_, path in save_paths_df.iterrows():
        min_path = path['min_edge_path']
        max_path = path['max_edge_path']

        gdf_edges.loc[gdf_edges['edge_id'].isin(min_path), min_ind_cols] += path[ch_min_ind_cols].values
        gdf_edges.loc[gdf_edges['edge_id'].isin(max_path), max_ind_cols] += path[ch_max_ind_cols].values


    for ind in industry_columns:
        gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)), axis = 1)
        gdf_edges[['min_{}'.format(ind),'max_{}'.format(ind)]] = gdf_edges['swap'].apply(pd.Series)
        gdf_edges.drop('swap', axis=1, inplace=True)

    if save_edges == True:
        gdf_edges.to_file(shape_output_path,encoding='utf-8')

    del gdf_edges, save_paths_df

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
    3. Range of usage factors for each mode to represent uncertainty in cost estimations
        List of tuples of float types 
    4. Names of all industry sector and crops in VITRANSS2 and IFPRI datasets
        List of string types
    5. Names of commodity/industry columns for which min-max tonnage column names already exist
        List of string types
    6. Percentage of OD flow we want to send along path
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
    usage_factors = [(0, 0), (0, 0), (0, 0), (0.2, 0.25), (0.2, 0.25)]
    ind_cols = ['cement', 'coal', 'constructi', 'fertilizer', 'fishery', 'manufactur', 'acof', 'cash', 'cass', 'maiz', 'pepp', 'rcof',
                'rubb', 'swpo', 'teas', 'meat', 'rice', 'petroluem', 'steel', 'sugar', 'wood', 'tons']
    min_max_exist = ['rice','tons']
    percentage = 100.0
    
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
    output_dir = os.path.join(output_path, 'flow_mapping_shapefiles')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    output_dir = os.path.join(output_path, 'flow_mapping_paths')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    if percentage == 100.0:
        flow_output_excel = os.path.join(
            output_dir, 'national_scale_flow_paths.xlsx')
    else:
        flow_output_excel = os.path.join(
            output_dir, 
            'national_scale_flow_paths_{}_percent.xlsx'.format(percentage))
    
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
        all_ods = pd.read_excel(od_output_excel, sheet_name=modes[m])
        all_ods = all_ods[all_ods['max_tons'] > 0.5]
        
        """Calculate mode OD paths
        """
        print ('* Calculating {} OD paths'.format(modes[m]))
        all_paths = network_od_paths_assembly_changing_tonnages(
            all_ods, G, veh_wt[m], usage_factors[m], modes[m], 
            percentage = percentage,excel_writer=excl_wrtr)

        """Create network shapefiles with flows
        """
        print ('* Creating {} network shapefiles with flows'.format(modes[m]))
        if percentage == 100.0:
            shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles','weighted_flows_national_{}.shp'.format(modes[m]))
        else:
            shp_output_path = os.path.join(output_path,'flow_mapping_shapefiles','weighted_flows_national_{}_{}_percent.shp'.format(modes[m],percentage))

        write_national_flow_paths_to_network_shapefile(all_paths,modes[m],
            ind_cols,min_max_exist,gdf_edges, save_edges=True, shape_output_path=shp_output_path)


if __name__ == '__main__':
    main()
