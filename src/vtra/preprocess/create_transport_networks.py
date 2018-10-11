"""
Purpose
-------

Creating post-processed transport networks with attributes 

From pre-processed input Shapefiles and collected network attributes data  

For all Province road networks: ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']

For all transport modes at national scale: ['road', 'rail', 'air', 'inland', 'coastal','multi']

Input data requirements
-----------------------

Correct paths to all files and correct input parameters

Edge and node shapefiles for all Province roads and national-scale networks

All Geometries in Edge Shapefiles should be valid LineStrings

All Geometries in Node Shapefiles should be valid Points
	
A. Node Shapefiles should contain following column names and attributes
    - node_id - String node ID
    - geometry - Shapely Point geometry of nodes
    - attributes - Multiple types depending upon sector and context 

B. Edge Shapefiles should contain following column names and attributes:
    - edge_id - String Edge ID
    - from_node - String node ID that should be present in node_id column
    - to_node - String node ID that should be present in node_id column 
    - geometry - Shapely LineString geometry of edges 


Results
-------
1. Excel sheets with post-processed network nodes and edges 
2. Shapefiles with post-processed network nodes and edges
	
3. All nodes have the following attributes:
	- node_id - String Node ID
	- name - String name in Vietnamese/English
	- tons - Float assigned cargo freight tonnage using node 
	- population - Float assigned passenger/population number using node 
	- capacity - Float assigned capacity in tons/passenger numbers/other units
	- geometry - Shapely Point geometry of node

4. Attributes only present in inland and coastal port nodes
	- port_type - String name of type of port: inland or sea 	
	- port_class - String name of class of port: class1A (international) or class1 (domestic hub)  

5. All edges have the following attributes:
	- edge_id - String edge ID
	- g_id - Interger edge ID
	- from_node - String node ID that should be present in node_id column
	- to_node - String node ID that should be present in node_id column
	- geometry - Shapely LineString geometry of edge
	- terrain - String name of terrain of edge	
	- level - Integer number for edge level: National, Provincial, Local, etc.
	- width - Float width in meters of edge
	- length - Float estimated length in kilometers of edge	
	- min_speed - Float estimated minimum speed in km/hr on edge
	- max_speed - Float estimated maximum speed in km/hr on edge
	- min_time - Float estimated minimum time of travel in hours on edge
	- max_time - Float estimated maximum time of travel in hours on edge	
	- min_time_cost - Float estimated minimum cost of time in USD on edge
	- max_time_cost - Float estimated maximum cost of time in USD on edge
	- min_tariff_cost - Float estimated minimum tariff cost in USD on edge	
	- max_tariff_cost - Float estimated maximum tariff cost in USD on edge
	- vehicle_co - Integer number of daily vehicle counts on edge

6. Attributes only present in Province and national roads edges
	- surface - String value for surface
	- road_class - Integer between 1 and 6
	- road_cond - String value: paved or unpaved 
	- asset_type - String name of type of asset

References
----------
1. Pant, R., Koks, E.E., Russell, T., Schoenmakers, R. & Hall, J.W. (2018).
   Analysis and development of model for addressing climate change/disaster risks in multi-modal transport networks in Vietnam.
   Final Report, Oxford Infrastructure Analytics Ltd., Oxford, UK.
2. All input data folders and files referred to in the code below.		 
"""
import os
import subprocess
import sys

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd
from shapely import wkt,wkb
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from vtra.preprocess.transport_network_inputs import *
from vtra.utils import *


def main():
    """
    Specify the paths from where you want to read and write:
    
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    
    1. Paths of the mode files
        - List of tuples of strings 
    2. Names of modes
        - List of strings
    3. Unit weight of vehicle assumed for each mode
        - List of float types 
    4. Range of usage factors for each mode to represent uncertainty in cost estimations
        - List of tuples of float types 
    5. Ranges of speeds for each mode to represent uncertainty in speeds
        - List of tuple of float types 
    6. Names of all industry sector and crops in VITRANSS2 and IFPRI datasets
        - List of string types
    7. Names of commodity/industry columns for which min-max tonnage column names already exist
        - List of string types
    8. Percentage of OD flow we want to send along path
        - Float type

    Give the paths to the input data files:
    
    1. Pre-proccesed network shapefiles
    2. Costs of modes Excel file 
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
    modes_file_paths = [('Roads', 'national_roads'), ('Railways', 'national_rail'),
                        ('Airports', 'airnetwork'), ('Waterways', 'waterways'), 
                        ('Waterways', 'waterways'),('Multi', 'multi_edges')]
    modes = ['road', 'rail', 'air', 'inland', 'coastal','multi']
    veh_wt = [20, 800, 0, 800, 1200]
    cost_uncertainty_factors = [(0, 0), (0, 0), (0, 0), (0.2, 0.25), (0.2, 0.25), (0,0)]
    speeds = [(0, 0), (40, 60), (700, 800), (9, 20), (9, 20), (0,0)]
    multi_md_len = 3.0
    province_results = 'Yes'
    national_results = 'Yes'
    
    """Give the paths to the input data files
    """
    network_data_path = os.path.join(
        data_path, 'pre_processed_networks_data')
    md_prop_file = os.path.join(network_data_path, 'mode_properties', 'mode_costs.xlsx')
    rd_prop_file = os.path.join(network_data_path, 'mode_properties', 'road_properties.xlsx')


    """Specify the output files and paths to be created 
    """
    output_dir = os.path.join(data_path, 'post_processed_networks')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    

    """Start the OD flow mapping process    
    """
    if province_results == 'Yes':
        edges_excel = os.path.join(
            output_dir, 'province_roads_edges.xlsx')
        edges_excl_wrtr = pd.ExcelWriter(edges_excel)

        nodes_excel = os.path.join(
            output_dir, 'province_roads_nodes.xlsx')
        nodes_excl_wrtr = pd.ExcelWriter(nodes_excel)        
        for prn in range(len(province_list)):
            province = province_list[prn]
            province_name = province.replace(' ', '').lower()

            """Load igraph network and GeoDataFrame
            """
            print ('* Creating and writing {} GeoDataFrame to Shapefiles and Excel Sheets'.format(province))
            province_path = os.path.join(network_data_path, 'Roads', '{}_roads'.format(province_name))
            nodes_in,edges_in = get_node_edge_files_in_path(province_path)
            gdf_edges = province_shapefile_to_dataframe(
                edges_in, province_terrian[prn], rd_prop_file,(0,0))
            
            gdf_nodes = read_setor_nodes(nodes_in,province_name)
            
            gdf_edges.to_file(os.path.join(output_dir, '{}_roads_edges.shp'.format(province_name)),encoding = 'utf-8')
            gdf_edges.drop('geometry', axis=1, inplace=True)
            gdf_edges.to_excel(edges_excl_wrtr,province_name, index=False)
            edges_excl_wrtr.save()
            
            gdf_nodes.to_file(os.path.join(output_dir, '{}_roads_nodes.shp'.format(province_name)),encoding = 'utf-8')
            gdf_nodes.drop('geometry', axis=1, inplace=True)
            gdf_nodes.to_excel(nodes_excl_wrtr,province_name, index=False)
            nodes_excl_wrtr.save()


    """Start the OD flow mapping process    
    """
    if national_results == 'Yes':
        edges_excel = os.path.join(
            output_dir, 'national_edges.xlsx')
        edges_excl_wrtr = pd.ExcelWriter(edges_excel)

        nodes_excel = os.path.join(
            output_dir, 'national_nodes.xlsx')
        nodes_excl_wrtr = pd.ExcelWriter(nodes_excel) 
        for m in range(len(modes)):
            print ('* Creating and writing {} GeoDataFrame to Shapefiles and Excel Sheets'.format(modes[m]))
            mode_data_path = os.path.join(
                network_data_path, modes_file_paths[m][0], modes_file_paths[m][1])
            nodes_in,edges_in = get_node_edge_files_in_path(mode_data_path)

            """Load mode igraph network and GeoDataFrame
            """
            if modes[m] == 'road':
                gdf_edges = national_road_shapefile_to_dataframe(edges_in, rd_prop_file,cost_uncertainty_factors[m])
            elif modes[m] == 'multi':
                gdf_edges = multi_modal_shapefile_to_dataframe(edges_in, md_prop_file, modes[m], multi_md_len,cost_uncertainty_factors[m])
            else:
                gdf_edges = network_shapefile_to_dataframe(
                    edges_in, md_prop_file, modes[m], speeds[m][0], speeds[m][1],cost_uncertainty_factors[m])

            if modes[m] in ['inland','coastal']:
                port_names = os.path.join(mode_data_path,'ports.shp')
                gdf_nodes = read_waterway_ports(nodes_in,port_names)
            else:
                gdf_nodes = read_setor_nodes(nodes_in,modes[m])
            
            gdf_edges.to_file(os.path.join(output_dir, '{}_edges.shp'.format(modes[m])),encoding = 'utf-8')
            gdf_edges.drop('geometry', axis=1, inplace=True)
            gdf_edges.to_excel(edges_excl_wrtr,modes[m], index=False)
            edges_excl_wrtr.save()
            
            gdf_nodes.to_file(os.path.join(output_dir, '{}_nodes.shp'.format(modes[m])),encoding = 'utf-8')
            gdf_nodes.drop('geometry', axis=1, inplace=True)
            gdf_nodes.to_excel(nodes_excl_wrtr,modes[m], index=False)
            nodes_excl_wrtr.save()

if __name__ == '__main__':
    main()
