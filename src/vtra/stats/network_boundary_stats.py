"""Summarise length of edges/number of nodes within each boundary (commune, district, province)

Purpose
-------

Collect network attributes
    - Combine with boundary Polygons to collect network-boundary intersection attributes
    - Write final results to an Excel sheet

Input data requirements
-----------------------

1. Correct paths to all files and correct input parameters

2. Shapefiles of networks with attributes:
    - edge_id or node_id - String/Integer/Float Edge ID or Node ID of network
    - length - Float length of edge intersecting with hazards
    - geometry - Shapely geometry of edges as LineString or nodes as Points

3. Shapefile of administrative boundaries of Vietnam with attributes:
    - province_i - String/Integer ID of Province
    - pro_name_e - String name of Province in English
    - district_i - String/Integer ID of District
    - dis_name_e - String name of District in English
    - commune_id - String/Integer ID of Commune
    - name_eng - String name of Commune in English
    - geometry - Shapely geometry of boundary Polygon

Results
-------

1. Excel sheet of network-hazard-boundary intersection with attributes:
    - edge_id/node_id - String name of intersecting edge ID or node ID
    - length - Float length of intersection of edge LineString and hazard Polygon: Only for edges
    - province_id - String/Integer ID of Province
    - province_name - String name of Province in English
    - district_id - String/Integer ID of District
    - district_name - String name of District in English
    - commune_id - String/Integer ID of Commune
    - commune_name - String name of Commune in English
"""
import itertools
import os
import sys

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from vtra.utils import *
from vtra.transport_flow_and_failure_functions import *

def main():
    """Summarise

    1. Specify the paths from where you to read and write:
        - Input data
        - Intermediate calcuations data
        - Output results

    2. Supply input data and parameters
        - Names of the three Provinces - List of string types
        - Names of modes - List of strings
        - Names of output modes - List of strings
        - Names of hazard bands - List of integers
        - Names of hazard thresholds - List of integers
        - Condition 'Yes' or 'No' is the users wants to process results

    3. Give the paths to the input data files:
        - Commune boundary and stats data shapefile
        - String name of sheet in hazard datasets description Excel file

    4. Specify the output files and paths to be created
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    # Supply input data and parameters
    provinces = ['Lao Cai','Binh Dinh','Thanh Hoa']
    modes = ['road','rail','air','inland','coastal']
    out_modes = ['national_roads', 'national_rail', 'air_ports', 'inland_ports', 'sea_ports']
    province_results = 'Yes'
    national_results = 'Yes'

    # Give the paths to the input data files
    commune_shp = os.path.join(data_path, 'Vietnam_boundaries',
                                'who_boundaries', 'who_communes.shp')


    # Specify the output files and paths to be created
    output_dir = os.path.join(output_path, 'network_stats')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    # Process province scale results
    if province_results == 'Yes':
        print ('* Processing province scale results')
        data_excel = os.path.join(
            output_dir,'province_roads_stats.xlsx')
        prov_excel_writer = pd.ExcelWriter(data_excel)
        for province in provinces:
            province_name = province.replace(' ', '').lower()
            network_shp = os.path.join(
                data_path,
                'post_processed_networks',
                '{}_roads_edges.shp'.format(province_name))

            data_dict = []
            data_dict = spatial_scenario_selection(
                            network_shp, commune_shp, {}, data_dict,
                            network_type = 'edges',name_province = province)

            data_df = pd.DataFrame(data_dict)
            data_df.to_excel(prov_excel_writer, province_name, index=False)
            prov_excel_writer.save()
            del data_df

    # Process national scale results
    if national_results == 'Yes':
        print ('* Processing national scale results')
        data_excel = os.path.join(
            output_dir,'national_scale_stats.xlsx')
        nat_excel_writer = pd.ExcelWriter(data_excel)
        for m in range(len(modes)):
            if modes[m] in ['road','rail']:
                ntype = 'edges'
                network_shp = os.path.join(
                    data_path,
                    'post_processed_networks',
                    '{}_edges.shp'.format(modes[m]))
            else:
                ntype = 'nodes'
                network_shp = os.path.join(
                    data_path,
                    'post_processed_networks',
                    '{}_nodes.shp'.format(modes[m]))

            data_dict = []
            data_dict = spatial_scenario_selection(
                            network_shp, commune_shp, {}, data_dict,
                            network_type = ntype,name_province = '')
            data_df = pd.DataFrame(data_dict)

            data_df.to_excel(nat_excel_writer, modes[m], index=False)
            nat_excel_writer.save()
            del data_df




if __name__ == "__main__":
    main()
