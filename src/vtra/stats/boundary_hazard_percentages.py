"""Summarise network-hazard intersections per-boundary (district, commune or province)

Purpose
-------

Collect network-hazard intersection attributes
    - Combine with boundary Polygons to collect network-boundary intersection attributes
    - Write final results to an Excel sheet

Input data requirements
-----------------------

1. Correct paths to all files and correct input parameters

2. Shapefiles of network-hazard intersections results with attributes:
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
    """Summarise intersections

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
        - Hazard datasets description Excel file
        - String name of sheet in hazard datasets description Excel file

    4. Specify the output files and paths to be created
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    # Supply input data and parameters
    provinces = ['Lao Cai','Binh Dinh','Thanh Hoa']
    modes = ['road','rail','air','inland','coastal']
    out_modes = ['national_roads', 'national_rail', 'air_ports', 'inland_ports', 'sea_ports']
    boundary_cols = ['commune_id','commune_name','district_id','district_name','province_id','province_name']
    hazard_cols = ['climate_scenario','hazard_type','model','probability','year']
    province_results = 'No'
    national_results = 'Yes'

    # Give the paths to the input data files
    province_file = os.path.join(output_path,
            'network_stats',
            'province_roads_stats.xlsx')

    national_file = os.path.join(output_path,
            'network_stats',
            'national_scale_stats.xlsx')

    province_hazard_file = os.path.join(output_path,
            'hazard_scenarios',
            'province_roads_hazard_intersections.xlsx')

    national_hazard_file = os.path.join(output_path,
            'hazard_scenarios',
            'national_scale_hazard_intersections.xlsx')
    # Specify the output files and paths to be created
    output_dir = os.path.join(output_path, 'network_stats')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    # Process province scale results
    if province_results == 'Yes':
        print ('* Processing province scale results')
        data_excel = os.path.join(
            output_dir,'province_roads_hazards_stats.xlsx')
        prov_excel_writer = pd.ExcelWriter(data_excel)
        for province in provinces:
            province_name = province.replace(' ', '').lower()
            province_roads_stats = pd.read_excel(province_file,sheet_name=province_name)
            province_roads_stats = province_roads_stats.groupby(boundary_cols)['length'].sum().reset_index()
            province_roads_stats.rename(columns={'length':'total_length'},inplace=True)

            hazard_stats = pd.read_excel(province_hazard_file,sheet_name=province_name)
            hazard_stats = hazard_stats.groupby(boundary_cols+hazard_cols)['length'].sum().reset_index()
            hazard_stats = pd.merge(hazard_stats,province_roads_stats,how='left', on=boundary_cols).fillna(0)
            hazard_stats['percentage'] = 100.0*hazard_stats['length']/hazard_stats['total_length']

            hazard_stats.to_excel(prov_excel_writer, province_name, index=False)
            prov_excel_writer.save()
            del hazard_stats

    # Process national scale results
    if national_results == 'Yes':
        print ('* Processing national scale results')
        data_excel = os.path.join(
            output_dir,'national_scale_hazards_stats.xlsx')
        nat_excel_writer = pd.ExcelWriter(data_excel)
        for m in range(len(modes)):
            if modes[m] in ['road','rail']:
                national_edges_stats = pd.read_excel(national_file,sheet_name=modes[m])
                national_edges_stats = national_edges_stats.groupby(boundary_cols)['length'].sum().reset_index()
                national_edges_stats.rename(columns={'length':'total_length'},inplace=True)

                hazard_stats = pd.read_excel(national_hazard_file,sheet_name=modes[m])
                hazard_stats = hazard_stats.groupby(boundary_cols+hazard_cols)['length'].sum().reset_index()
                hazard_stats = pd.merge(hazard_stats,national_edges_stats,how='left', on=boundary_cols).fillna(0)
                hazard_stats['percentage'] = 100.0*hazard_stats['length']/hazard_stats['total_length']

                hazard_stats.to_excel(nat_excel_writer, modes[m], index=False)
                nat_excel_writer.save()
                del hazard_stats


if __name__ == "__main__":
    main()
