"""
Purpose
-------

Collect network-hazard intersection attributes 
    - Combine with boundary Polygons to collect network-hazard-boundary intersection attributes 
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

4. Excel sheet of hazard attributes with attributes:
    - azard_type - String name of hazard type
    - odel - String name of hazard model
    - ear - String name of hazard year
    - limate_scenario - String name of hazard scenario
    - robability - Float/String value of hazard probability
    - and_num - Integer value of hazard band
    - in_val - Integer value of minimum value of hazard threshold
    - ax_val - Integer value of maximum value of hazard threshold

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
    - sector - String name of transport mode
    - hazard_type - String name of hazard type
    - model - String name of hazard model
    - year - String name of hazard year
    - climate_scenario - String name of hazard scenario
    - probability - Float/String value of hazard probability
    - band_num - Integer value of hazard band
    - min_val - Integer value of minimum value of hazard threshold
    - max_val - Integer value of maximum value of hazard threshold  
"""
import itertools
import os
import sys

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from vtra.utils import *


def spatial_scenario_selection(network_shapefile, polygon_shapefile, hazard_dictionary, data_dictionary, network_type ='nodes',name_province =''):
    """
    Intersect network edges/nodes and boundary Polygons to collect boundary and hazard attributes  

    Parameters
        - network_shapefile - Shapefile of edge LineStrings or node Points 
        - polygon_shapefile - Shapefile of boundary Polygons
        - hazard_dictionary - Dictionary of hazard attributes
        - data_dictionary - Dictionary of network-hazard-boundary intersection attributes
        - network_type - String value -'edges' or 'nodes' - Default = 'nodes'
        - name_province - String name of province if needed - Default = ''    

    Outputs
        data_dictionary - Dictionary of network-hazard-boundary intersection attributes:
            - edge_id/node_id - String name of intersecting edge ID or node ID
            - length - Float length of intersection of edge LineString and hazard Polygon: Only for edges 
            - province_id - String/Integer ID of Province
            - province_name - String name of Province in English
            - district_id - String/Integer ID of District
            - district_name - String name of District in English
            - commune_id - String/Integer ID of Commune
            - commune_name - String name of Commune in English
            - hazard_attributes - Dictionary of all attributes from hazard dictionary   
    """
    line_gpd = gpd.read_file(network_shapefile)
    poly_gpd = gpd.read_file(polygon_shapefile)


    if len(line_gpd.index) > 0 and len(poly_gpd.index) > 0:
        line_gpd.columns = map(str.lower, line_gpd.columns)
        poly_gpd.columns = map(str.lower, poly_gpd.columns)
        if name_province != '':
            poly_gpd = poly_gpd[poly_gpd['pro_name_e'] == name_province]

        # create spatial index
        poly_sindex = poly_gpd.sindex

        poly_sindex = poly_gpd.sindex
        for l_index, lines in line_gpd.iterrows():
            intersected_polys = poly_gpd.iloc[list(
                poly_sindex.intersection(lines.geometry.bounds))]
            for p_index, poly in intersected_polys.iterrows():
                if (lines['geometry'].intersects(poly['geometry']) is True) and (poly.geometry.is_valid is True):
                    if network_type == 'edges':
                        value_dictionary = {'edge_id': lines['edge_id'], 'length': lines['length'],
                                            'province_id': poly['province_i'], 'province_name': poly['pro_name_e'],
                                            'district_id': poly['district_i'], 'district_name': poly['dis_name_e'],
                                            'commune_id': poly['commune_id'], 'commune_name': poly['name_eng']}
                    elif network_type == 'nodes':
                        value_dictionary = {'node_id': lines['node_id'],
                                            'province_id': poly['province_i'], 'province_name': poly['pro_name_e'],
                                            'district_id': poly['district_i'], 'district_name': poly['dis_name_e'],
                                            'commune_id': poly['commune_id'], 'commune_name': poly['name_eng']}

                    data_dictionary.append({**value_dictionary, **hazard_dictionary})

    del line_gpd, poly_gpd
    return data_dictionary

def create_hazard_attributes_for_network(intersection_dir,sector,hazard_files,hazard_df,bands,thresholds,commune_shape,network_type='',name_province=''):
    """
    Extract results of network edges/nodes and hazard intersections to collect network-hazard intersection attributes 
        - Combine with boundary Polygons to collect network-hazard-boundary intersection attributes 
        - Write final results to an Excel sheet 

    Parameters
        - intersection_dir - String Path to Directory where the network-hazard shapefile results are stored
        - sector - String name of transport mode
        - hazard_files - List of string names of all hazard files
        - hazard_df - Pandas DataFrame of hazard attributes   
        - bands - List of intergers values of hazard bands
        - thresholds - List of intergers values of hazard thresholds
        - commune_shahe - Shapefile of commune boundaries and attributes
        - network_type - String value -'edges' or 'nodes': Default = 'nodes'
        - name_province - String name of province if needed: Default = ''    

    Outputs
        data_df - Pandas DataFrame of network-hazard-boundary intersection attributes:
            - edge_id/node_id - String name of intersecting edge ID or node ID
            - length - Float length of intersection of edge LineString and hazard Polygon: Only for edges 
            - province_id - String/Integer ID of Province
            - province_name - String name of Province in English
            - district_id - String/Integer ID of District
            - district_name - String name of District in English
            - commune_id - String/Integer ID of Commune
            - commune_name - String name of Commune in English
            - sector - String name of transport mode
            - hazard_type - String name of hazard type
            - model - String name of hazard model
            - year - String name of hazard year
            - climate_scenario - String name of hazard scenario
            - probability - Float/String value of hazard probability
            - band_num - Integer value of hazard band
            - min_val - Integer value of minimum value of hazard threshold
            - max_val - Integer value of maximum value of hazard threshold
            - length - Float length of intersection of edge LineString and hazard Polygon: Only for edges    
    """
    data_dict = []
    for root, dirs, files in os.walk(intersection_dir):
        for file in files:
            if file.endswith(".shp"):
                hazard_dict = {}
                hazard_dict['sector'] = sector
                hazard_shp = os.path.join(root, file)
                hz_file = [h for h in hazard_files if h in file][0]
                hazard_dict['hazard_type'] = hazard_df.loc[hazard_df.file_name ==
                                                            hz_file].hazard_type.values[0]
                hazard_dict['model'] = hazard_df.loc[hazard_df.file_name ==
                                                        hz_file].model.values[0]
                hazard_dict['year'] = hazard_df.loc[hazard_df.file_name ==
                                                    hz_file].year.values[0]
                hazard_dict['climate_scenario'] = hazard_df.loc[hazard_df.file_name ==
                                                                hz_file].climate_scenario.values[0]
                hazard_dict['probability'] = hazard_df.loc[hazard_df.file_name ==
                                                            hz_file].probability.values[0]
                band_type = hazard_df.loc[hazard_df.file_name == hz_file].banded.values[0]
                if str(band_type) == 'True':
                    hazard_dict['band_num'] = [
                        b for b in bands if '{}_band'.format(b) in file][0]
                    band_names = hazard_df.loc[hazard_df.file_name ==
                                                hz_file].bands.values[0].split(',')
                    hazard_dict['band_name'] = [
                        b for b in band_names if str(hazard_dict['band_num']) in b][0]
                    hazard_dict['min_val'] = 0
                    hazard_dict['max_val'] = 0
                else:
                    hazard_dict['band_num'] = 0
                    hazard_dict['band_name'] = 'none'
                    hazard_thrs = [(thresholds[t], thresholds[t+1]) for t in range(len(thresholds)-1)
                                           if '{0}m-{1}m'.format(thresholds[t], thresholds[t+1]) in file][0]
                    hazard_dict['min_val'] = hazard_thrs[0]
                    hazard_dict['max_val'] = hazard_thrs[1]

                data_dict = spatial_scenario_selection(
                            hazard_shp, commune_shape, hazard_dict, data_dict,
                            network_type = network_type,name_province = name_province)

    data_df = pd.DataFrame(data_dict)
    data_df_cols = data_df.columns.values.tolist()
    if 'length' in data_df_cols:
        selected_cols = [cols for cols in data_df_cols if cols != 'length']
        data_df = data_df.groupby(selected_cols)['length'].sum().reset_index()

    return data_df

def main():
    """
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

    """Supply input data and parameters
    """
    provinces = ['Lao Cai','Binh Dinh','Thanh Hoa']
    modes = ['road','rail','air','inland','coastal']
    out_modes = ['national_roads', 'national_rail', 'air_ports', 'inland_ports', 'sea_ports']
    bands = [3, 4, 5]
    thresholds = [1, 2, 3, 4, 999]
    province_results = 'Yes'
    national_results = 'Yes'

    """Give the paths to the input data files
    """
    commune_shp = os.path.join(data_path, 'Vietnam_boundaries',
                                'who_boundaries', 'who_communes.shp')

    hazard_description_file = os.path.join(
        data_path, 'Hazard_data', 'hazard_data_folder_data_info.xlsx')
    hazard_sheet = 'file_contents'
    
    """Specify the output files and paths to be created
    """
    output_dir = os.path.join(output_path, 'hazard_scenarios')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    """Read hazard datasets desciptions
    """
    print ('* Reading hazard datasets desciptions')
    hazard_df = pd.read_excel(hazard_description_file, sheet_name=hazard_sheet)
    hazard_files = hazard_df['file_name'].values.tolist()

    """Process province scale results
    """
    if province_results == 'Yes':
        print ('* Processing province scale results') 
        data_excel = os.path.join(
            output_dir,'province_roads_hazard_intersections.xlsx')
        prov_excel_writer = pd.ExcelWriter(data_excel)
        for province in provinces:
            province_name = province.replace(' ', '').lower()
            intersection_dir = os.path.join(
                output_path, 
                'networks_hazards_intersection_shapefiles',
                '{}_roads_hazards_intersections'.format(province_name))

            data_df = create_hazard_attributes_for_network(
                intersection_dir,'Roads',hazard_files,hazard_df,
                bands,thresholds,commune_shp,network_type='edges',
                name_province=province)
            data_df.to_excel(prov_excel_writer, province_name, index=False)
            prov_excel_writer.save()
            del data_df

    """Process national scale results
    """
    if national_results == 'Yes':
        print ('* Processing national scale results') 
        data_excel = os.path.join(
            output_dir,'national_scale_hazard_intersections.xlsx')
        nat_excel_writer = pd.ExcelWriter(data_excel)
        for m in range(len(modes)):
            intersection_dir = os.path.join(
                output_path, 
                'networks_hazards_intersection_shapefiles',
                '{}_hazard_intersections'.format(out_modes[m]))

            if modes[m] in ['road','rail']:
                ntype = 'edges'
            else:
                ntype = 'nodes'
            data_df = create_hazard_attributes_for_network(
                intersection_dir,modes[m],hazard_files,hazard_df,
                bands,thresholds,commune_shp,network_type=ntype,
                name_province='')

            data_df.to_excel(nat_excel_writer, modes[m], index=False)
            nat_excel_writer.save()
            del data_df




if __name__ == "__main__":
    main()
