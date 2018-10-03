
"""
DOES NOT WORK WITH GEOPANDAS v0.3.0 WHICH GIVES ERROR IF SHAPEFILE IS EMPTY!!

Intersect hazards and network line and point geometries 
With hazatd polygons

Write final results to Shapefiles

Input data requirements
-----------------------
1. Correct paths to all files and correct input parameters 

2. Shapefiles of network edges or nodes
    Should contain following column names and attributes:
        edge_id or node_id - String/Integer/Float Edge ID or Node ID of network
        geometry - Shapely geometry of edges as LineStrings or nodes as Points

3. Shapefile of hazards
    Should contain following column names and attributes:
        geometry - Shapely geometry of hazard Polygon

Results
-------
Edge shapefiles
    Contains:
        edge_id - String name of intersecting edge ID
        length - Float length of intersection of edge LineString and hazard Polygon
        geometry - LineString geometry of intersection of edge LineString and hazard Polygon

Node Shapefile
    Contains:
        node_id - String name of intersecting node ID
        geometry - Point geometry of intersecting node ID

"""
import itertools
import os
import sys

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from vtra.utils import line_length, load_config


def networkedge_hazard_intersection(edge_shapefile, hazard_shapefile, output_shapefile):
    """
    Intersect network edges and hazards and write results to shapefiles

    Parameters
    ---------
    edge_shapefile - Shapefile of network LineStrings 
    hazard_shapefile - Shapefile of hazard Polygons
    output_shapefile - String name of edge-hazard shapefile for storing results   

    Outputs
    -------
    output_shapefile - Shapefile
        Contains:
            edge_id - String name of intersecting edge ID
            length - Float length of intersection of edge LineString and hazard Polygon
            geometry - LineString geometry of intersection of edge LineString and hazard Polygon
    """
    print ('* Starting {} and {} intersections'.format(edge_shapefile,hazard_shapefile))
    line_gpd = gpd.read_file(edge_shapefile)
    poly_gpd = gpd.read_file(hazard_shapefile)

    if len(line_gpd.index) > 0 and len(poly_gpd.index) > 0:
        line_gpd.columns = map(str.lower, line_gpd.columns)
        poly_gpd.columns = map(str.lower, poly_gpd.columns)

        line_bounding_box = line_gpd.total_bounds
        line_bounding_box_coord = list(itertools.product([line_bounding_box[0], line_bounding_box[2]], [
                                       line_bounding_box[1], line_bounding_box[3]]))
        line_bounding_box_geom = Polygon(line_bounding_box_coord)
        line_bounding_box_gpd = gpd.GeoDataFrame(pd.DataFrame(
            [[1], [line_bounding_box_geom]]).T, crs='epsg:4326')
        line_bounding_box_gpd.columns = ['ID', 'geometry']

        poly_bounding_box = poly_gpd.total_bounds
        poly_bounding_box_coord = list(itertools.product([poly_bounding_box[0], poly_bounding_box[2]], [
                                       poly_bounding_box[1], poly_bounding_box[3]]))
        poly_bounding_box_geom = Polygon(poly_bounding_box_coord)
        poly_bounding_box_gpd = gpd.GeoDataFrame(pd.DataFrame(
            [[1], [poly_bounding_box_geom]]).T, crs='epsg:4326')
        poly_bounding_box_gpd.columns = ['ID', 'geometry']

        poly_sindex = poly_bounding_box_gpd.sindex

        selected_polys = poly_bounding_box_gpd.iloc[list(
            poly_sindex.intersection(line_bounding_box_gpd['geometry'].iloc[0].bounds))]
        if len(selected_polys.index) > 0:
            data = []
            poly_sindex = poly_gpd.sindex
            for l_index, lines in line_gpd.iterrows():
                intersected_polys = poly_gpd.iloc[list(
                    poly_sindex.intersection(lines.geometry.bounds))]
                for p_index, poly in intersected_polys.iterrows():
                    if (lines['geometry'].intersects(poly['geometry']) is True) and (poly.geometry.is_valid is True):
                        data.append({'edge_id': lines['edge_id'], 'length': 1000.0*line_length(lines['geometry'].intersection(
                            poly['geometry'])), 'geometry': lines['geometry'].intersection(poly['geometry'])})
            if data:
                intersections_data = gpd.GeoDataFrame(
                    data, columns=['edge_id', 'length', 'geometry'], crs='epsg:4326')
                intersections_data.to_file(output_shapefile)

                del intersections_data

    del line_gpd, poly_gpd


def networknode_hazard_intersection(node_shapefile, hazard_shapefile, output_shapefile):
    """
    Intersect network nodes and hazards and write results to shapefiles

    Parameters
    ---------
    node_shapefile - Shapefile of network Points 
    hazard_shapefile - Shapefile of hazard Polygons
    output_shapefile - String name of node-hazard shapefile for storing results   

    Outputs
    -------
    output_shapefile - Shapefile
        Contains:
            node_id - String name of intersecting node ID
            geometry - Point geometry of intersecting node ID
    """
    print ('* Starting {} and {} intersections'.format(node_shapefile,hazard_shapefile))
    point_gpd = gpd.read_file(node_shapefile)
    poly_gpd = gpd.read_file(hazard_shapefile)

    if len(point_gpd.index) > 0 and len(poly_gpd.index) > 0:
        point_gpd.columns = map(str.lower, point_gpd.columns)
        poly_gpd.columns = map(str.lower, poly_gpd.columns)
        data = []
        # create spatial index
        poly_sindex = poly_gpd.sindex
        for pt_index, points in point_gpd.iterrows():
            intersected_polys = poly_gpd.iloc[list(
                poly_sindex.intersection(points.geometry.bounds))]
            if len(intersected_polys.index) > 0:
                data.append({'node_id': points['node_id'], 'geometry': points['geometry']})
        if data:
            intersections_data = gpd.GeoDataFrame(
                data, columns=['node_id', 'geometry'], crs='epsg:4326')
            intersections_data.to_file(output_shapefile)

            del intersections_data

    del point_gpd, poly_gpd

def intersect_networks_and_all_hazards(hazard_dir,network_file_path,network_file_name,output_file_path,network_type = ''):
    """
    Walk through all hazard files and select network-hazard intersection criteria
    Call other functions accordingly 

    Parameters
    ---------
    hazard_dir - String name of directory where all hazard shapefiles are stored 
    network_file_path - String name of directory where network shapefile is stored 
    network_file_name - String name network shapefile
    output_file_path - String name of directory where network-hazard instersection result shapefiles will be stored
    network_type - String values
        'edges' or 'nodes'    
    
    Outputs
    -------
    Edge or Node shapefiles
    """
    for root, dirs, files in os.walk(hazard_dir):
        for file in files:
            if file.endswith(".shp"):
                hazard_file = os.path.join(root, file)
                out_shp_name = network_file_name[:-4] + '_' + file
                output_file = os.path.join(output_file_path,out_shp_name)
                if network_type == 'edges':
                    networkedge_hazard_intersection(network_file_path, hazard_file, output_file)
                elif network_type == 'nodes':
                    networknode_hazard_intersection(network_file_path, hazard_file, output_file)


def main():
    """
    Specify the paths from where you to read and write:
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    1. Names of the three Provinces
        List of string types 
    2. Paths of the mode files
        List of tuples of strings
    3. Names of modes
        List of strings
    4. Names of output modes
        List of strings
    5. Condition 'Yes' or 'No' is the users wants to process
        Province scale results
        National scale results 

    Give the paths to the input data files:
    1. Hazard directory
    
    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    provinces = ['Lao Cai','Binh Dinh','Thanh Hoa']

    modes_file_paths = [('Roads','national_roads'), ('Railways','national_rail'), ('Airports','airnetwork'), ('Waterways','waterways'), ('Waterways','waterways')]
    modes = ['road', 'rail', 'air', 'inland', 'coastal']
    out_modes = ['national_roads', 'national_rail', 'air_ports', 'inland_ports', 'sea_ports']
    province_results = 'Yes'
    national_results = 'Yes'

    """Give the paths to the input data files
    """
    hazard_dir = os.path.join(data_path, 'Hazard_data')

    """Specify the output files and paths to be created 
    """
    output_dir = os.path.join(output_path, 'networks_hazards_intersection_shapefiles')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)
    """Start province roads and all hazards intersections
    """
    if province_results == 'Yes':
        for province in provinces:
            province_name = province.replace(' ','').lower()
            road_shp_path = os.path.join(data_path,'Roads','{}_roads'.format(province_name),'vietbando_{}_edges.shp'.format(province_name))
            road_shp = 'vietbando_{}_edges.shp'.format(province_name)

            output_dir = os.path.join(output_path, 'networks_hazards_intersection_shapefiles','{}_roads_hazards_intersections'.format(province_name))
            if os.path.exists(output_dir) == False:
                os.mkdir(output_dir)

            print ('* Starting {} roads and all hazards intersections'.format(province))
            intersect_networks_and_all_hazards(hazard_dir,road_shp_path,road_shp,output_dir,network_type = 'edges')

    if national_results == 'Yes':
        for m in range(len(modes)):
            mode_data_path = os.path.join(
                data_path, modes_file_paths[m][0], modes_file_paths[m][1])
            if modes[m] in ['road', 'rail']:
                for mode_file in os.listdir(mode_data_path):
                    try:
                        if mode_file.endswith(".shp") and 'edges' in mode_file.lower().strip():
                            edges_in = os.path.join(mode_data_path, mode_file)
                            edges_name = mode_file
                    except:
                        return ('Network edge file necessary')

                output_dir = os.path.join(output_path, 'networks_hazards_intersection_shapefiles','{}_hazard_intersections'.format(out_modes[m]))
                if os.path.exists(output_dir) == False:
                    os.mkdir(output_dir)

                print ('* Starting national {} and all hazards intersections'.format(modes[m]))
                intersect_networks_and_all_hazards(hazard_dir,edges_in,edges_name,output_dir,network_type = 'edges')

            elif modes[m] in ['air', 'inland', 'coastal']:
                for mode_file in os.listdir(mode_data_path):
                    try:
                        if mode_file.endswith(".shp") and 'nodes' in mode_file.lower().strip():
                            nodes_in = os.path.join(mode_data_path, mode_file)
                            nodes_name = mode_file
                    except:
                        return ('Network node file necessary')
                
                output_dir = os.path.join(output_path, 'networks_hazards_intersection_shapefiles','{}_hazard_intersections'.format(out_modes[m]))
                if os.path.exists(output_dir) == False:
                    os.mkdir(output_dir)

                print ('* Starting national {} and all hazards intersections'.format(modes[m]))
                intersect_networks_and_all_hazards(hazard_dir,nodes_in,nodes_name,output_dir,network_type = 'nodes')


if __name__ == "__main__":
    main()
