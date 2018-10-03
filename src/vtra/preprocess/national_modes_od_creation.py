"""Create national scale OD matrices at node and province levels from: 
1. VITRANSS2 province-scale OD data
2. IFPRI crop data at 1km resolution

Input data requirements
-----------------------
1. Correct paths to all files and correct input parameters 

Results
------- 
"""
import os
import subprocess
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from vtra.transport_network_creation import *
from vtra.utils import *


def vitranss2_od_split(vitranss2_od_data_file,modes,o_id_col,d_id_col):
    """Create VITRANSS2 OD modal split estimates
    
    By combining the commodity-wise OD values with the mode-wise OD values
    And estimating the OD modal-split from the mode-wsie OD values  

    Parameters
    ----------
    vitranss2_od_data_file - Excel file with VITRANSS 2 OD data
    modes - List of strings of mode types, e.g. ['road', 'rail', 'air', 'inland', 'coastal']
    o_id_col - String name of name of Origin column 
    d_id_col - String name of name of Destination column  

    Outputs
    -------
    od_fracs - Pandas dataframe of VITRANSS 2 commodity OD data with modal splits 
    """
    od_data_modes = pd.read_excel(vitranss2_od_data_file, sheet_name='mode').fillna(0)
    od_data_com = pd.read_excel(vitranss2_od_data_file, sheet_name='goods').fillna(0)

    od_data_modes.columns = map(str.lower, od_data_modes.columns)
    od_data_modes['total'] = od_data_modes[modes].sum(axis=1)
    for m in modes:
        od_data_modes[m] = od_data_modes[m]/od_data_modes['total'].replace(np.inf, 0)

    od_data_modes = od_data_modes.fillna(0)
    od_fracs = od_data_modes[[o_id_col,d_id_col] + modes]
    od_fracs = pd.merge(od_fracs, od_data_com, how='left', on=[o_id_col,d_id_col]).fillna(0)

    del od_data_modes, od_data_com
    return od_fracs

def ifpri_crop_od_split(vitranss2_od_data_file,modes,o_id_col,d_id_col,crop_cols):
    """Create IFPRI crop OD modal split estimates from VITRANSS2 OD values
    
    By combining the crop-wise OD values with the mode-wise OD values
    And estimating the OD modal-split from the mode-wsie OD values  

    Parameters
    ----------
    vitranss2_od_data_file - Excel file with VITRANSS2 OD data
    modes - List of strings of mode types, e.g. ['road', 'rail', 'air', 'inland', 'coastal']
    o_id_col - String name of name of Origin column 
    d_id_col - String name of name of Destination column
    crop_cols - List of strings of crop names in VITRANSS 2 OD data  

    Outputs
    -------
    od_fracs_crops - Pandas dataframe of VITRANSS2 crop OD data with modal splits 
    """
    od_data_modes = pd.read_excel(vitranss2_od_data_file, sheet_name='mode').fillna(0)
    od_data_com = pd.read_excel(vitranss2_od_data_file, sheet_name='goods').fillna(0)

    od_data_modes.columns = map(str.lower, od_data_modes.columns)
    od_data_modes['total'] = od_data_modes[modes].sum(axis=1)
    for m in modes:
        od_data_modes[m] = od_data_modes[m]/od_data_modes['total'].replace(np.inf, 0)

    od_data_modes = od_data_modes.fillna(0)
    od_fracs_crops = od_data_modes[[o_id_col,d_id_col] + modes]
    for cr in crop_cols:
        od_data_com_sums = od_data_com.groupby([o_id_col,d_id_col]).agg({cr: 'sum'})
        od_com_frac = od_data_com_sums.groupby(level=0).apply(lambda x: x/float(x.sum()))
        od_com_frac = od_com_frac.reset_index(level=[o_id_col,d_id_col])
        od_fracs_crops = pd.merge(od_fracs_crops, od_com_frac,
                                  how='left', on=[o_id_col,d_id_col]).fillna(0)
    
    del od_data_com, od_data_com_sums, od_com_frac
    return od_fracs_crops

def riceatlas_crop_minmax(riceatlas_crop_file,crop_month_fields):
    """Create MIN_MAX fractions of rice production estimates for provinces
    
    By reading data from the RiceAtlas data
    And estimating the minimum and maximum monthly values > 0 

    Parameters
    ----------
    riceatlas_crop_file - Shapefile with RiceAtlas data
    crop_month_fields - List of strings of names of columns indicating rice monthly production 

    Outputs
    -------
    rice_prod_months - Geopandas dataframe of RiceAtlas crop values with minimum and maximum monthly production as fraction of annual production
    """
    rice_prod_months = gpd.read_file(riceatlas_crop_file)
    rice_prod_months['total_prod'] = rice_prod_months[crop_month_fields].sum(axis=1)
    rice_prod_months['min_tons'] = rice_prod_months[rice_prod_months[crop_month_fields] > 0].min(
        axis=1)
    rice_prod_months['max_tons'] = rice_prod_months[rice_prod_months[crop_month_fields] > 0].max(
        axis=1)

    rice_prod_months['min_frac'] = rice_prod_months['min_tons']/rice_prod_months['total_prod']
    rice_prod_months['max_frac'] = rice_prod_months['max_tons']/rice_prod_months['total_prod']
    return rice_prod_months

def assign_province_name_id_to_nodes(province_path,nodes_in,province_name_col,province_id_col):
    """Match the nodes to their province names and province IDs
    
    By finding the province that contains or is nearest to the node 

    Parameters
    ----------
    province_path - Path of province shapefile 
    nodes_in - Path of nodes shapefile
    province_name_col - String name of column containing province names
    province_id_col - String name of column containing province ID's that match VITRANSS 2 OD ids

    Outputs
    -------
    nodes - Geopandas dataframe of nodes with new columns called province_name and od_id
    """

    # load provinces and get geometry of the right province
    provinces = gpd.read_file(province_path)
    provinces = provinces.to_crs({'init': 'epsg:4326'})
    sindex_provinces = provinces.sindex
    # load nodes of the network
    nodes = gpd.read_file(nodes_in)
    nodes = nodes.to_crs({'init': 'epsg:4326'})
    nodes.columns = map(str.lower, nodes.columns)

    nodes['province_name'] = nodes.apply(lambda x: extract_gdf_values_containing_nodes(
        x, sindex_provinces, provinces,province_name_col), axis=1)
    nodes['od_id'] = nodes.apply(lambda x: extract_gdf_values_containing_nodes(
        x, sindex_provinces, provinces, province_id_col), axis=1)

    del provinces

    return nodes

def assign_road_weights(nodes,edges_in,aadt_column):
    """Assign weights to nodes on the road network
    
    By finding the total AADT counts converging on the node 

    Parameters
    ----------
    nodes - Geopandas dataframe of nodes  
    edges_in - Path of edges shapefile
    aadt_column - String name of column containing AADT values

    Outputs
    -------
    nodes - Geopandas dataframe of nodes with new column called weight
    """

    edges_df = gpd.read_file(edges_in)
    edges_df.columns = map(str.lower, edges_df.columns)
    nodes_vehs = list(zip(edges_df['from_node'].values.tolist(
    ), edges_df['to_node'].values.tolist(), edges_df[aadt_column].values.tolist()))
    nd_veh_list = []
    for nd in nodes['node_id'].values.tolist():
        veh = 0.5*sum([int(v) for (f, t, v) in nodes_vehs if nd == f or nd == t])
        nd_veh_list.append((nd, veh))

    gdf_pops = pd.DataFrame(nd_veh_list, columns=['node_id', 'weight'])
    del nd_veh_list
    nodes = pd.merge(nodes, gdf_pops, how='left', on=['node_id']).fillna(0)
    del gdf_pops, edges_df

    return nodes

def assign_node_weights_by_commune_population_proximity(commune_path,nodes,commune_pop_col):
    """Assign weights to nodes based on their nearest commune populations 
    
    By finding the communes that intersect with the Voronoi extents of nodes 

    Parameters
    ----------
    commune_path - Path of commune shapefile 
    nodes_in - Path of nodes shapefile
    commune_pop_col - String name of column containing commune population values

    Outputs
    -------
    nodes - Geopandas dataframe of nodes with new column called weight
    """

    # load provinces and get geometry of the right communes within the provinces
    communes = gpd.read_file(commune_path)
    communes = communes.to_crs({'init': 'epsg:4326'})
    sindex_communes = communes.sindex

    # create Voronoi polygons for the nodes
    xy_list = []
    for iter_, values in nodes.iterrows():
        xy = list(values.geometry.coords)
        xy_list += [list(xy[0])]

    vor = Voronoi(np.array(xy_list))
    regions, vertices = voronoi_finite_polygons_2d(vor)
    min_x = vor.min_bound[0] - 0.1
    max_x = vor.max_bound[0] + 0.1
    min_y = vor.min_bound[1] - 0.1
    max_y = vor.max_bound[1] + 0.1

    mins = np.tile((min_x, min_y), (vertices.shape[0], 1))
    bounded_vertices = np.max((vertices, mins), axis=0)
    maxs = np.tile((max_x, max_y), (vertices.shape[0], 1))
    bounded_vertices = np.min((bounded_vertices, maxs), axis=0)

    box = Polygon([[min_x, min_y], [min_x, max_y], [max_x, max_y], [max_x, min_y]])

    poly_list = []
    for region in regions:
        polygon = vertices[region]
        # Clipping polygon
        poly = Polygon(polygon)
        poly = poly.intersection(box)
        poly_list.append(poly)

    poly_index = list(np.arange(0, len(poly_list), 1))
    poly_df = pd.DataFrame(list(zip(poly_index, poly_list)),
                                   columns=['gid', 'geometry'])
    gdf_voronoi = gpd.GeoDataFrame(poly_df, crs='epsg:4326')
    gdf_voronoi['node_id'] = gdf_voronoi.apply(
        lambda x: extract_nodes_within_gdf(x, nodes, 'node_id'), axis=1)
    
    gdf_voronoi[commune_pop_col] = 0
    gdf_voronoi = assign_value_in_area_proportions(communes, gdf_voronoi, commune_pop_col)

    gdf_voronoi.rename(columns={'population': 'weight'}, inplace=True)
    gdf_pops = gdf_voronoi[['node_id', 'weight']]
    del gdf_voronoi, poly_list, poly_df

    nodes = pd.merge(nodes, gdf_pops, how='left', on=['node_id']).fillna(0)
    del gdf_pops, communes

    return nodes

def assign_industry_od_flows_to_nodes(national_ods_df,ind_cols,modes_df,modes,od_fracs,o_id_col,d_id_col):
    """Assign VITRANSS 2 OD flows to nodes 

    Parameters
    ----------
    national_ods_df - List of lists of Pandas dataframes 
    ind_cols - List of strings of names of indsutry columns
    modes_df - List of Geopnadas dataframes with nodes of each transport mode
    modes - List of strings of names of transport modes
    od_fracs - Pandas dataframe of Industry OD flows and modal splits
    o_id_col - String name of Origin province ID column 
    d_id_col - String name of Destination province ID column 

    Outputs
    -------
    national_ods_df - List of Lists of Pandas dataframes 
        Each dataframe has columns:
            origin - Origin node ID
            o_region - Origin province name
            destination - Destination node ID
            d_region - Destination province ID
            ind - Tonnage values for the named industry 

    """
    for ind in ind_cols:
        national_ods_modes_df = []
        for m in range(len(modes_df)):
            nodes = modes_df[m]
            od_nodes_regions = list(zip(nodes['node_id'].values.tolist(), nodes['province_name'].values.tolist(
            ), nodes['od_id'].values.tolist(), nodes['weight'].values.tolist()))
            ind_mode = modes[m] + '_' + ind
            od_fracs[ind_mode] = od_fracs[modes[m]]*od_fracs[ind]

            od_fracs_ind = od_fracs[[o_id_col, d_id_col, ind_mode]]
            od_fracs_ind = od_fracs_ind[od_fracs_ind[ind_mode] > 0]
            od_flows = list(zip(od_fracs_ind[o_id_col].values.tolist(
            ), od_fracs_ind[d_id_col].values.tolist(), od_fracs_ind[ind_mode].values.tolist()))
            origins = list(set(od_fracs_ind[o_id_col].values.tolist()))
            destinations = list(set(od_fracs_ind[d_id_col].values.tolist()))

            # print (od_flows)
            od_list = []
            for o in origins:
                for d in destinations:
                    fval = [fl for (org, des, fl) in od_flows if org == o and des == d]
                    if len(fval) == 1 and fval[0] > 0:
                        o_matches = [(item[0], item[1], item[3])
                                     for item in od_nodes_regions if item[2] == o]
                        if len(o_matches) > 0:
                            for o_vals in o_matches:
                                o_val = 1.0*fval[0]*(1.0*o_vals[2])
                                o_node = o_vals[0]
                                o_region = o_vals[1]
                                d_matches = [(item[0], item[1], item[3])
                                             for item in od_nodes_regions if item[2] == d]
                                if len(d_matches) > 0:
                                    for d_vals in d_matches:
                                        od_val = 1.0*o_val*(1.0*d_vals[2])
                                        d_node = d_vals[0]
                                        d_region = d_vals[1]
                                        if od_val > 0 and o_node != d_node:
                                            od_list.append(
                                                (o_node, o_region, d_node, d_region, od_val))


            national_ods_modes_df.append(pd.DataFrame(
                od_list, columns=['origin', 'o_region', 'destination', 'd_region', ind]))
            del od_list, nodes

        national_ods_df.append(national_ods_modes_df)

    return national_ods_df

def assign_daily_min_max_tons_rice(crop_df, rice_prod_df):
    """Estimate minimum and maximum daily rice tonnages

    Parameters
    ----------
    crop_df - Geopandas dataframe of crop points with annual tonnages 
    rice_prod_df - Geopandas dataframe of RiceAtlas crop values with minimum and maximum monthly production as fraction of annual production

    Outputs
    -------
    crop_df - Geopandas dataframe of crop points 
        With new columns:
            min_rice - Minimum daily rice tonnages
            max_rice - Maximum daily rice tonnages
       
    """
    sindex_rice_prod_df = rice_prod_df.sindex

    crop_df['min_frac'] = crop_df.geometry.apply(
        lambda x: get_nearest_node(x, sindex_rice_prod_df, rice_prod_df, 'min_frac'))
    crop_df['max_frac'] = crop_df.geometry.apply(
        lambda x: get_nearest_node(x, sindex_rice_prod_df, rice_prod_df, 'max_frac'))

    crop_df['min_rice'] = 1.0*crop_df['min_frac']*crop_df['tons']/30.0
    crop_df['max_rice'] = 1.0*crop_df['max_frac']*crop_df['tons']/30.0

    return crop_df

def assign_crop_od_flows_to_nodes(national_ods_df,province_path,province_name_col,calc_path,
    crop_data_path,crop_names,rice_prod_months,modes_df,modes,od_fracs_crops,o_id_col,d_id_col):
    """Assign IFPRI crop values to OD nodes based on VITRANSS 2 OD distributions
    
    Based on VITRANSS 2 OD distributions

    Parameters
    ----------
    national_ods_df - List of lists of Pandas dataframes
    province_path - Path of province shapefile
    province_name_col - String name of column containing province names 
    calc_path - Path to store intermediary calculations 
    crop_data_path - Path to crop datasets
    crop_names - List of string of crop names in IFPRI datasets
    rice_prod_months - Geopandas dataframe of RiceAtlas crop values with minimum and maximum monthly production as fraction of annual production
    modes_df - List of Geopnadas dataframes with nodes of each transport mode
    modes - List of strings of names of transport modes
    od_fracs_crops - Pandas dataframe of crop OD distributions and modal splits as per VITRANSS 2 data
    o_id_col - String name of Origin province ID column 
    d_id_col - String name of Destination province ID column 

    Outputs
    -------
    national_ods_df - List of Lists of Pandas dataframes 
        Each dataframe has columns:
            origin - Origin node ID
            o_region - Origin province name
            destination - Destination node ID
            d_region - Destination province ID
            min_crop - Minimum Tonnage values for the named crop
            max_crop - Maximum Tonnage values for the named crop 

    """
    # load provinces and get geometry of the right province
    provinces = gpd.read_file(province_path)
    provinces = provinces.to_crs({'init': 'epsg:4326'})
    sindex_provinces = provinces.sindex

    for file in os.listdir(crop_data_path):
        if file.endswith(".tif") and ('spam_p' in file.lower().strip()):
            fpath = os.path.join(crop_data_path, file)
            crop_name = [cr for cr in crop_names if cr in file.lower().strip()][0]
            outCSVName = os.path.join(calc_path, 'crop_concentrations.csv')
            subprocess.run(["gdal2xyz.py", '-csv', fpath, outCSVName])

            # Load points and convert to geodataframe with coordinates
            load_points = pd.read_csv(outCSVName, header=None, names=[
                                      'x', 'y', 'tons'], index_col=None)
            load_points = load_points[load_points['tons'] > 0]

            geometry = [Point(xy) for xy in zip(load_points.x, load_points.y)]
            load_points = load_points.drop(['x', 'y'], axis=1)
            crop_points = gpd.GeoDataFrame(load_points, crs={'init': 'epsg:4326'}, geometry=geometry)

            del load_points

            if crop_name == 'rice':
                crop_points = assign_daily_min_max_tons_rice(crop_points, rice_prod_months)
            else:
                crop_points['min_{}'.format(crop_name)] = 1.0*crop_points['tons']/365.0
                crop_points['max_{}'.format(crop_name)] = 1.0*crop_points['tons']/365.0


            crop_points['province_name'] = crop_points.apply(lambda x: extract_gdf_values_containing_nodes(
                x, sindex_provinces, provinces, province_name_col), axis=1)
            
            national_ods_modes_df = []
            for m in range(len(modes_df)):
                nodes = modes_df[m]
                crop_pts = crop_points.copy(deep=True)
                crop_pts['node_id'] = crop_pts.apply(lambda x: get_nearest_node_within_region(
                    x, nodes, 'node_id', 'province_name'), axis=1)
                
                crop_pts = crop_pts[crop_pts['node_id'] != '']
                crop_pts = crop_pts[['node_id', 'min_{}'.format(
                    crop_name), 'max_{}'.format(crop_name)]]
                crop_nodes = crop_pts.groupby(['node_id'])['min_{}'.format(
                    crop_name), 'max_{}'.format(crop_name)].sum().reset_index()
                crop_nodes = crop_nodes.reset_index()


                del crop_pts
                nodes = pd.merge(nodes, crop_nodes, how='left', on=['node_id']).fillna(0)
                del crop_nodes

                crop_mode = modes[m] + '_' + crop_name
                if crop_name in ('rice', 'cereal', 'wheat'):
                    od_fracs_crops[crop_mode] = od_fracs_crops[modes[m]]*od_fracs_crops['rice']
                else:
                    od_fracs_crops[crop_mode] = od_fracs_crops[modes[m]] * \
                        od_fracs_crops['indust-cro']

                od_nodes_regions = list(zip(nodes['node_id'].values.tolist(), nodes['province_name'].values.tolist(), nodes['od_id'].values.tolist(
                ), nodes['min_{}'.format(crop_name)].values.tolist(), nodes['max_{}'.format(crop_name)].values.tolist(), nodes['weight'].values.tolist()))

                od_fracs_ind = od_fracs_crops[[o_id_col, d_id_col, crop_mode]]
                od_fracs_ind = od_fracs_ind[od_fracs_ind[crop_mode] > 0]
                od_flows = list(zip(od_fracs_ind[o_id_col].values.tolist(
                ), od_fracs_ind[d_id_col].values.tolist(), od_fracs_ind[crop_mode].values.tolist()))
                origins = list(set(od_fracs_ind[o_id_col].values.tolist()))
                destinations = list(set(od_fracs_ind[d_id_col].values.tolist()))

                od_list = []
                for o in origins:
                    for d in destinations:
                        fval = [fl for (org, des, fl) in od_flows if org == o and des == d]
                        if len(fval) == 1 and fval[0] > 0:
                            o_matches = [(item[0], item[1], item[3], item[4])
                                         for item in od_nodes_regions if item[2] == o]
                            if len(o_matches) > 0:
                                for o_vals in o_matches:
                                    o_val_min = 1.0*fval[0]*o_vals[2]
                                    o_val_max = 1.0*fval[0]*o_vals[3]
                                    o_node = o_vals[0]
                                    o_region = o_vals[1]
                                    d_matches = [(item[0], item[1], item[5])
                                                 for item in od_nodes_regions if item[2] == d]
                                    if len(d_matches) > 0:
                                        for d_vals in d_matches:
                                            od_val_min = 1.0*o_val_min*d_vals[2]
                                            od_val_max = 1.0*o_val_max*d_vals[2]
                                            d_node = d_vals[0]
                                            d_region = d_vals[1]
                                            if od_val_max > 0 and o_node != d_node:
                                                od_list.append(
                                                    (o_node, o_region, d_node, d_region, od_val_min, od_val_max))


                national_ods_modes_df.append(pd.DataFrame(od_list, columns=[
                                             'origin', 'o_region', 'destination', 'd_region', 'min_{}'.format(crop_name), 'max_{}'.format(crop_name)]))
                del od_list, nodes

            del crop_points
            national_ods_df.append(national_ods_modes_df)

    del provinces
    return national_ods_df

def main():
    """
    Specify the paths from where you to read and write:
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    1. Names of modes
        List of strings
    2. OD column names in OD file
        String types 
    3. Names of industry columns in VITRANSS2 data
        List of string types
    4. Names of crop columns in VITRANSS2 data
        List of string types
    5. Names of crops in IFPRI crop data
        List of string types
    6. Names of months in Rice Atlas data
        List of string types
    
    Give the paths to the input data files:
    1. Network nodes files
    2. VITRANSS 2 OD file
    3. IFPRI crop data files
    4. Rice Altas data shapefile
    5. Province boundary and stats data shapefile
    6. Commune boundary and stats data shapefile

    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    # modes_file_paths = [('Roads','national_roads'), ('Railways','national_rail'), ('Airports','airnetwork'), ('Waterways','waterways'), ('Waterways','waterways')]
    modes = ['road', 'rail', 'air', 'inland', 'coastal']
    o_id_col = 'o'
    d_id_col = 'd'
    ind_cols = ['sugar', 'wood', 'steel', 'constructi', 'cement',
                'fertilizer', 'coal', 'petroluem', 'manufactur', 'fishery', 'meat']
    crop_cols = ['rice', 'indust-cro']
    crop_names = ['rice', 'cash', 'cass', 'teas',
                  'maiz', 'rubb', 'swpo', 'acof', 'rcof', 'pepp']
    crop_month_fields = ['P_Jan', 'P_Feb', 'P_Mar', 'P_Apr', 'P_May',
                         'P_Jun', 'P_Jul', 'P_Aug', 'P_Sep', 'P_Oct', 'P_Nov', 'P_Dec']

    """Give the paths to the input data files:
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    od_data_file = os.path.join(data_path, 'OD_data', 'OD_transport_data.xlsx')
    crop_data_path = os.path.join(data_path, 'Agriculture_crops', 'crop_data')
    rice_month_file = os.path.join(data_path, 'rice_atlas_vietnam', 'rice_production.shp')
    province_path = os.path.join(data_path, 'Vietnam_boundaries',
                                 'boundaries_stats', 'province_level_stats.shp')
    commune_path = os.path.join(data_path, 'Vietnam_boundaries',
                                'boundaries_stats', 'commune_level_stats.shp')

    """Specify the output files and paths to be created 
    """
    output_dir = os.path.join(output_path, 'flow_ods')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    flow_output_excel = os.path.join(output_dir,'national_scale_flow_ods.xlsx')
    excl_wrtr = pd.ExcelWriter(flow_output_excel)

    flow_output_excel = os.path.join(output_dir, 'national_scale_od_matrix.xlsx')
    excl_wrtr_reg = pd.ExcelWriter(flow_output_excel)

    """Get values for:
    1. VITRANSS2 OD by commidity and modal split
    2. IFPRI OD by crop type and modal split derived from VITRANSS 2 
    3. RiceAtlas MIN-MAX  (> 0) seasonality estimates  
    """
    print ('* Creating OD dataframes')
    od_fracs = vitranss2_od_split(od_data_file,modes,o_id_col,d_id_col)
    od_fracs_crops = ifpri_crop_od_split(od_data_file,modes,o_id_col,d_id_col,crop_cols)
    rice_prod_months = riceatlas_crop_minmax(rice_month_file,crop_month_fields)
   
    """Assign weights to nodes to distribute OD values
    """
    print ('* Assinging weights to nodes to distribute OD values')
    modes_df = []
    for m in range(len(modes)):
        nodes_in,edges_in = get_node_edge_files(network_data_path,modes[m])
        
        nodes = assign_province_name_id_to_nodes(province_path,nodes_in,'name_eng','od_id')

        if modes[m] == 'road':
            nodes = assign_road_weights(nodes,edges_in,'vehicle_co')

        elif modes[m] in ('inland', 'coastal'):
            nodes['weight'] = nodes['tons']

        else:
            nodes = assign_node_weights_by_commune_population_proximity(commune_path,nodes,'population')

        nodes_sums = nodes.groupby(['od_id', 'node_id']).agg({'weight': 'sum'})
        nodes_frac = nodes_sums.groupby(level=0).apply(lambda x: x/float(x.sum()))
        nodes_frac = nodes_frac.reset_index(level=['od_id', 'node_id'])

        nodes.drop('weight', axis=1, inplace=True)
        nodes = pd.merge(nodes, nodes_frac[['node_id', 'weight']],
                         how='left', on=['node_id']).fillna(0)

        modes_df.append(nodes)

        del nodes_frac, nodes_sums, nodes
    

    national_ods_df = []
    """Assign the industry OD flows to nodes 
    """ 
    print ('* Assinging Industry OD flows to nodes')
    national_ods_df = assign_industry_od_flows_to_nodes(national_ods_df,ind_cols,modes_df,modes,od_fracs,o_id_col,d_id_col)

    """Assign the crop OD flows to nodes 
    """
    print ('* Assinging Crop OD flows to nodes')
    national_ods_df = assign_crop_od_flows_to_nodes(national_ods_df,province_path,'name_eng',
        calc_path,crop_data_path,crop_names,rice_prod_months,modes_df,modes,od_fracs_crops,o_id_col,d_id_col)
    
    """Store OD values in Excel Sheets to create:
    1. Node level OD flows
    2. Province level OD matrices  
    """
    print ('* Storing OD flows Excel sheets')
    national_ods_df = list(map(list, zip(*national_ods_df)))
    region_total = []
    for m in range(len(modes_file_paths)):
        all_ods = pd.concat(national_ods_df[m], axis=0,
                            sort='False', ignore_index=True).fillna(0)

        all_min_cols = ind_cols + ['min_{}'.format(c) for c in crop_names]
        all_ods['min_tons'] = all_ods[all_min_cols].sum(axis=1)
        all_max_cols = ind_cols + ['max_{}'.format(c) for c in crop_names]
        all_ods['max_tons'] = all_ods[all_max_cols].sum(axis=1)
        crops_norice = [cr for cr in crop_names if cr != 'rice']
        for cr in crops_norice:
            all_ods.drop('min_{}'.format(cr), axis=1, inplace=True)
            all_ods.rename(columns={'max_{}'.format(cr): cr}, inplace=True)

        all_ods_val_cols = [c for c in all_ods.columns.values.tolist(
        ) if c not in ('origin', 'o_region', 'destination', 'd_region')]
        all_ods = all_ods.groupby(['origin', 'o_region', 'destination', 'd_region'])[
            all_ods_val_cols].sum().reset_index()

        all_ods_regions = all_ods[['o_region', 'd_region'] + all_ods_val_cols]
        all_ods_regions = all_ods_regions.groupby(['o_region', 'd_region'])[
            all_ods_val_cols].sum().reset_index()
        all_ods_regions.to_excel(excl_wrtr_reg, modes[m], index=False)
        excl_wrtr_reg.save()

        region_total.append(all_ods_regions)
        del all_ods_regions

        all_ods = all_ods[all_ods['max_tons'] > 0.5]
        all_ods.to_excel(excl_wrtr, modes[m], index=False)
        excl_wrtr.save()
        del all_ods

    all_ods = pd.concat(region_total, axis=0, sort = 'False', ignore_index=True).fillna(0)
    all_ods_val_cols = [c for c in all_ods.columns.values.tolist() if c not in ('o_region','d_region')]
    all_ods_regions = all_ods.groupby(['o_region','d_region'])[all_ods_val_cols].sum().reset_index()
    all_ods_regions.to_excel(excl_wrtr_reg,'total', index = False)
    excl_wrtr_reg.save()


if __name__ == '__main__':
    main()
