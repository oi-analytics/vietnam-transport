"""
Purpose
-------

Create province scale OD matrices between roads connecting villages to nearest communes: 
    - Net revenue estimates of commune villages
    - IFPRI crop data at 1km resolution

Input data requirements
-----------------------

1. Correct paths to all files and correct input parameters
2. Geotiff files with IFPRI crop data:
    - tons - Float values of production tonnage at each grid cell
    - geometry - Raster grid cell geometry

3. Shapefile of RiceAtlas data:
    - month production columns - tonnage of rice for each month
    - geometry - Shapely Polygon geometry of Provinces

4. Shapefile of Provinces
    - od_id - Integer Province ID corresponding to OD ID
    - name_eng - String name of Province in English
    - geometry - Shapely Polygon geometry of Provinces

5. Shapefile of Communes
    - population - Float values of populations in Communes
    - nfrims - Float values of number of firms in Provinces
    - netrevenue - Float values of Net Revenue in Provinces
    - argi_prop - Float values of proportion of agrivculture firms in Provinces
    - geometry - Shapely Polygon geometry of Communes

6. Shapefiles of network nodes
    - node_id - String node ID
    - geometry - Shapely point geometry of nodes

7. Shapefiles of network edges
    - vehicle_co - Count of vehiles only for roads
    - geometry - Shapely LineString geometry of edges 

8. Shapefiles of Commune center points
    - object_id - Integer ID of point
    - geometry - Shapely point geometry of points

9. Shapefiles of Village center points
    - object_id - Integer ID of points
    - geometry - Shapely point geometry of points

Results
-------

1. Excel workbook with sheet of mode-wise and total OD flows
    - origin - String node ID of origin node
    - destination - String node ID of destination node
    - crop_names - Float values of daily tonnages of IFPRI crops (except rice) between OD nodes
    - min_rice - Float values of minimum daily tonnages of rice between OD nodes
    - max_rice - Float values of maximum daily tonnages of rice between OD nodes
    - min_croptons - Float values of minimum daily tonnages of crops between OD nodes
    - max_croptons - Float values of maximum daily tonnages of crops between OD nodes
    - min_agrirev - Float value of Minimum daily revenue of agriculture firms between OD nodes
    - max_agrirev - Float value of Maximum daily revenue of agriculture firms between OD nodes
    - min_noagrirev - Float value of Minimum daily revenue of non-agriculture firms between OD nodes
    - max_noagrirev - Float value of Maximum daily revenue of non-agriculture firms between OD nodes
    - min_netrev - Float value of Minimum daily revenue of all firms between OD nodes
    - max_netrev - Float value of Maximum daily revenue of all firms between OD nodes

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
from shapely.geometry import Point
from vtra.utils import *


def netrev_od_pairs(start_points, end_points):
    """
    Assign crop tonnages to OD pairs 

    Parameters
    ---------
    - start_points - GeoDataFrame of start points for Origins
    - end_points - GeoDataFrame of potential end points for Destinations

    Outputs
    -------
    od_pairs_df - Pandas DataFrame with columns:
        - origin - Origin node ID
        - destination - Destination node ID
        - netrev_argi - Net revenue of agriculture firms
        - netrev_noargi - Net revenue of non-agriculture firms
    """
    save_paths = []
    for iter_, place in start_points.iterrows():
        try:
            closest_center = end_points.loc[end_points['OBJECTID']
                                            == place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]

            save_paths.append(
                (closest_center, place['NEAREST_G_NODE'], place['netrev_agri'], place['netrev_noagri']))
        except:
            print(iter_)

    od_pairs_df = pd.DataFrame(
        save_paths, columns=['origin', 'destination', 'netrev_agri', 'netrev_noagri'])
    od_pairs_df = od_pairs_df.groupby(['origin', 'destination'])[
        'netrev_agri', 'netrev_noagri'].sum().reset_index()

    return od_pairs_df


def crop_od_pairs(start_points, end_points, crop_name):
    """
    Assign crop tonnages to OD pairs 

    Parameters
    ----------
    - start_points - GeoDataFrame of start points for Origins
    - end_points - GeoDataFrame of potential end points for Destinations
    - crop_name - String name of crop

    Outputs
    -------
    od_pairs_df - Pandas DataFrame wit columns:
        - origin - Origin node ID
        - destination - Destination node ID
        - crop - Tonnage values for the named crop
        - netrev_argi - Daily Net revenue of agriculture firms in USD
        - netrev_noargi - Daily Net revenue of non-agriculture firms in USD
    """
    save_paths = []
    for iter_, place in start_points.iterrows():
        try:
            closest_center = end_points.loc[end_points['OBJECTID']
                                            == place['NEAREST_C_CENTER']]['NEAREST_G_NODE'].values[0]

            save_paths.append((closest_center, place['NEAREST_G_NODE'], place['tons']))
        except:
            print(iter_)

    
    od_pairs_df = pd.DataFrame(save_paths, columns=['origin', 'destination', crop_name])
    od_pairs_df = od_pairs_df.groupby(['origin', 'destination'])[crop_name].sum().reset_index()

    return od_pairs_df

def assign_monthly_tons_crops(x,rice_prod_file,crop_month_fields,province,x_cols):
    """
    Assign crop tonnages to OD pairs 

    Parameters
    ---------
    - x - Pandas DataFrame of values
    - rice_prod_file - Shapefile of RiceAltas monthly production value
    - crop_month_fields - Lsit of strings of month columns in Rice Atlas shapefile
    - province - Stirng name of province
    - x_cols - List of string names of crops

    Outputs
    -------
    - min_croptons - Float value of Minimum daily tonnages of crops 
    - max_croptons - Float value of Maximum daily tonnages of crops 
    """
    # find the crop production months for the province
    rice_prod_months = gpd.read_file(rice_prod_file)
    rice_prod_months = rice_prod_months.loc[rice_prod_months.SUB_REGION == province]
    rice_prod_months = rice_prod_months[crop_month_fields].values.tolist()
    rice_prod_months = np.array(rice_prod_months[0])/sum(rice_prod_months[0])
    rice_prod_months = rice_prod_months[rice_prod_months > 0]
    rice_prod_months = rice_prod_months.tolist()
    
    min_croptons = 0
    max_croptons = 0
    for x_name in x_cols:
        if x_name == 'rice':
            min_croptons += (1.0*min(rice_prod_months)*x[x_name])/30.0
            max_croptons += (1.0*max(rice_prod_months)*x[x_name])/30.0
        else:
            min_croptons += (1.0*x[x_name])/365.0
            max_croptons += (1.0*x[x_name])/365.0

    return min_croptons, max_croptons


def assign_io_rev_costs_crops(x, cost_dataframe,rice_prod_file,crop_month_fields,province, x_cols, ex_rate):
    """
    Assign crop tonnages to daily net revenues

    Parameters
    ---------
    - x - Pandas DataFrame of values
    - cost_dataframe - Pandas DataFrame of conversion of tonnages to net revenues
    - rice_prod_file - Shapefile of RiceAltas monthly production value
    - province - Stirng name of province
    - x_cols - List of string names of crops
    - ex_rate - Exchange rate from VND millions to USD

    Outputs
    -------
    - min_croprev - Float value of Minimum daily revenue of crops 
    - max_croprev - Float value of Maximum daily revenue of crops 
    """
    # find the crop production months for the province
    rice_prod_months = gpd.read_file(rice_prod_file)
    rice_prod_months = rice_prod_months.loc[rice_prod_months.SUB_REGION == province]
    rice_prod_months = rice_prod_months[crop_month_fields].values.tolist()
    rice_prod_months = np.array(rice_prod_months[0])/sum(rice_prod_months[0])
    rice_prod_months = rice_prod_months[rice_prod_months > 0]
    rice_prod_months = rice_prod_months.tolist()

    min_croprev = 0
    max_croprev = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if cost_param.crop_code in x_cols:
            if cost_param.crop_code == 'rice':
                min_croprev += (1.0*min(rice_prod_months)*ex_rate*cost_param.est_net_rev *
                                (x[cost_param.crop_code]/cost_param.tot_tons))/30.0
                max_croprev += (1.0*max(rice_prod_months)*ex_rate*cost_param.est_net_rev *
                                (x[cost_param.crop_code]/cost_param.tot_tons))/30.0
            else:
                min_croprev += 1.0/365.0 * \
                    (ex_rate*cost_param.est_net_rev *
                     (x[cost_param.crop_code]/cost_param.tot_tons))
                max_croprev += 1.0/365.0 * \
                    (ex_rate*cost_param.est_net_rev *
                     (x[cost_param.crop_code]/cost_param.tot_tons))

    return min_croprev, max_croprev

def netrevenue_values_to_province_od_nodes(province_ods_df,prov_communes,commune_sindex,netrevenue,
    n_firms,agri_prop,prov_pop,prov_pop_sindex,nodes,sindex_nodes,prov_commune_center,
    sindex_commune_center,node_id,object_id,exchange_rate):
    """
    Assign commune level netrevenue values to OD nodes in provinces
        - Based on finding nearest nodes to village points with netrevenues as Origins
        - And finding nearest commune centers as Destinations

    Parameters
    ----------
    - province_ods_df - List of lists of Pandas dataframes
    - prov_communes - GeoDataFrame of commune level statistics
    - commune_sindex - Spatial index of communes
    - netrevenue - String name of column for netrevenue of communes in VND millions
    - nfirm - String name of column for numebr of firms in communes
    - agri_prop - Stirng name of column for proportion of agriculture firms in communes
    - prov_pop - GeoDataFrame of population points in Province
    - prov_pop_sindex - Spatial index of population points in Province
    - nodes - GeoDataFrame of province road nodes
    - sindex_nodes - Spatial index of province road nodes 
    - prov_commune_center - GeoDataFrame of province commune center points 
    - sindex_commune_center - Spatial index of commune center points
    - node_id - String name of Node ID column
    - object_id - String name of commune ID column
    - exchange_rate - Float value for exchange rate from VND million to USD 

    Outputs
    -------
    province_ods_df - List of Lists of Pandas dataframes with columns:
        - origin - Origin node ID
        - destination - Destination node ID
        - netrev_argi - Net revenue of agriculture firms
        - netrev_noargi - Net revenue of non-agriculture firms
    """

    # create new column in prov_communes with amount of villages
    prov_communes['n_villages'] = prov_communes.geometry.apply(
        lambda x: count_points_in_polygon(x, prov_pop_sindex))
    prov_communes['netrev_village'] = exchange_rate * \
        (prov_communes[netrevenue]*prov_communes[n_firms])/prov_communes['n_villages']
    # also get the net revenue of the agriculture sector which is called nongnghiep
    prov_communes['netrev_village_agri'] = 1.0/365.0 * \
        (prov_communes[agri_prop]*prov_communes['netrev_village'])
    prov_communes['netrev_village_noagri'] = 1.0/365.0 * \
        (prov_communes['netrev_village'] - prov_communes['netrev_village_agri'])

    # give each village a net revenue based on average per village in commune
    prov_pop['netrev_agri'] = prov_pop.geometry.apply(lambda x: extract_value_from_gdf(
        x, commune_sindex, prov_communes, 'netrev_village_agri'))
    prov_pop['netrev_noagri'] = prov_pop.geometry.apply(lambda x: extract_value_from_gdf(
        x, commune_sindex, prov_communes, 'netrev_village_noagri'))

    # get nearest node in network for all start and end points
    prov_pop['NEAREST_G_NODE'] = prov_pop.geometry.apply(
        lambda x: get_nearest_node(x, sindex_nodes, nodes, node_id))

    prov_pop['NEAREST_C_CENTER'] = prov_pop.geometry.apply(
        lambda x: get_nearest_node(x, sindex_commune_center, prov_commune_center, object_id))

    # find all OD pairs of the revenues
    netrev_ods = netrev_od_pairs(prov_pop, prov_commune_center)
    province_ods_df.append(netrev_ods)

    return province_ods_df

def crop_values_to_province_od_nodes(province_ods_df,province_geom,calc_path,
    crop_data_path,crop_names,nodes,sindex_nodes,prov_commune_center,sindex_commune_center,node_id,object_id):
    """
    Assign IFPRI crop values to OD nodes in provinces
        - Based on finding nearest nodes to crop production sites as Origins
        - And finding nearest commune centers as Destinations

    Parameters
    ----------
    - province_ods_df - List of lists of Pandas dataframes
    - province_geom - Shapely Geometry of province 
    - calc_path - Path to store intermediary calculations 
    - crop_data_path - Path to crop datasets
    - crop_names - List of string of crop names in IFPRI datasets
    - nodes - GeoDataFrame of province road nodes
    - sindex_nodes - Spatial index of province road nodes 
    - prov_commune_center - GeoDataFrame of province commune center points 
    - sindex_commune_center - Spatial index of commune center points
    - node_id - String name of Node ID column
    - object_id - String name of commune ID column 

    Outputs
    -------
    province_ods_df - List of Lists of Pandas dataframes with columns:
        - origin - Origin node ID
        - destination - Destination node ID
        - crop - Tonnage values for the named crop
    """
    # all the crop OD pairs
    for file in os.listdir(crop_data_path):
        if file.endswith(".tif") and 'spam_p' in file.lower().strip():
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

            # clip all to province
            prov_crop = gdf_geom_clip(crop_points, province_geom)

            if len(prov_crop.index) > 0:
                prov_crop_sindex = prov_crop.sindex
                prov_crop['NEAREST_G_NODE'] = prov_crop.geometry.apply(
                    lambda x: get_nearest_node(x, sindex_nodes, nodes, node_id))
                prov_crop['NEAREST_C_CENTER'] = prov_crop.geometry.apply(
                    lambda x: get_nearest_node(x, sindex_commune_center, prov_commune_center, object_id))

                crop_ods = crop_od_pairs(prov_crop, prov_commune_center, crop_name)
                province_ods_df.append(crop_ods)

    return province_ods_df

def main():
    """
    1. Specify the paths from where you to read and write:
        - Input data
        - Intermediate calcuations data
        - Output results
    
    2. Supply input data and parameters
        - Names of the three Provinces: ['Lao Cai', 'Binh Dinh', 'Thanh Hoa'] 
        - Exchange rate to convert 2012 Net revenue in million VND values to USD in 2016
        - Names of crops in IFPRI crop data
        - Names of months in Rice Atlas data
        - Name of column for netrevenue of communes in VND millions
        - Name of column for numebr of firms in communes
        - Name of column for proportion of agriculture firms in communes
        - Name of Node ID column
        - Name of commune ID column

    3. Give the paths to the input data files:
        - Network nodes files
        - IFPRI crop data files
        - Rice Altas data shapefile
        - Province boundary and stats data shapefile 
        - Commune boundary and stats data shapefile
        - Population points shapefile for locations of villages
        - Commune center points shapefile 
    
    4. Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    province_list = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
    exchange_rate = 1.05*(1000000/21000)
    crop_names = ['rice', 'cash', 'cass', 'teas',
                  'maiz', 'rubb', 'swpo', 'acof', 'rcof', 'pepp']
    crop_month_fields = ['P_Jan', 'P_Feb', 'P_Mar', 'P_Apr', 'P_May',
                        'P_Jun', 'P_Jul', 'P_Aug', 'P_Sep', 'P_Oct', 'P_Nov', 'P_Dec']
    netrevenue = 'netrevenue'
    n_firms = 'nfirm'
    agri_prop = 'nongnghiep'
    node_id = 'NODE_ID'
    object_id = 'OBJECTID'

    """Give the paths to the input data files:
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    crop_data_path = os.path.join(data_path, 'Agriculture_crops', 'crop_data')
    rice_month_file = os.path.join(data_path, 'rice_atlas_vietnam', 'rice_production.shp')
    province_path = os.path.join(data_path, 'Vietnam_boundaries',
                                'boundaries_stats', 'province_level_stats.shp')
    commune_path = os.path.join(data_path, 'Vietnam_boundaries',
                                'boundaries_stats', 'commune_level_stats.shp')
    population_points_in = os.path.join(
        data_path, 'Points_of_interest', 'population_points.shp')
    commune_center_in = os.path.join(
            data_path, 'Points_of_interest', 'commune_committees_points.shp')

    """Specify the output files and paths to be created  
    """
    output_dir = os.path.join(output_path, 'flow_ods')
    if os.path.exists(output_dir) == False:
        os.mkdir(output_dir)

    flow_output_excel = os.path.join(
        output_dir, 'province_roads_commune_center_flow_ods.xlsx')
    excl_wrtr = pd.ExcelWriter(flow_output_excel)
   
    """Start the province OD allocations 
    """
    for prn in range(len(province_list)):
        province = province_list[prn]
        province_name = province.replace(' ', '').lower()

        # load provinces and get geometry of the right province
        provinces = gpd.read_file(province_path)
        provinces = provinces.to_crs({'init': 'epsg:4326'})
        province_geom = provinces.loc[provinces.name_eng == province].geometry.values[0]

        # clip all the populations to the province
        prov_pop = gdf_clip(population_points_in, province_geom)
        # create sindex of all villages to count number of villages in commune
        prov_pop_sindex = prov_pop.sindex

        # clip all the commune centers to the province
        prov_commune_center = gdf_clip(commune_center_in, province_geom)
        if object_id not in prov_commune_center.columns.values.tolist():
            prov_commune_center[object_id] = prov_commune_center.index

        sindex_commune_center = prov_commune_center.sindex

        # clip all the communes to the province
        prov_communes = gdf_clip(commune_path, province_geom)
        commune_sindex = prov_communes.sindex

        # load nodes of the network
        nodes_in = os.path.join(network_data_path, '{}_roads_nodes.shp'.format(province_name))
        nodes = gpd.read_file(nodes_in)
        nodes = nodes.to_crs({'init': 'epsg:4326'})
        sindex_nodes = nodes.sindex

        province_ods_df = []
        prov_commune_center['NEAREST_G_NODE'] = prov_commune_center.geometry.apply(
            lambda x: get_nearest_node(x, sindex_nodes, nodes, node_id))
        
        """Assign revenue values for each village to nearest road nodes
        And commune center point to nearest road nodes
        For Net Revenue OD pairs
        """
        print ('* Assigning revenue OD values for each village in {}'.format(province))
        province_ods_df = netrevenue_values_to_province_od_nodes(
                                province_ods_df,prov_communes,commune_sindex,netrevenue,n_firms,
                                agri_prop,prov_pop,prov_pop_sindex,nodes,sindex_nodes,
                                prov_commune_center,sindex_commune_center,
                                node_id,object_id,exchange_rate)

        """Get crop values and assign to the nearest road nodes
        And assign commune centers to nearest road nodes
        For crop OD pairs 
        """
        print ('* Getting crop OD values in {}'.format(province))
        province_ods_df = crop_values_to_province_od_nodes(
                                province_ods_df,province_geom,calc_path,
                                crop_data_path,crop_names,nodes,sindex_nodes,
                                prov_commune_center,sindex_commune_center,
                                node_id,object_id)

        """Combine the Net Revenue abd Crop OD results
        """
        print ('* Combining OD values in {}'.format(province))
        # Get totals across all crops
        all_ods = pd.concat(province_ods_df, axis=0, sort='False', ignore_index=True).fillna(0)

        all_ods_crop_cols = [c for c in all_ods.columns.values.tolist() if c in crop_names]
        all_ods['crop_tot'] = all_ods[all_ods_crop_cols].sum(axis=1)


        all_ods_val_cols = [c for c in all_ods.columns.values.tolist()
                            if c not in ('origin', 'destination')]
        all_ods = all_ods.groupby(['origin', 'destination'])[
            all_ods_val_cols].sum().reset_index()

        # Find minimum and maximum crop daily tonnages
        all_ods['croptons'] = all_ods.apply(lambda x: assign_monthly_tons_crops(
            x, rice_month_file,crop_month_fields,province, all_ods_crop_cols), axis=1)
        all_ods[['min_croptons', 'max_croptons']] = all_ods['croptons'].apply(pd.Series)
        all_ods.drop('croptons', axis=1, inplace=True)

        # Translate crop tonnages to netrevenues and compared with max netrevenue of firms
        cost_values_df = pd.read_excel(os.path.join(
            crop_data_path, 'crop_unit_costs.xlsx'), sheet_name='io_rev')
        all_ods['croprev'] = all_ods.apply(lambda x: assign_io_rev_costs_crops(
            x, cost_values_df,rice_month_file,crop_month_fields,province, 
            all_ods.columns.values.tolist(), exchange_rate), axis=1)
        all_ods[['min_agrirev', 'max_croprev']] = all_ods['croprev'].apply(pd.Series)
        all_ods.drop('croprev', axis=1, inplace=True)
        all_ods['max_agrirev'] = all_ods[['max_croprev', 'netrev_agri']].max(axis=1)
        all_ods.drop(['max_croprev', 'netrev_agri'], axis=1, inplace=True)

        all_ods['min_netrev'] = all_ods['min_agrirev'] + all_ods['netrev_noagri']
        all_ods['max_netrev'] = all_ods['max_agrirev'] + all_ods['netrev_noagri']

        print ('* Writing {} values to Excel'.format(province))
        all_ods.to_excel(excl_wrtr, province_name, index=False)
        excl_wrtr.save()

if __name__ == '__main__':
    main()
