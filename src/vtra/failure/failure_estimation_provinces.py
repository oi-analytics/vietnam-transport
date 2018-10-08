"""Failure analysis of province-scale road networks
To estimate changing accessibility to commune points

Input data requirements
-----------------------
1. Correct paths to all files and correct input parameters
2. Excel sheets with results of flow mapping based on MIN-MAX generalised costs estimates:
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
3. Shapefiles
    edge_id - String/Integer/Float Edge ID
    geometry - Shapely LineString geomtry of edges
        
Results
-------
Csv sheets with results of failure analysis:

1. All failure scenarios
    edge_id - String name or list of failed edges
    origin - String node ID of Origin of disrupted OD flow
    destination - String node ID of Destination of disrupted OD flow
    o_region - String name of Province of Origin node ID of disrupted OD flow
    d_region - String name of Province of Destination node ID of disrupted OD flow
    no_access - Boolean 1 (no reroutng) or 0 (rerouting)     
    min/max_distance - Float value of estimated distance of OD journey before disruption   
    min/max_time - Float value of estimated time of OD journey before disruption     
    min/max_gcost - Float value of estimated travel cost of OD journey before disruption    
    min/max_vehicle_nums - Float value of estimated vehicles of OD journey before disruption       
    new_cost - Float value of estimated cost of OD journey after disruption 
    new_distance - Float value of estimated distance of OD journey after disruption    
    new_path - List of string edge ID's of estimated new route of OD journey after disruption   
    new_time - Float value of estimated time of OD journey after disruption      
    dist_diff - Float value of Post disruption minus per-disruption distance  
    time_diff - Float value Post disruption minus per-disruption timee  
    min/max_tr_loss - Float value of estimated change in rerouting cost
    min/max_netrev - Float values of total daily net revenues along disrupted OD pairs
    min/max_tons - Float values of total daily crop tonnages along disrupted OD pairs
    min_max_econ_impact - Float values of total daily economic impact of disrupted OD pairs

2. Min-max combined scenarios - Combined min-max results of total network impacts of each edge
    edge_id - String name or list of failed edges
    no_access - Boolean 1 (no reroutng) or 0 (rerouting)
    min/max_tr_loss - Float values of estimated change in rerouting cost
    min/max_tons - Float values of total daily tonnages along edge
    min/max_netrev - Float values of total daily net revenues along edge
    min/max_econ_impact - Float value of total daily economic impact of edge 

Shapefiles

5. Min-max combined scenarios - Combined min-max reults of total network impacts of each edge
    edge_id - String name or list of failed edges
    no_access - Boolean 1 (no reroutng) or 0 (rerouting)
    min/max_tr_loss - Float values of estimated change in rerouting cost
    min/max_tons - Float values of total daily tonnages along edge
    min/max_netrev - Float values of total daily net revenues along edge
    min/max_econ_impact - Float value of total daily economic impact of edge     
    geometry - Shapely LineString geomtry of edges

"""
import ast
import copy
import csv
import itertools
import math
import operator
import os
import sys

import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from vtra.utils import *
from vtra.transport_flow_and_failure_functions import *


def main():
    """
    Specify the paths from where you want to read and write:
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    1. Names of Provinces
        List of strings
    2. Unit weight of vehicle assumed for each mode
        List of float types 
    3. Min-max names of names of different types of attributes - paths, distance, time, cost, vehicles, tons, revenue   
        List of string types
    4. Percentage of OD flows that are assumed disrupted
        List of float type
    5. Condition on whether analysis is single failure or multiple failure
        Boolean condition True or False 

    Give the paths to the input data files:
    1. Network edges Excel and shapefiles
    2. OD flows Excel file
    3. Failure scenarios Excel file
    
    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    province_list = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
    truck_unit_wt = [5.0]

    types = ['min', 'max']
    path_types = ['min_edge_path', 'max_edge_path']
    dist_types = ['min_distance', 'max_distance']
    time_types = ['min_time', 'max_time']
    cost_types = ['min_gcost', 'max_gcost']
    vehicle_types = ['min_vehicle_nums', 'max_vehicle_nums']
    tons_types = ['min_croptons', 'max_croptons']
    rev_types = ['min_netrev', 'max_netrev']
    percentage = [100.0]
    single_edge = True


    """Give the paths to the input data files
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    network_data_excel = os.path.join(data_path,'post_processed_networks','province_roads_edges.xlsx')
    flow_paths_data = os.path.join(output_path, 'flow_mapping_paths')
    fail_scenarios_data = os.path.join(
        output_path, 'hazard_scenarios', 'province_roads_hazard_intersections.xlsx')

    """Specify the output files and paths to be created 
    """
    shp_output_path = os.path.join(output_path, 'failure_shapefiles')
    if os.path.exists(shp_output_path) == False:
        os.mkdir(shp_output_path)

    fail_output_path = os.path.join(output_path, 'failure_results')
    if os.path.exists(fail_output_path) == False:
        os.mkdir(fail_output_path)

    all_fail_scenarios = os.path.join(fail_output_path,'all_fail_scenarios')
    if os.path.exists(all_fail_scenarios) == False:
        os.mkdir(all_fail_scenarios)

    isolated_ods = os.path.join(fail_output_path,'isolated_od_scenarios')
    if os.path.exists(isolated_ods) == False:
        os.mkdir(isolated_ods)

    rerouting = os.path.join(fail_output_path,'rerouting_scenarios')
    if os.path.exists(rerouting) == False:
        os.mkdir(rerouting)

    minmax_combine = os.path.join(fail_output_path,'minmax_combined_scenarios')
    if os.path.exists(minmax_combine) == False:
        os.mkdir(minmax_combine)

    """
    Path OD flow disruptions
    """
    for prn in range(len(province_list)):
        province = province_list[prn]
        province_name = province.replace(' ', '').lower()

        """Load igraph network and GeoDataFrame
        """
        print ('* Loading {} igraph network and GeoDataFrame'.format(province))
        G_df = pd.read_excel(network_data_excel,sheet_name = province_name,encoding='utf-8')
        gdf_edges = gpd.read_file(os.path.join(network_data_path,'{}_roads_edges.shp'.format(province_name)),encoding='utf-8')
        gdf_edges = gdf_edges[['edge_id','geometry']]

        """Create failure scenarios
        """
        print ('* Creating {} failure scenarios'.format(province))
        fail_df = pd.read_excel(fail_scenarios_data, sheet_name=province_name)
        ef_sc_list = edge_failure_sampling(fail_df,'edge_id')

        for perct in percentage:
            for tr_wt in truck_unit_wt:
                """Load flow paths
                """
                print ('* Loading {} flow paths'.format(province))
                flow_excel_file = 'province_roads_commune_center_access_flow_paths_{}_percent.xlsx'.format(int(perct))
                flow_df = pd.read_excel(
                    os.path.join(flow_paths_data,flow_excel_file), 
                    sheet_name = province_name + '_{}_tons'.format(int(tr_wt)))

                """Perform failure analysis
                """
                edge_fail_ranges = []
                for t in range(len(types)):
                    print ('* Performing  {} {} failure analysis'.format(types[t],province))
                    ef_list = []
                    for fail_edge in ef_sc_list:
                        if isinstance(fail_edge,list) == False:
                            fail_edge = [fail_edge]
                        ef_dict = igraph_scenario_edge_failures(
                            G_df, fail_edge, flow_df, tr_wt,path_types[t], tons_types[t], cost_types[t], time_types[t])

                        if ef_dict:
                            ef_list += ef_dict

                        print('Done with province {0} edge {1} type {2}'.format(
                            province, fail_edge, types[t]))

                    df = pd.DataFrame(ef_list)

                    print ('* Assembling {} {} failure results'.format(types[t],province))
                    ic_cols = [rev_types[t],tons_types[t]]
                    select_cols = ['origin', 'destination', dist_types[t], time_types[t],
                               cost_types[t], vehicle_types[t]] + ic_cols
                    flow_df_select = flow_df[select_cols]
                    flow_df_select = merge_failure_results(flow_df_select,df,tons_types[t],
                        dist_types[t],time_types[t],cost_types[t],vehicle_types[t],changing_tonnages=False)

                    del df
                    
                    tr_loss = '{}_tr_loss'.format(types[t])
                    econ_impact = '{}_econ_impact'.format(types[t])

                    flow_df_select.rename(columns={'tr_loss': tr_loss}, inplace=True)
                    flow_df_select[econ_impact] = flow_df_select['no_access']*flow_df_select[rev_types[t]] + \
                        flow_df_select[tr_loss]
                    
                    if single_edge == True:
                        file_name = 'single_edge_failures_all_{0}_{1}_tons_{2}_percent_disrupt.csv'.format(province_name, int(tr_wt),int(perct))
                    else:
                        file_name = 'multiple_edge_failures_all_{0}_{1}_tons_{2}_percent_disrupt.csv'.format(province_name, int(tr_wt),int(perct))

                    df_path = os.path.join(all_fail_scenarios,file_name)
                    flow_df_select.to_csv(df_path, index=False)
                    
                    print ('* Combining {} {} rerouting and isolation failure results'.format(types[t],province))
                    edge_impact = flow_df_select[['edge_id','no_access',tr_loss,tons_types[t],rev_types[t],econ_impact]]
                    edge_impact = edge_impact.groupby(
                        ['edge_id','no_access'])[tr_loss,tons_types[t],rev_types[t],econ_impact].sum().reset_index()
                    edge_fail_ranges.append(edge_impact)

                print ('* Assembling {} min-max failure results'.format(province))
                edge_impact = edge_fail_ranges[0]
                edge_impact = pd.merge(edge_impact, edge_fail_ranges[1], how='left', on=[
                                       'edge_id','no_access']).fillna(0)
                
                del edge_fail_ranges
                if single_edge == True:
                    file_name = 'single_edge_failures_minmax_{0}_{1}_tons_{2}_percent_disrupt'.format(province_name, int(tr_wt),int(perct))
                else:
                    file_name = 'multiple_edge_failures_minmax_{0}_{1}_tons_{2}_percent_disrupt'.format(province_name, int(tr_wt),int(perct))

                df_path = os.path.join(minmax_combine,file_name + '.csv')
                edge_impact.to_csv(df_path, index=False)

                """Create network shapefiles with flows
                """
                print ('* Creating {} network shapefiles with flows'.format(province))
                shp_path = os.path.join(
                    shp_output_path,file_name + '.shp')
                network_failure_assembly_shapefiles(edge_impact,gdf_edges, save_edges=True, shape_output_path=shp_path)


if __name__ == "__main__":
    main()
