# -*- coding: utf-8 -*-
"""
Python script to assign commodity flows on the road network

TODO DONT THINK THIS SCRIPT IS NEEDED
"""
import ast
import copy
import csv
import itertools
import math
import operator
import os
import sys
from collections import Counter

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

    1. Names of modes
        List of strings
    2. Unit weight of vehicle assumed for each mode
        List of float types 
    3. Range of usage factors for each mode to represent uncertainty in cost estimations
        List of tuples of float types 
    4. Min-max names of names of different types of attributes - paths, distance, time, cost, vehicles, tons   
        List of string types
    5. Names of commodity/industry columns for which min-max tonnage column names already exist
        List of string types
    6. Percentage of OD flows that are assumed disrupted
        List of float type
    7. Condition on whether analysis is single failure or multiple failure
        Boolean condition True or False 

    Give the paths to the input data files:
    
    1. Network edges Excel and shapefiles
    2. OD flows Excel file
    3. Costs of modes Excel file 
    4. Road properties Excel file
    5. Failure scenarios Excel file
    
    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    modes = ['road', 'rail','inland']

    types = ['min', 'max']
    percentage = [10,100.0]
    single_edge = True
    
    """Give the paths to the input data files
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    all_fail_scenarios = os.path.join(fail_output_path,'all_fail_scenarios')

    """Specify the output files and paths to be created 
    """
    shp_output_path = os.path.join(output_path, 'failure_shapefiles')
    if os.path.exists(shp_output_path) == False:
        os.mkdir(shp_output_path)

    minmax_combine = os.path.join(fail_output_path,'minmax_combined_scenarios')
    if os.path.exists(minmax_combine) == False:
        os.mkdir(minmax_combine)

    """Create theee multi-modal networks
    """
    print ('* Creating multi-modal networks') 
    G_multi_df = []
    for m in range(len(modes)):
        """Load mode igraph network and GeoDataFrame
        """
        print ('* Loading {} igraph network and GeoDataFrame'.format(modes[m]))
        gdf_edges = gpd.read_file(os.path.join(network_data_path,'{}_edges.shp'.format(modes[m])),encoding='utf-8')
        G_multi_df.append(gdf_edges[['edge_id','geometry']])

        del gdf_edges

    """
    Get the modal shares
    """
    fail_mode = ['road']

    for m in range(len(fail_mode)):
        for perct in percentage:
            edge_fail_ranges = []
            for t in range(len(types)):
                mode_dict = {'road': {}, 'rail': {}, 'waterways': {}}
                if single_edge == True:
                    file_name = 'single_edge_failures_all_national_{0}_{1}_{2}_percent_disrupt_multi_modal.csv'.format(modes[m], types[t],int(perct))
                else:
                    file_name = 'multiple_edge_failures_all_national_{0}_{1}_{2}_percent_disrupt_multi_modal.csv'.format(modes[m], types[t],int(perct))

                df_path = os.path.join(all_fail_scenarios,file_name)
                flow_df_select.read_csv(df_path).fillna(0)

                flow_df_select = flow_df_select[[
                    'origin', 'destination', 'new_path', '{}_tons'.format(types[t])]]
                print('Old number of paths', len(flow_df_select.index))
                flow_df_select = flow_df_select.drop_duplicates(
                    subset=['origin', 'destination', 'new_path'], keep='first')
                print('New Number of paths', len(flow_df_select.index))
                for iter_, row in flow_df_select.iterrows():
                    new_path = ast.literal_eval(row['new_path'])
                    if new_path:
                        for p in new_path:
                            if 'road' in p:
                                if p in (mode_dict['road'].keys()):
                                    mode_dict['road'][p].append(
                                        (row['origin'], row['destination'], row['{}_tons'.format(types[t])]))
                                else:
                                    mode_dict['road'][p] = [
                                        (row['origin'], row['destination'], row['{}_tons'.format(types[t])])]
                            elif 'rail' in p:
                                if p in (mode_dict['rail'].keys()):
                                    mode_dict['rail'][p].append(
                                        (row['origin'], row['destination'], row['{}_tons'.format(types[t])]))
                                else:
                                    mode_dict['rail'][p] = [
                                        (row['origin'], row['destination'], row['{}_tons'.format(types[t])])]
                            elif 'water' in p:
                                if p in (mode_dict['waterways'].keys()):
                                    mode_dict['waterways'][p].append(
                                        (row['origin'], row['destination'], row['{}_tons'.format(types[t])]))
                                else:
                                    mode_dict['waterways'][p] = [
                                        (row['origin'], row['destination'], row['{}_tons'.format(types[t])])]
                del flow_df_select
                edge_impact = []
                for key, values in mode_dict.items():
                    values_tup_list = []
                    for a, b in values.items():
                        c = list(set(b))
                        v = sum([z for (x, y, z) in c])
                        values_tup_list.append((a, v))
                    edge_impact.append(pd.DataFrame(values_tup_list, columns=[
                                       'edge_id', '{}_tons'.format(types[t])]))

                edge_impact = pd.concat(edge_impact, axis=0, sort='False', ignore_index=True)

                edge_fail_ranges.append(edge_impact)
                del edge_impact

            edge_impact = edge_fail_ranges[0]
            edge_impact = pd.merge(edge_impact, edge_fail_ranges[1], how='left', on=[
                                   'edge_id']).fillna(0)
            df_path = os.path.join(minmax_combine,
                                   'single_edge_failures_transfers_national_{0}_{1}_percent_shift.csv'.format(fail_mode[m],int(perct)))
            edge_impact.to_csv(df_path, index=False)

            # network_failure_assembly_shapefiles(edge_impact, fail_mode[m], [
            #                          'road', 'rail', 'water'], G_multi_df, save_edges=True, output_path=shp_output_path)


if __name__ == "__main__":
    main()
