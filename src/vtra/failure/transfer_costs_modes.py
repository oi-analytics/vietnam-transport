# -*- coding: utf-8 -*-
"""Assign commodity flows on the road network
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
import geopandas as gpd
from vtra.utils import *
from vtra.transport_flow_and_failure_functions import *

def main():
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    modes = ['road', 'rail']
    types = ['min', 'max']
    percentage = [10]
    single_edge = True

    network_data_path = os.path.join(data_path,'post_processed_networks')
    shp_output_path = os.path.join(output_path, 'failure_shapefiles')
    csv_data_path = os.path.join(output_path, 'failure_results','minmax_combined_scenarios')

    """
    Get the modal shares
    """

    for m in range(len(modes)):
        """Load mode igraph network and GeoDataFrame
        """
        print ('* Loading {} network GeoDataFrame'.format(modes[m]))
        gdf_edges = gpd.read_file(os.path.join(network_data_path,'{}_edges.shp'.format(modes[m])),encoding='utf-8')
        gdf_edges = gdf_edges[['edge_id','geometry']]

        for perct in percentage:
            """Load flow paths
            """
            if perct < 100:
                if single_edge == True:
                    file_name_1 = 'single_edge_failures_minmax_national_{0}_{1}_percent_disrupt_multi_modal'.format(modes[m],int(perct))
                    file_name_2 = 'single_edge_failures_minmax_national_{0}_{1}_percent_disrupt'.format(modes[m],int(100-perct))
                    file_name_3 = 'single_edge_failures_minmax_national_{0}_{1}_percent_disrupt_multi_modal'.format(modes[m],int(100))
                else:
                    file_name_1 = 'multiple_edge_failures_minmax_national_{0}_{1}_percent_disrupt_multi_modal'.format(modes[m],int(perct))
                    file_name_2 = 'multiple_edge_failures_minmax_national_{0}_{1}_percent_disrupt'.format(modes[m],int(100-perct))
                    file_name_3 = 'multiple_edge_failures_minmax_national_{0}_{1}_percent_disrupt_multi_modal'.format(modes[m],int(100))

                if file_name_1 + '.csv' in os.listdir(csv_data_path) and file_name_2 + '.csv'  in os.listdir(csv_data_path) and file_name_3 + '.csv' in os.listdir(csv_data_path):
                    edge_impact_1 = pd.read_csv(os.path.join(csv_data_path,file_name_1 + '.csv')).fillna(0)
                    edge_impact_1 = edge_impact_1[['edge_id','min_tr_loss','max_tr_loss']]
                    edge_impact_1.rename(columns={'min_tr_loss': 'min_tr_loss_1','max_tr_loss': 'max_tr_loss_1'}, inplace=True)

                    edge_impact_2 = pd.read_csv(os.path.join(csv_data_path,file_name_2 + '.csv')).fillna(0)
                    edge_impact_2 = edge_impact_2[['edge_id','min_tr_loss','max_tr_loss']]
                    edge_impact_2.rename(columns={'min_tr_loss': 'min_tr_loss_2','max_tr_loss': 'max_tr_loss_2'}, inplace=True)

                    edge_impact = pd.merge(edge_impact_1, edge_impact_2, on=[
                        'edge_id'], how='left').fillna(0)

                    edge_impact['min_tr_loss'] = edge_impact['min_tr_loss_1'] + edge_impact['min_tr_loss_2']
                    edge_impact['max_tr_loss'] = edge_impact['max_tr_loss_1'] + edge_impact['max_tr_loss_2']
                    edge_impact.drop(['min_tr_loss_1', 'min_tr_loss_2','max_tr_loss_1', 'max_tr_loss_2'], axis=1, inplace=True)

                    edge_impact_3 = pd.read_csv(os.path.join(csv_data_path,file_name_3 + '.csv')).fillna(0)
                    edge_impact_3 = edge_impact_3[['edge_id','min_econ_loss','max_econ_loss']]
                    edge_impact = pd.merge(edge_impact, edge_impact_3, on=[
                        'edge_id'], how='left').fillna(0)

                    edge_impact['min_econ_impact'] = edge_impact['min_tr_loss'] + edge_impact['min_econ_loss']
                    edge_impact['max_econ_impact'] = edge_impact['max_tr_loss'] + edge_impact['max_econ_loss']

                    if single_edge == True:
                        file_name = 'single_edge_failures_minmax_national_{0}_{1}_percent_modal_shift'.format(modes[m],int(perct))    
                    else:
                        file_name = 'multiple_edge_failures_minmax_national_{0}_{1}_percent_modal_shift'.format(modes[m],int(perct))

                    edge_impact = rearrange_minmax_values(edge_impact)
                    
                    edge_impact.to_csv(os.path.join(csv_data_path,file_name + '.csv'),index=False)

                    print ('* Creating {} network shapefiles with failure results'.format(modes[m]))
                    
                    shp_path = os.path.join(shp_output_path,file_name + '.shp')
                    network_failure_assembly_shapefiles(edge_impact,gdf_edges, save_edges=True, shape_output_path=shp_path)


if __name__ == "__main__":
    main()
