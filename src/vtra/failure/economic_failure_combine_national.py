"""Combine national-scale macroeconomic loss estimates with rerouting losses
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

import geopandas as gpd
import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from vtra.utils import *
from vtra.transport_flow_and_failure_functions import *


def main():
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    modes = ['road', 'rail']
    types = ['min', 'max']
    percentage = [10,90,100]
    single_edge = True
    multi_modal = [False,True]

    network_data_path = os.path.join(data_path,'post_processed_networks')
    shp_output_path = os.path.join(output_path, 'failure_shapefiles')
    csv_data_path = os.path.join(output_path, 'failure_results','minmax_combined_scenarios')
    econ_data_path = os.path.join(
        output_path, 'economic_failure_results', 'summarized')

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
            print ('* Loading {} transport loss and tonnage estimates'.format(modes[m]))
            for m_mod in multi_modal:
                if m_mod == True:
                    if single_edge == True:
                        file_name = 'single_edge_failures_minmax_national_{0}_{1}_percent_disrupt_multi_modal'.format(modes[m],int(perct))
                    else:
                        file_name = 'multiple_edge_failures_minmax_national_{0}_{1}_percent_disrupt_multi_modal'.format(modes[m],int(perct))
                else:
                    if single_edge == True:
                        file_name = 'single_edge_failures_minmax_national_{0}_{1}_percent_disrupt'.format(modes[m],int(perct))
                    else:
                        file_name = 'multiple_edge_failures_minmax_national_{0}_{1}_percent_disrupt'.format(modes[m],int(perct))

                if file_name + '.csv' in os.listdir(csv_data_path):
                    all_results = pd.read_csv(os.path.join(csv_data_path,file_name + '.csv'))
                    print ('* Loading {} economic loss estimates'.format(modes[m]))
                    edge_fail_ranges = []
                    for t in range(len(types)):
                        if m_mod == True:
                            if single_edge == True:
                                f_name = 'single_edge_failures_od_losses_national_{0}_{1}_{2}_percent_disrupt_multi_modal_summarized.csv'.format(modes[m], types[t],int(perct))
                            else:
                                f_name = 'multiple_edge_failures_od_losses_national_{0}_{1}_{2}_percent_disrupt_multi_modal_summarized.csv'.format(modes[m], types[t],int(perct))
                        else:
                            if single_edge == True:
                                f_name = 'single_edge_failures_od_losses_national_{0}_{1}_{2}_percent_disrupt_summarized.csv'.format(modes[m], types[t],int(perct))
                            else:
                                f_name = 'multiple_edge_failures_od_losses_national_{0}_{1}_{2}_percent_disrupt_summarized.csv'.format(modes[m], types[t],int(perct))

                        if f_name in os.listdir(econ_data_path):
                            df_path = os.path.join(
                                econ_data_path, f_name)
                            df = pd.read_csv(df_path, index_col=0)
                            df.index.names = ['edge_id']
                            df = df.reset_index()

                            df['total_losses'] = -1e6*df['total_losses']
                            df = df[df['total_losses'] > 0]
                            df.rename(columns={'total_losses': '{}_econ_loss'.format(types[t])}, inplace=True)
                        else:
                            df = all_results[['edge_id']]
                            df['{}_econ_loss'.format(types[t])] = 0

                        edge_fail_ranges.append(df)
                        del df

                    edge_impact = edge_fail_ranges[0]
                    edge_impact = pd.merge(edge_impact, edge_fail_ranges[1], how='left', on=[
                                           'edge_id']).fillna(0)
                    
                    print ('* Loading {} merging transport loss and economic loss estimates'.format(modes[m]))
                    
                    if 'min_econ_loss' in all_results.columns.values.tolist():
                        all_results.drop('min_econ_loss', axis=1, inplace=True)

                    if 'max_econ_loss' in all_results.columns.values.tolist():
                        all_results.drop('max_econ_loss', axis=1, inplace=True)

                    all_results = pd.merge( all_results, edge_impact, how='left', on=[
                                           'edge_id']).fillna(0)

                    all_results = rearrange_minmax_values(all_results)

                    all_results['min_econ_impact'] = all_results['min_tr_loss'] + all_results['min_econ_loss']
                    all_results['max_econ_impact'] = all_results['max_tr_loss'] + all_results['max_econ_loss']

                    all_results.to_csv(os.path.join(csv_data_path,file_name + '.csv'),index=False)
                    """Create network shapefiles with flows
                    """
                    print ('* Creating {} network shapefiles with failure results'.format(modes[m]))
                    shp_path = os.path.join(
                        shp_output_path,file_name + '.shp')
                    network_failure_assembly_shapefiles(all_results,gdf_edges, save_edges=True, shape_output_path=shp_path)


if __name__ == "__main__":
    main()
