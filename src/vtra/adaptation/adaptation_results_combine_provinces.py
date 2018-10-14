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

import geopandas as gpd
import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from vtra.utils import *


def main():
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    types = ['min', 'max']

    shp_output_path = os.path.join(output_path, 'failure_shapefiles')

    province_list = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']

    """
    Path OD flow disruptions
    """
    for prn in range(len(province_list)):
        # for prn in range(0, 1):
        province = province_list[prn]
        # set all paths for all input files we are going to use
        province_name = province.replace(' ', '').lower()
        mode_data_file = os.path.join(
            shp_output_path, 'weighted_edges_commune_center_failures_{0}_5_tons.shp'.format(province_name))
        G_df = gpd.read_file(mode_data_file)

        print(G_df)
        df_path = os.path.join(data_path, 'Results', 'Adaptation_results',
                               'single_edge_failures_commune_access_scenarios_{}_5_tons_adapt_options.xlsx'.format(province_name))
        # df_path = os.path.join(econ_paths_data,'single_edge_failures_totals_national_{0}_{1}_summarized.csv'.format(modes[m], types[t]))
        df = pd.read_excel(df_path, sheet_name='forecast_10')
        df = df[['edge_id', 'min_adapt_npv', 'max_adapt_npv', 'min_bc_ratio', 'max_bc_ratio']]
        df = df.groupby(['edge_id'])['min_adapt_npv', 'max_adapt_npv',
                                     'min_bc_ratio', 'max_bc_ratio'].min().reset_index()
        print(df)
        network_failure_assembly(df, province_name, G_df,
                                 save_edges=True, output_path=shp_output_path)


if __name__ == "__main__":
    main()
