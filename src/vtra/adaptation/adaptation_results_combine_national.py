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

    """
    Get the modal shares
    """
    # modes = ['road','rail','inland','coastal']
    modes = ['road']

    for m in range(len(modes)):
        mode_data_file = os.path.join(
            shp_output_path, 'weighted_edges_failures_national_{0}_2.shp'.format(modes[m]))
        G_df = gpd.read_file(mode_data_file)

        print(G_df)
        df_path = os.path.join(data_path, 'Results', 'Adaptation_results',
                               'single_edge_failures_scenarios_national_{}_adapt_options.xlsx'.format(modes[m]))
        # df_path = os.path.join(econ_paths_data,'single_edge_failures_totals_national_{0}_{1}_summarized.csv'.format(modes[m], types[t]))
        df = pd.read_excel(df_path, sheet_name='forecast_10')
        df = df[['edge_id', 'min_adapt_npv', 'max_adapt_npv', 'min_bc_ratio', 'max_bc_ratio']]
        df = df.groupby(['edge_id'])['min_adapt_npv', 'max_adapt_npv',
                                     'min_bc_ratio', 'max_bc_ratio'].min().reset_index()
        print(df)
        network_failure_assembly(
            df, modes[m], G_df, save_edges=True, output_path=shp_output_path)


if __name__ == "__main__":
    main()
