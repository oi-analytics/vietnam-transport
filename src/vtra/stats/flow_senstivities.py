"""Summarise hazard data

Get OD data and process it
"""
import ast
import itertools
import math
import operator
import os
import subprocess
import sys

import geopandas as gpd
import igraph as ig
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi
from shapely.geometry import Point, Polygon
from vtra.utils import *


def main():
    data_path, calc_path, output_path, figure_path = load_config()['paths']['data'], load_config(
    )['paths']['calc'], load_config()['paths']['output'], load_config()['paths']['figures']

    # Get the modal shares
    modes_cols = ['road', 'rail', 'inland', 'coastal']

    for m in range(len(modes_cols)):
        mode_data_path = os.path.join(output_path, 'flow_mapping_paths',
                                      'national_scale_flow_paths_100_percent.xlsx')

        flow = pd.read_excel(mode_data_path,sheet_name=modes_cols[m])
        flow = flow[['min_edge_path','max_edge_path']]
        diff = 0
        for iter_,row in flow.iterrows():
            if row[0] != row[1]:
                diff += 1

        print ('Percentage of changing paths in {} OD flows {}'.format(modes_cols[m],100.0*diff/len(flow.index)))

    provinces = ['Lao Cai','Binh Dinh','Thanh Hoa']

    for province in provinces:
        province_name = province.replace(' ', '').lower()

        mode_data_path = os.path.join(output_path, 'flow_mapping_paths',
                                      'province_roads_commune_center_access_flow_paths_100_percent.xlsx')

        flow = pd.read_excel(mode_data_path,sheet_name='{}_5_tons'.format(province_name))
        flow = flow[['min_edge_path','max_edge_path']]
        diff = 0
        for iter_,row in flow.iterrows():
            if row[0] != row[1]:
                diff += 1

        print ('Percentage of changing paths in {} OD flows {}'.format(province,100.0*diff/len(flow.index)))




        


if __name__ == '__main__':
    main()
