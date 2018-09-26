# -*- coding: utf-8 -*-
"""
Python script to assign commodity flows on the road network in Tanzania
Created on Wed March 06 2018

@author: Raghav Pant
"""

import pandas as pd
import os
import psycopg2
import networkx as nx
import csv
import itertools
import operator
import ast
from sqlalchemy import create_engine
import numpy as np
import igraph as ig
import copy
from collections import Counter
import sys
import math
import copy
import geopandas as gpd

from vtra.utils import *
from vtra.transport_network_creation import *

def swap_min_max(x,min_col,max_col):
	'''
	'''
	if x[min_col] > x[max_col]:
		return x[max_col],x[min_col]
	else:
		return x[min_col],x[max_col]

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	# provinces to consider

	growth_scenarios = [(5,'low'),(6.5,'forecast'),(10,'high')]
	base_year = 2016
	types = ['min','max']
	# modes = ['road','rail','inland','coastal']
	modes = ['road','rail']
	for m in range(len(modes)):
		flow_output_excel = os.path.join(output_path,'failure_results','single_edge_failures_totals_national_{0}_projections.xlsx'.format(modes[m]))
		excl_wrtr_1 = pd.ExcelWriter(flow_output_excel)
		for grth in growth_scenarios:
			edge_fail_ranges = []
			for t in range(len(types)):
				df_path = os.path.join(output_path,'failure_shapefiles','weighted_edges_failures_national_{}.shp'.format(modes[m]))
				df = gpd.read_file(df_path)
				df = df[['edge_id','{}_loss'.format(types[t])]]
				edge_impact = df[['edge_id']]
				# print (edge_impact)
				cols = []
				for year in range(2016,2050):
					edge_impact['{0}_loss_{1}'.format(types[t],year)] = math.pow((1+grth[0]/100),year - base_year)*df['{}_loss'.format(types[t])]

					cols += ['{0}_loss_{1}'.format(types[t],year)]

				edge_impact = edge_impact.groupby(['edge_id'])[cols].sum().reset_index()
				edge_fail_ranges.append(edge_impact)

				print ('Done with national {0} {1} with growth rate {2}'.format(modes[m],types[t],grth))

			edge_impact = edge_fail_ranges[0]
			edge_impact = pd.merge(edge_impact,edge_fail_ranges[1],how='left', on=['edge_id']).fillna(0)

			edge_impact.to_excel(excl_wrtr_1,grth[1],index = False)
			excl_wrtr_1.save()



if __name__ == "__main__":
	main()
