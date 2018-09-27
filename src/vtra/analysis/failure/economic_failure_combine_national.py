# -*- coding: utf-8 -*-
"""
Python script to assign commodity flows on the road network
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
	if abs(x[min_col]) > abs(x[max_col]):
		return x[max_col],x[min_col]
	else:
		return x[min_col],x[max_col]

def network_failure_assembly(edge_failure_dataframe,transport_mode,gdf_edges,save_edges = True,output_path =''):
	"""
	Assign net revenue to roads assets in Vietnam

	Inputs are:
	start_points - GeoDataFrame of start points for shortest path analysis.
	end_points - GeoDataFrame of potential end points for shorest path analysis.
	G - iGraph network of the province.
	save_edges -

	Outputs are:
	Shapefile with all edges and the total net reveneu transferred along each edge
	GeoDataFrame of total net revenue transferred along each edge
	"""
	failure_columns = edge_failure_dataframe.columns.values.tolist()
	failure_columns = [f for f in failure_columns if f != 'edge_id']


	drop_cols = [c for c in gdf_edges.columns.values.tolist() if ('econ' in c) or (c in ('min_loss','max_loss'))]
	if drop_cols:
		for d in drop_cols:
			gdf_edges.drop(d,axis=1,inplace=True)

	for fc in failure_columns:
		gdf_edges[fc] = 0

	for iter_, row in edge_failure_dataframe.iterrows():
		# print (row[1:])
		gdf_edges.loc[gdf_edges['edge_id'] == row['edge_id'],failure_columns] = row[failure_columns].values

	# gdf_edges[min_ind_cols] = gdf_edges['min_vals'].apply(pd.Series)
	# gdf_edges[max_ind_cols] = gdf_edges['max_vals'].apply(pd.Series)
	# gdf_edges.drop('min_vals',axis=1,inplace=True)
	# gdf_edges.drop('max_vals',axis=1,inplace=True)
	gdf_edges['min_loss'] = gdf_edges['min_tr_los'] + gdf_edges['min_econ_loss']
	gdf_edges['max_loss'] = gdf_edges['max_tr_los'] + gdf_edges['max_econ_loss']

	failure_columns += ['min_loss','max_loss']
	industry_columns = list(set([f.split('min_')[1] for f in failure_columns if 'min' in f]))

	for ind in industry_columns:
		gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)),axis = 1)
		gdf_edges[['min_{}'.format(ind),'max_{}'.format(ind)]] = gdf_edges['swap'].apply(pd.Series)
		gdf_edges.drop('swap',axis=1,inplace=True)

	if save_edges == True:
		gdf_edges.to_file(os.path.join(output_path,'weighted_edges_failures_national_{0}_multi_modal_options.shp'.format(transport_mode)))
		# gdf_edges.to_file(os.path.join(output_path,'weighted_edges_failures_national_{0}_2.shp'.format(transport_mode)))

	del gdf_edges, edge_failure_dataframe

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	types = ['min','max']

	shp_output_path = os.path.join(output_path,'failure_shapefiles')
	econ_paths_data = os.path.join(data_path,'Results','Economic_Failure_Results','summarized')

	'''
	Get the modal shares
	'''
	# modes = ['road','rail','inland','coastal']
	modes = ['road','rail']

	for m in range(len(modes)):
		mode_data_file = os.path.join(shp_output_path,'weighted_edges_failures_national_{0}_multi_modal_options.shp'.format(modes[m]))
		# mode_data_file = os.path.join(shp_output_path,'weighted_edges_failures_national_{0}_2.shp'.format(modes[m]))
		G_df = gpd.read_file(mode_data_file)

		print (G_df)
		edge_fail_ranges = []
		for t in range(len(types)):
			df_path = os.path.join(econ_paths_data,'single_edge_failures_totals_national_{0}_{1}_multi_modal_options_summarized.csv'.format(modes[m],types[t]))
			# df_path = os.path.join(econ_paths_data,'single_edge_failures_totals_national_{0}_{1}_summarized.csv'.format(modes[m],types[t]))
			df = pd.read_csv(df_path,index_col = 0)
			df.index.names = ['edge_id']
			df = df.reset_index()

			df['total_losses'] = -1e6*df['total_losses']
			df = df[df['total_losses'] > 0]
			df.rename(columns={'total_losses': '{}_econ_loss'.format(types[t])}, inplace=True)

			edge_fail_ranges.append(df)
			del df

		edge_impact = edge_fail_ranges[0]
		edge_impact = pd.merge(edge_impact,edge_fail_ranges[1],how='left', on=['edge_id']).fillna(0)
		print (edge_impact)

		network_failure_assembly(edge_impact,modes[m],G_df,save_edges = True,output_path =shp_output_path)




if __name__ == "__main__":
	main()
