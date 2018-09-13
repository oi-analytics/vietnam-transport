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
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import *
from scripts.transport_network_creation import *

def swap_min_max(x,min_col,max_col):
	'''
	'''
	if abs(x[min_col]) > abs(x[max_col]):
		return x[max_col],x[min_col]
	else:
		return x[min_col],x[max_col]


def network_failure_assembly(edge_failure_dataframe,vehicle_wt,transport_mode,gdf_edges,save_edges = True,output_path =''):
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
	# min_ind_cols = []
	# max_ind_cols = []
	# ch_min_ind_cols = []
	# ch_max_ind_cols = []
	# for ind in industry_columns:
	# 	min_ind_cols.append('min_{}'.format(ind))
	# 	max_ind_cols.append('max_{}'.format(ind))
	# 	if ind in ('rice','tons'):
	# 		ch_min_ind_cols.append('min_{}'.format(ind))
	# 		ch_max_ind_cols.append('max_{}'.format(ind))
	# 	else:
	# 		ch_min_ind_cols.append(ind)
	# 		ch_max_ind_cols.append(ind)

	# print (len(ch_min_ind_cols))
	# print (len(ch_max_ind_cols))
	# print (len(min_ind_cols))
	# print (len(max_ind_cols))

	failure_columns = edge_failure_dataframe.columns.values.tolist()
	failure_columns = [f for f in failure_columns if f != 'edge_id']

	for fc in failure_columns:
		gdf_edges[fc] = 0

	for iter_, row in edge_failure_dataframe.iterrows():
		# print (row[1:])
		gdf_edges.loc[gdf_edges['edge_id'] == row['edge_id'],failure_columns] = row[failure_columns].values
	
	# gdf_edges[min_ind_cols] = gdf_edges['min_vals'].apply(pd.Series)
	# gdf_edges[max_ind_cols] = gdf_edges['max_vals'].apply(pd.Series)
	# gdf_edges.drop('min_vals',axis=1,inplace=True)
	# gdf_edges.drop('max_vals',axis=1,inplace=True)

	industry_columns = list(set([f.split('min_')[1] for f in failure_columns if 'min' in f]))

	for ind in industry_columns:
		gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)),axis = 1)
		gdf_edges[['min_{}'.format(ind),'max_{}'.format(ind)]] = gdf_edges['swap'].apply(pd.Series)
		gdf_edges.drop('swap',axis=1,inplace=True)

	if save_edges == True:
		gdf_edges.to_file(os.path.join(output_path,'weighted_edges_failures_national_{0}_10_percent_shift.shp'.format(transport_mode)))

	del gdf_edges, edge_failure_dataframe

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	types = ['min','max']
	path_types = ['min_edge_path','max_edge_path']
	tons_types = ['min_tons','max_tons']
	dist_types = ['min_distance','max_distance']
	time_types = ['min_time','max_time']
	cost_types = ['min_gcost','max_gcost']
	vehicle_types = ['min_vehicle_nums','max_vehicle_nums']
	rice_type = ['min_rice','max_rice']
	ind_crop_cols =['sugar','wood','steel','constructi','cement','fertilizer','coal','petroluem','manufactur',
				'fishery','meat','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp'] 

	shp_output_path = os.path.join(output_path,'failure_shapefiles')
	flow_paths_data = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_paths_2.xlsx')
	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','national_scale_hazard_intersections.xlsx')

	cols = ['origin','destination','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
			'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']

	
	'''
	Get the modal shares
	'''
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail'),('Waterways','waterways'),('Waterways','waterways')]
	# modes_file_paths = [('Roads','national_roads'),('Railways','national_rail')]
	modes_file_paths = [('Roads','national_roads')]
	modes = ['road','rail','inland','coastal']
	veh_wt = [20,800,800,1200]
	usage_factors = [(0,0),(0,0),(0.2,0.25),(0.2,0.25)]
	speeds = [(0,0),(40,60),(9,20),(9,20)]

	md_prop_file = os.path.join(data_path,'mode_properties','mode_costs.xlsx')
	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	for m in range(len(modes_file_paths)):
		mode_data_path = os.path.join(data_path,modes_file_paths[m][0],modes_file_paths[m][1])
		for file in os.listdir(mode_data_path):
			try:
				if file.endswith(".shp") and 'edges' in file.lower().strip():
					edges_in = os.path.join(mode_data_path, file)
			except:
				return ('Network nodes and edge files necessary')
			
		if modes[m] == 'road': 
			G_df =  national_road_shapefile_to_dataframe(edges_in,rd_prop_file)
		else:
			G_df =  network_shapefile_to_dataframe(edges_in,md_prop_file,modes[m],speeds[m][0],speeds[m][1])

		
		perct = 90
		edge_fail_ranges = []
		for t in range(len(types)):
			# df.to_csv(os.path.join(output_path,'failure_results','single_edge_failures_all_paths_national_{0}_{1}_2.csv'.format(modes[m],types[t])),index = False)
			flow_paths_data = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_paths_roads_{}_percent.xlsx'.format(perct))
			flow_df = pd.read_excel(flow_paths_data,sheet_name = modes[m])
			select_cols = ['origin','destination',tons_types[t],cost_types[t]]
			flow_df_select = flow_df[select_cols]
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_national_{0}_{1}_{2}_shift.csv'.format(modes[m],types[t],100-perct))
			df = pd.read_csv(df_path).fillna(0)
			df = df[['origin','destination','edge_id','no_access','new_cost']]
			flow_df_select = pd.merge(flow_df_select,df,on = ['origin','destination'],how = 'left').fillna(0)
			flow_df_select = flow_df_select[(flow_df_select[tons_types[t]] > 0) & (flow_df_select['edge_id'] != 0)]
			
			flow_df_select['tr_loss'] = (1 - flow_df_select['no_access'])*(flow_df_select['new_cost'] - flow_df_select[cost_types[t]])
			select_cols = ['edge_id','no_access','tr_loss']
			edge_impact = flow_df_select[select_cols]
			edge_impact_1 = edge_impact.groupby(['edge_id','no_access'])[select_cols[2:]].sum().reset_index()
			edge_impact_1.rename(columns={'tr_loss': 'tr_loss_1'}, inplace=True)

			flow_paths_data = os.path.join(output_path,'flow_mapping_paths','national_scale_flow_paths_roads_{}_percent.xlsx'.format(100-perct))
			flow_df = pd.read_excel(flow_paths_data,sheet_name = modes[m])
			select_cols = ['origin','destination',tons_types[t],cost_types[t]]
			flow_df_select = flow_df[select_cols]
			df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_national_{0}_{1}_multi_modal_options_{2}_shift.csv'.format(modes[m],types[t],100-perct))
			df = pd.read_csv(df_path).fillna(0)
			df = df[['origin','destination','edge_id','no_access','new_cost']]
			flow_df_select = pd.merge(flow_df_select,df,on = ['origin','destination'],how = 'left').fillna(0)
			flow_df_select = flow_df_select[(flow_df_select[tons_types[t]] > 0) & (flow_df_select['edge_id'] != 0)]
			
			flow_df_select['tr_loss'] = (1 - flow_df_select['no_access'])*(flow_df_select['new_cost'] - flow_df_select[cost_types[t]])
			select_cols = ['edge_id','no_access','tr_loss']
			edge_impact = flow_df_select[select_cols]
			edge_impact_2 = edge_impact.groupby(['edge_id','no_access'])[select_cols[2:]].sum().reset_index()
			edge_impact_2.rename(columns={'tr_loss': 'tr_loss_2'}, inplace=True)

			edge_impact = pd.merge(edge_impact_1,edge_impact_2,on = ['edge_id','no_access'],how = 'left').fillna(0)
			edge_impact['tr_loss'] = edge_impact['tr_loss_1'] + edge_impact['tr_loss_2']
			edge_impact.drop(['tr_loss_1','tr_loss_2'],axis=1,inplace=True)
			edge_impact.rename(columns={'tr_loss': '{}_tr_loss'.format(types[t])}, inplace=True)

			
			# df_path = os.path.join(output_path,'failure_results','single_edge_failures_all_path_impacts_national_{0}_{1}_{2}_shift.csv'.format(modes[m],types[t],100-perct))
			# flow_df_select.to_csv(df_path,index = False)
			# flow_df_select = pd.read_csv(df_path).fillna(0)
			# flow_df_select.rename(columns={'transport_loss': 'tr_loss'}, inplace=True)

			# select_cols = ['edge_id','o_region','d_region','no_access'] + ind_crop_cols + [rice_type[t],tons_types[t]]
			# edge_impact = flow_df_select[select_cols]
			# edge_impact = edge_impact[edge_impact['no_access'] == 1]
			# edge_impact = edge_impact.groupby(['edge_id', 'o_region','d_region'])[ind_crop_cols + [rice_type[t],tons_types[t]]].sum().reset_index()
			# df_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_national_{0}_{1}_2.csv'.format(modes[m],types[t]))
			# edge_impact.to_csv(df_path,index = False)
			# edge_fail_ranges.append(edge_impact)
			# edge_impact = flow_df_select[select_cols+['tr_loss']]
			# edge_impact = edge_impact[edge_impact['no_access'] == 0]
			# edge_impact = edge_impact.groupby(['edge_id', 'o_region','d_region'])['tr_loss'].sum().reset_index()
			# df_path = os.path.join(output_path,'failure_results','single_edge_failures_tr_loss_national_{0}_{1}_2.csv'.format(modes[m],types[t]))
			# edge_impact.to_csv(df_path,index = False)

			# select_cols = ['edge_id','no_access','tr_loss'] + ind_crop_cols + [rice_type[t],tons_types[t]]
			# select_cols = ['edge_id','no_access','tr_loss'] + [tons_types[t]]
			# edge_impact = flow_df_select[select_cols]
			# edge_impact = edge_impact.groupby(['edge_id','no_access'])[select_cols[2:]].sum().reset_index()
			# for col_name in [c for c in select_cols if c not in ['edge_id','no_access',rice_type[t],tons_types[t]]]:
			# 	edge_impact.rename(columns={col_name: '{}_'.format(types[t])+col_name}, inplace=True)

			edge_fail_ranges.append(edge_impact)
			del edge_impact

		edge_impact = edge_fail_ranges[0]
		edge_impact = pd.merge(edge_impact,edge_fail_ranges[1],how='left', on=['edge_id','no_access']).fillna(0)
		df_path = os.path.join(output_path,'failure_results','single_edge_failures_transport_losses_national_{0}_{1}_shift.csv'.format(modes[m],100 - perct))
		edge_impact.to_csv(df_path,index = False)

		network_failure_assembly(edge_impact,veh_wt[m],modes[m],G_df,save_edges = True,output_path =shp_output_path)


			

if __name__ == "__main__":
	main()