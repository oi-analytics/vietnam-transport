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

	# truck_unit_wt = [5.0,20.0]
	truck_unit_wt = [5.0]
	# provinces to consider
	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']
	province_terrian = ['mountain','flat','flat']

	growth_scenarios = [(5,'low'),(6.5,'forecast'),(10,'high')]
	base_year = 2016
	types = ['min','max']
	# path_types = ['min_edge_path','max_edge_path']
	# rev_types = ['min_netrev','max_netrev']
	# tons_types = ['min_croptons','max_croptons']
	# dist_types = ['min_distance','max_distance']
	# time_types = ['min_time','max_time']
	# cost_types = ['min_gcost','max_gcost']
	# vechicle_types = ['min_vehicle_nums','max_vehicle_nums']

	# flow_paths_data = os.path.join(output_path,'flow_mapping_paths','province_roads_district_center_flow_paths.xlsx')
	# flow_paths_data = os.path.join(output_path,'flow_mapping_paths','province_roads_commune_center_flow_paths.xlsx')
	# fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')

	# rd_prop_file = os.path.join(data_path,'Roads','road_properties','road_properties.xlsx')

	# cols = ['origin','destination','min_edge_path','max_edge_path','min_netrev','max_netrev','min_croptons','max_croptons',
	# 		'min_distance','max_distance','min_time','max_time','min_gcost','max_gcost']

	'''
	Path OD flow disruptions
	'''
	for prn in range(len(province_list)):
	# for prn in range(0,1):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()
		for tr_wt in truck_unit_wt:
			flow_output_excel = os.path.join(output_path,'failure_results','single_edge_failures_commune_access_totals_{0}_{1}_tons_projections.xlsx'.format(province_name,int(tr_wt)))
			excl_wrtr_1 = pd.ExcelWriter(flow_output_excel)
			for grth in growth_scenarios:
				edge_fail_ranges = []
				for t in range(len(types)):
					df_path = os.path.join(output_path,'failure_results','single_edge_failures_commune_access_all_path_impacts_{0}_{1}_{2}_tons.csv'.format(province_name,types[t],int(tr_wt)))
					df = pd.read_csv(df_path).fillna(0)
					edge_impact = df[['edge_id']]
					# print (edge_impact)
					cols = []
					for year in range(2016,2050):
						edge_impact['{0}_econ_value_{1}'.format(types[t],year)] = math.pow((1+grth[0]/100),year - base_year)*df['econ_value']
						edge_impact['{0}_tons_{1}'.format(types[t],year)] = math.pow((1+grth[0]/100),year - base_year)*df['tons']
						edge_impact['{0}_econ_loss_{1}'.format(types[t],year)] = df['no_access']*edge_impact['{0}_econ_value_{1}'.format(types[t],year)] + (1 - df['no_access'])*np.maximum(1,np.ceil(edge_impact['{0}_tons_{1}'.format(types[t],year)]/tr_wt))*(df['new_cost'] - df['old_cost'])

						cols += ['{0}_econ_value_{1}'.format(types[t],year),'{0}_tons_{1}'.format(types[t],year),'{0}_econ_loss_{1}'.format(types[t],year)]

					edge_impact = edge_impact.groupby(['edge_id'])[cols].sum().reset_index()
					edge_fail_ranges.append(edge_impact)

					print ('Done with province {0} {1} tons {2} with growth rate {3}'.format(province_name,tr_wt,types[t],grth))

				edge_impact = edge_fail_ranges[0]
				edge_impact = pd.merge(edge_impact,edge_fail_ranges[1],how='left', on=['edge_id']).fillna(0)
				for year in range(2016,2050):
					edge_impact['swap'] = edge_impact.apply(lambda x: swap_min_max(x,'min_econ_value_{0}'.format(year),'max_econ_value_{0}'.format(year)),axis = 1)
					edge_impact[['min_econ_value_{0}'.format(year),'max_econ_value_{0}'.format(year)]] = edge_impact['swap'].apply(pd.Series)
					edge_impact.drop('swap',axis=1,inplace=True)


					edge_impact['swap'] = edge_impact.apply(lambda x: swap_min_max(x,'min_tons_{0}'.format(year),'max_tons_{0}'.format(year)),axis = 1)
					edge_impact[['min_tons_{0}'.format(year),'max_tons_{0}'.format(year)]] = edge_impact['swap'].apply(pd.Series)
					edge_impact.drop('swap',axis=1,inplace=True)

					edge_impact['swap'] =  edge_impact.apply(lambda x: swap_min_max(x,'min_econ_loss_{0}'.format(year),'max_econ_loss_{0}'.format(year)),axis = 1)
					edge_impact[['min_econ_loss_{0}'.format(year),'max_econ_loss_{0}'.format(year)]] = edge_impact['swap'].apply(pd.Series)
					edge_impact.drop('swap',axis=1,inplace=True)

				edge_impact.to_excel(excl_wrtr_1,grth[1],index = False)
				excl_wrtr_1.save()



if __name__ == "__main__":
	main()
