# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 13:52:00 2017

@author: cenv0574
"""

# Import
import pandas as pd
import geopandas as gpd
import numpy as np
import sys
import os


from vtra.utils import load_config
from vtra.analysis.mria.table import io_basic
from vtra.analysis.mria.model import MRIA_IO as MRIA
from vtra.analysis.mria.disruption import create_disruption

from pathos.multiprocessing import Pool,cpu_count

def estimate_losses(input_file):

	print('{} started!'.format(input_file))

	data_path = load_config()['paths']['data']

	''' Set booleans'''
	if 'min' in input_file:
		min_rice = True
	elif 'max' in input_file:
		min_rice = False

	if 'single' in input_file:
		single_point = True
	elif 'multiple' in input_file:
		single_point = False

	''' Specify file path '''
	if min_rice == True:
		filepath =  os.path.join(data_path,'input_data','IO_VIETNAM_MIN.xlsx')
	else:
		filepath =  os.path.join(data_path,'input_data','IO_VIETNAM_MAX.xlsx')

	'''Create data input'''
	DATA = io_basic('Vietnam',filepath,2010)
	DATA.prep_data()

	'''Run model and create some output'''
	output = pd.DataFrame()

	'''Specify disruption'''
	output_dir = os.path.join(
		data_path,
		'Results',
		'Economic_Failure_Results',
		os.path.basename(os.path.splitext(input_file)[0])
	)

	'''Create output folders'''
	if os.path.exists(output_dir) == False:
		os.mkdir(output_dir)

	event_dict = create_disruption(input_file,output_dir,min_rice=min_rice,single_point=single_point)

	# print (event_dict.keys())
	collect_outputs = {}
	for iter_,event in enumerate(list(event_dict.keys())):

		if np.average(1 - np.array(list(event_dict[event].values()))) < 0.001:
			print('Event {} will cause no impacts'.format(event))
			continue

		print('Event {} started!'.format(event))

		try:
			disr_dict_fd = {}
			disr_dict_sup = event_dict[event]

			'''Get direct losses '''
			disrupt = pd.DataFrame.from_dict(disr_dict_sup,orient='index')
			disrupt.reset_index(inplace=True)
			disrupt[['region', 'sector']] = disrupt['index'].apply(pd.Series)


			'''Create model'''
			MRIA_RUN = MRIA(DATA.name,DATA.countries,DATA.sectors,EORA=False,list_fd_cats=['FinDem'])

			'''Define sets and alias'''
			# CREATE SETS
			MRIA_RUN.create_sets()

			# CREATE ALIAS
			MRIA_RUN.create_alias()

			''' Define tables and parameters'''
			MRIA_RUN.baseline_data(DATA,disr_dict_sup,disr_dict_fd)
			MRIA_RUN.impact_data(DATA,disr_dict_sup,disr_dict_fd)

			'''Get base line values'''
			output['x_in'] = pd.Series(MRIA_RUN.X.get_values())*43
			output.index.names = ['region','sector']

			'''Get direct losses '''
			disrupt = pd.DataFrame.from_dict(disr_dict_sup,orient='index')
			disrupt.reset_index(inplace=True)
			disrupt[['region', 'sector']] = disrupt['index'].apply(pd.Series)
			disrupt.drop('index',axis=1,inplace=True)
			disrupt = 1- disrupt.groupby(['region', 'sector']).sum()
			disrupt.columns = ['shock']

			output['dir_losses'] = (disrupt['shock']*output['x_in']).fillna(0)*-1

			MRIA_RUN.run_impactmodel()
			output['x_out'] = pd.Series(MRIA_RUN.X.get_values())*43
			output['total_losses'] = (output['x_out'] - output['x_in'])
			output['ind_losses'] = (output['total_losses'] - output['dir_losses'])

			output = output/365

			output = output.drop(['x_in','x_out'],axis=1)

			output.to_csv(os.path.join(output_dir,'{}.csv'.format(event)))

			prov_impact = output.groupby(level=0,axis=0).sum()
			collect_outputs[event] = prov_impact

		except Exception as e:
				print('Failed to finish {} because of {}!'.format(event,e))

	if collect_outputs:
		'''Specify disruption'''
		output_dir = os.path.join(data_path,
			'Results',
			'Economic_failure_results',
			'provincial'
		)

		'''Create output folders'''
		if os.path.exists(output_dir) == False:
			os.mkdir(output_dir)

		pd.concat(collect_outputs).to_csv(os.path.join(
			data_path,
			'Results',
			'Economic_failure_results',
			'provincial',
			'{}_provincial.csv'.format(os.path.basename(os.path.splitext(input_file)[0]))))


		get_sums = {}
		for event in collect_outputs:
			get_sums[event] = collect_outputs[event]['total_losses'].sum()

		sums = pd.DataFrame.from_dict(get_sums,orient='index')
		sums.columns = ['total_losses']

		'''Specify disruption'''
		output_dir = os.path.join(data_path,
			'Results',
			'Economic_failure_results',
			'summarized'
		)

		'''Create output folders'''
		if os.path.exists(output_dir) == False:
			os.mkdir(output_dir)

		sums.to_csv(os.path.join(data_path,
			'Results',
			'Economic_failure_results',
			'summarized',
			'{}_summarized.csv'.format(os.path.basename(os.path.splitext(input_file)[0]))))

		return pd.concat(collect_outputs),sums

if __name__ == '__main__':

	data_path = load_config()['paths']['data']

	# input_file = os.path.join(data_path,'Results','Failure_results','single_edge_failures_totals_national_road_max.csv')

	# get_all_input_files = [os.path.join(data_path,'Results','Failure_results',x) for x in os.listdir(os.path.join(data_path,'Results','Failure_results')) if x.endswith(".csv")]
	get_all_input_files = [os.path.join(data_path,'Results','Failure_results','roads',x) for x in os.listdir(os.path.join(data_path,'Results','Failure_results','roads')) if x.endswith(".csv")]

	# with Pool(int(cpu_count())-2) as pool:
	#     pool.map(estimate_losses,get_all_input_files,chunksize=1)

	for gi in get_all_input_files:
		estimate_losses(gi)
