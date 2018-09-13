# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 12:00:00 2018
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


def net_present_value(probability_wt,duration,adaptation_options_dataframe,strategy_parameter,parameter_value_list,
				min_economic_loss,max_economic_loss,edge_options_dictionary,start_year,end_year,growth_rate,discount_rate,edge_width = 1.0, edge_length = 1.0):
	'''
	strategy	
	asset_type	
	asset_cond	
	location	
	height_m	
	adapt_cost_min	
	adapt_cost_max	
	maintain_cost_min	
	maintain_cost_max	
	rehab_cost_min	
	rehab_cost_max	
	height_incr_min	
	height_incr_max	
	maintenance_times_min	
	maintenance_times_max	
	cost_unit
	'''
	for param_val in parameter_value_list:
		
		st = adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'strategy'].values[0]
		
		cost_unit = adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'cost_unit'].values[0]
		if cost_unit == 'million USD/km':
			cost_fact = 1000000.0
		else:
			cost_fact = 1.0

		adapt_cost_min = cost_fact*edge_width*edge_length*adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'min_initial_cost'].values[0]
		adapt_cost_max = cost_fact*edge_width*edge_length*adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'max_initial_cost'].values[0]
		

		st_min_rehab_benefit = cost_fact*probability_wt*edge_width*edge_length*adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'min_rehab_benefit'].sum()
		st_max_rehab_benefit = cost_fact*probability_wt*edge_width*edge_length*adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'max_rehab_benefit'].sum()

		st_min_cost = cost_fact*edge_width*edge_length*adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'min_cost'].sum()
		st_max_cost = cost_fact*edge_width*edge_length*adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter] == param_val,'max_cost'].sum()
		
		total_discount_ratio = []
		for year in range(start_year,end_year):
			# total_discount_ratio += 1.0/math.pow(1.0 + 1.0*discount_rate/100.0,year - start_year)
			total_discount_ratio.append(1.0*math.pow(1.0 + 1.0*growth_rate/100.0,year - start_year)/math.pow(1.0 + 1.0*discount_rate/100.0,year - start_year))

		min_eal = probability_wt*duration*sum(min_economic_loss*np.array(total_discount_ratio))
		max_eal = probability_wt*duration*sum(max_economic_loss*np.array(total_discount_ratio))
		
		st_min_benefit = st_min_rehab_benefit + min_eal
		st_max_benefit = st_max_rehab_benefit + max_eal

		min_npv = st_min_benefit - st_max_cost 
		max_npv = st_max_benefit - st_min_cost

		if st_max_cost > 0:
			min_bc_ratio = 1.0*st_min_benefit/st_max_cost
		else:
			min_bc_ratio = 0

		if st_min_cost > 0:
			max_bc_ratio = 1.0*st_max_benefit/st_min_cost
		else:
			max_bc_ratio = 0

		edge_options_dictionary.append({'adapt_strategy':st,'min_initial_cost':adapt_cost_min,'max_initial_cost':adapt_cost_max,
								'min_benefit_npv':st_min_benefit,'max_benefit_npv':st_max_benefit,
								'min_cost_npv':st_min_cost, 'max_cost_npv':st_max_cost,
								'min_adapt_npv':min_npv,'max_adapt_npv':max_npv,'min_bc_ratio':min_bc_ratio,'max_bc_ratio':max_bc_ratio})

	return edge_options_dictionary

def main():
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	'''
	cols = ['band_name','band_num','climate_scenario','commune_id','commune_name','district_id','district_name',
			'edge_id','hazard_type','max_val','min_val','model','probability','province_id','province_name','sector',
			'year','length','road_cond','asset_type','width','min_econ_loss','max_econ_loss']
	'''

	start_year = 2016
	end_year = 2050

	discount_rate = 12.0
	total_discount_ratio = []
	for year in range(start_year,end_year):
		# total_discount_ratio += 1.0/math.pow(1.0 + 1.0*discount_rate/100.0,year - start_year)
		total_discount_ratio.append(1.0/math.pow(1.0 + 1.0*discount_rate/100.0,year - start_year))

	# total_discount_ratio = sum(total_discount_ratio_list)

	df_path =  os.path.join(data_path,'Adaptation_options','adaptation_options.xlsx')
	adaptation_df = pd.read_excel(df_path,sheet_name = 'adapt_options')

	adaptation_df['total_discount_ratio'] = sum(total_discount_ratio)

	min_maintain_discount_ratio_list = []
	max_maintain_discount_ratio_list = []

	for iter_, row in adaptation_df.iterrows():
		min_maintain_schedule = row['maintenance_times_min']
		max_maintain_schedule = row['maintenance_times_max']

		min_maintain_discount_ratio = 0
		max_maintain_discount_ratio = 0

		max_maintain_discount_years = np.arange(start_year,end_year,min_maintain_schedule)
		min_maintain_discount_years = np.arange(start_year,end_year,max_maintain_schedule)
		for year in max_maintain_discount_years[1:]:
			max_maintain_discount_ratio += 1.0/math.pow(1.0 + 1.0*discount_rate/100.0,year - start_year)

		for year in min_maintain_discount_years[1:]:
			min_maintain_discount_ratio += 1.0/math.pow(1.0 + 1.0*discount_rate/100.0,year - start_year)

		min_maintain_discount_ratio_list.append(min_maintain_discount_ratio)
		max_maintain_discount_ratio_list.append(max_maintain_discount_ratio)


	adaptation_df['min_maintain_discount_ratio'] = min_maintain_discount_ratio_list
	adaptation_df['max_maintain_discount_ratio'] = max_maintain_discount_ratio_list

	adaptation_df['min_rehab_benefit'] = adaptation_df['rehab_cost_min']*adaptation_df['total_discount_ratio']
	adaptation_df['max_rehab_benefit'] = adaptation_df['rehab_cost_max']*adaptation_df['total_discount_ratio']

	adaptation_df['min_cost'] = adaptation_df['adapt_cost_min']+ adaptation_df['maintain_cost_min']*adaptation_df['min_maintain_discount_ratio']
	adaptation_df['max_cost'] = adaptation_df['adapt_cost_max'] + adaptation_df['maintain_cost_max']*adaptation_df['max_maintain_discount_ratio']
	adaptation_df['min_initial_cost'] = adaptation_df['adapt_cost_min']
	adaptation_df['max_initial_cost'] = adaptation_df['adapt_cost_max']
	
	print (adaptation_df)

	cols = ['strategy','asset_type','asset_cond','location','height_m','cost_unit','min_initial_cost',
			'max_initial_cost','min_cost','max_cost','min_rehab_benefit','max_rehab_benefit']

	ad_opt = []
	st_desc = ['road change','roads','unpaved-paved','all',0,'$/m2']
	st_min_ini_cost = adaptation_df.loc[adaptation_df['asset_cond'] == 'paved','min_initial_cost'].sum()
	st_max_ini_cost = adaptation_df.loc[adaptation_df['asset_cond'] == 'paved','max_initial_cost'].sum()

	st_min_cost = adaptation_df.loc[adaptation_df['asset_cond'] == 'paved','min_cost'].sum()
	st_max_cost = adaptation_df.loc[adaptation_df['asset_cond'] == 'paved','max_cost'].sum()

	st_min_benefit = adaptation_df.loc[adaptation_df['asset_cond'] == 'unpaved','rehab_cost_min'].sum()*total_discount_ratio[0] + adaptation_df.loc[adaptation_df['asset_cond'] == 'paved','rehab_cost_min'].sum()*sum(total_discount_ratio[1:])
	st_max_benefit = adaptation_df.loc[adaptation_df['asset_cond'] == 'unpaved','rehab_cost_max'].sum()*total_discount_ratio[0] + adaptation_df.loc[adaptation_df['asset_cond'] == 'paved','rehab_cost_max'].sum()*sum(total_discount_ratio[1:])
	
	ad_opt.append(st_desc + [st_min_ini_cost,st_max_ini_cost,st_min_cost,st_max_cost,st_min_benefit,st_max_benefit])
	new_ht = 6
	st_desc = ['dyke building','dyke','rural','sea',new_ht,'million USD/km']
	st_min_ini_cost = adaptation_df.loc[adaptation_df['height_m'] == 4,'min_initial_cost'].sum() + (new_ht - 4)*adaptation_df.loc[adaptation_df['height_m'] == 4,'height_incr_min'].sum()
	st_max_ini_cost = adaptation_df.loc[adaptation_df['height_m'] == 4,'max_initial_cost'].sum() + (new_ht - 4)*adaptation_df.loc[adaptation_df['height_m'] == 4,'height_incr_max'].sum()


	st_min_cost = adaptation_df.loc[adaptation_df['height_m'] == 4,'min_cost'].sum() + (new_ht - 4)*adaptation_df.loc[adaptation_df['height_m'] == 4,'height_incr_min'].sum()
	st_max_cost = adaptation_df.loc[adaptation_df['height_m'] == 4,'max_cost'].sum() + (new_ht - 4)*adaptation_df.loc[adaptation_df['height_m'] == 4,'height_incr_max'].sum()

	st_min_benefit = adaptation_df.loc[adaptation_df['height_m'] == 4,'min_rehab_benefit'].sum()
	st_max_benefit = adaptation_df.loc[adaptation_df['height_m'] == 4,'max_rehab_benefit'].sum()
	ad_opt.append(st_desc + [st_min_ini_cost,st_max_ini_cost,st_min_cost,st_max_cost,st_min_benefit,st_max_benefit])

	ad_opt_df = pd.DataFrame(ad_opt,columns = cols)

	adaptation_options = adaptation_df[cols]
	adaptation_options = adaptation_options.append(ad_opt_df, ignore_index=True)

	del ad_opt_df
	print (adaptation_options)

	cols = ['band_num','climate_scenario',
			'edge_id','hazard_type','max_val','min_val','model','probability',
			'year','exposure_length']

	# index_cols = ['edge_id','hazard_type','model','climate_scenario','year','road_cond','asset_type','width','road_length']

	# provinces to consider 
	growth_scenarios = [(5,'low'),(6.5,'forecast'),(10,'high')]
	base_year = 2016
	types = ['min','max']

	fail_scenarios_data = os.path.join(output_path,'hazard_scenarios','national_scale_hazard_intersections.xlsx')
	rd_prop_file = os.path.join(data_path,'mode_properties','road_properties.xlsx')

	duration_max = [10,15,20,25,30]
	length_thr = 500.0

	modes = ['road']
	for m in range(len(modes)):
		index_cols = ['edge_id','hazard_type','model','climate_scenario','year','road_cond','width','road_length','min_daily_loss_{}'.format(start_year),'max_daily_loss_{}'.format(start_year)]
		all_edge_fail_scenarios = pd.read_excel(fail_scenarios_data,sheet_name = modes[m])
		all_edge_fail_scenarios.loc[all_edge_fail_scenarios['probability'] == 'none', 'probability'] = 1.0
		all_edge_fail_scenarios['probability'] = pd.to_numeric(all_edge_fail_scenarios['probability'])
		all_edge_fail_scenarios.rename(columns={'length': 'exposure_length'}, inplace=True)
		all_edge_fail_scenarios = all_edge_fail_scenarios[cols]
		all_edge_fail_scenarios = all_edge_fail_scenarios.drop_duplicates(subset=cols, keep='first')

		edges_in = os.path.join(data_path,'Results','Failure_shapefiles','weighted_edges_failures_national_{}_2.shp'.format(modes[m]))
		edges = gpd.read_file(edges_in)
		edges = edges[['edge_id','road_cond','width','length','min_loss','max_loss']]

		all_edge_fail_scenarios = pd.merge(all_edge_fail_scenarios,edges,on = ['edge_id'],how = 'left').fillna(0)
		del edges

		all_edge_fail_scenarios.rename(columns={'length': 'road_length','min_loss':'min_daily_loss_{}'.format(start_year),'max_loss':'max_daily_loss_{}'.format(start_year)}, inplace=True)
		all_edge_fail_scenarios = all_edge_fail_scenarios[all_edge_fail_scenarios['max_daily_loss_{}'.format(start_year)] > 0]
		all_edge_fail_scenarios['road_length'] = 1000.0*all_edge_fail_scenarios['road_length']

		# all_edge_fail_scenarios = all_edge_fail_scenarios.groupby(index_cols + ['probability'])['exposure_length'].sum().reset_index()
		all_edge_fail_scenarios['percent_exposure'] = 100.0*all_edge_fail_scenarios['exposure_length']/all_edge_fail_scenarios['road_length']
		df_path = os.path.join(output_path,'hazard_scenarios','roads_hazard_intersections_national_{}.csv'.format(modes[m]))
		all_edge_fail_scenarios.to_csv(df_path,index = False)
		
		all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(index_cols)
		scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
		print ('Number of failure scenarios',len(scenarios))
		scenarios_list = []
		for sc in scenarios:
			# edge = sc[0]
			# road_cond = sc[-4]
			# width = sc[-2]
			# print (road_cond,width)

			min_height = max(all_edge_fail_scenarios.loc[[sc],'min_val'].values.tolist())
			max_height = max(all_edge_fail_scenarios.loc[[sc],'max_val'].values.tolist())
			min_band_num = min(all_edge_fail_scenarios.loc[[sc],'band_num'].values.tolist())
			max_band_num = max(all_edge_fail_scenarios.loc[[sc],'band_num'].values.tolist())
			prob = all_edge_fail_scenarios.loc[[sc],'probability'].values
			# print (sc,prob, len(prob))
			if len(list(set(prob))) > 1:
				exposure_len = all_edge_fail_scenarios.loc[[sc],'exposure_length'].values
				per = all_edge_fail_scenarios.loc[[sc],'percent_exposure'].values

				prob_tup = list(zip(prob,exposure_len,per))
				u_pr = sorted(list(set(prob.tolist())))
				exposure_len = []
				per = []
				r_wt = []
				for pr in u_pr:
					per_exp = sum([z for (x,y,z) in prob_tup if x == pr])
					exp_len = sum([y for (x,y,z) in prob_tup if x == pr])
					if per_exp > 100.0:
						exposure_len.append(100.0*exp_len/per_exp)
						per.append(100.0)
						r_wt.append(1.0)
					else:
						exposure_len.append(exp_len)
						if exp_len < length_thr:
							per.append(per_exp)
							r_wt.append(0.01*per_exp)
						else:
							per.append(100.0)
							r_wt.append(1.0)

				max_exposure_len = max(exposure_len)
				min_exposure_len = min(exposure_len)

				min_per = min(per)
				max_per = max(per)
				min_dur = 0.01*min_per
				max_dur = 0.01*max_per
				risk_wt = 0
				for p in range(len(u_pr)-1):
					risk_wt += 0.5*(u_pr[p+1]-u_pr[p])*(r_wt[p+1]+r_wt[p])

			else:
				prob_wt = prob[0]
				min_exposure_len = sum(all_edge_fail_scenarios.loc[[sc],'exposure_length'].values.tolist())
				min_per = sum(all_edge_fail_scenarios.loc[[sc],'percent_exposure'].values.tolist())
				if min_per > 100.0:
					min_exposure_len = 100.0*min_exposure_len/min_per
					min_per = 100.0

				max_per = min_per
				max_exposure_len = min_exposure_len
				
				min_dur = 0.01*min_per
				if max_exposure_len < length_thr:
					max_dur = 0.01*max_per
					risk_wt = 0.01*max_per*prob_wt
				else:
					max_dur = 1.0
					risk_wt = prob_wt



			sc_list = list(sc) + [min_band_num,max_band_num,min_height,max_height,min_per,max_per,min_dur,max_dur,min_exposure_len,max_exposure_len]
			scenarios_list.append(list(sc) + [min_band_num,max_band_num,min_height,max_height,min_per,max_per,min_dur,max_dur,min_exposure_len,max_exposure_len,risk_wt])

		
		new_cols = ['min_band','max_band','min_height','max_height','min_exposure_percent','max_exposure_percent',
					'min_duration','max_duration','min_exposure_length','max_exposure_length','risk_wt']
		scenarios_df = pd.DataFrame(scenarios_list,columns = index_cols + new_cols)
		df_path = os.path.join(output_path,'hazard_scenarios','roads_hazard_intersections_{}_risks.csv'.format(modes[m]))
		scenarios_df.to_csv(df_path,index = False)
		del all_edge_fail_scenarios

		print ('done with scenarios creation')

		# df_path = os.path.join(output_path,'hazard_scenarios','roads_hazard_intersections_{}_risks.csv'.format(province_name))
		# scenarios_df = pd.read_csv(df_path)

		index_cols = scenarios_df.columns.values.tolist()
		adapt_output_excel = os.path.join(output_path,'failure_results','single_edge_failures_scenarios_national_{0}_adapt_options.xlsx'.format(modes[m]))
		excl_wrtr = pd.ExcelWriter(adapt_output_excel)
		for gr in range(len(growth_scenarios)):
		# for gr in range(0,1):
			grth = growth_scenarios[gr]
			for dur in range(len(duration_max)):
			# for dur in range(0,1):
				scenarios_list = []
				scenarios_df['min_duration'] = duration_max[dur]*scenarios_df['min_duration']
				scenarios_df['max_duration'] = duration_max[dur]*scenarios_df['max_duration']
				# scenarios_df['risk_wt'] = duration_max[dur]*scenarios_df['risk_wt']
				for iter_, row in scenarios_df.iterrows():
					edge_options = []
					max_height = row['max_height']
					prob_wt = row['risk_wt']
					min_loss = row['min_daily_loss_{}'.format(start_year)]
					max_loss = row['max_daily_loss_{}'.format(start_year)]
					road_cond = row['road_cond']
					width = row['width']
					max_exposure_len = row['max_exposure_length']

					if max_height > 4:
						max_height = 6
					if max_height > 0:
						edge_options = net_present_value(prob_wt,duration_max[dur],adaptation_options,'height_m',[max_height],min_loss,max_loss,edge_options,
														start_year,end_year,grth[0],discount_rate,edge_width = 1.0, edge_length = 0.001*max_exposure_len)
				
					if road_cond == 'unpaved':
						edge_options = net_present_value(prob_wt,duration_max[dur],adaptation_options,'asset_cond',['unpaved','unpaved-paved'],min_loss,max_loss,
														edge_options,start_year,end_year,grth[0],discount_rate,edge_width = width, edge_length = max_exposure_len)
							
					if road_cond == 'paved':
						edge_options = net_present_value(prob_wt,duration_max[dur],adaptation_options,'asset_cond',['paved'],min_loss,max_loss,edge_options,start_year,
														end_year,grth[0],discount_rate,edge_width = width, edge_length = max_exposure_len)

					for options in edge_options:
						# scenarios_list.append(row.values.tolist() + [min_edge_ead,max_edge_ead,options['strategy'],options['min_benefit_npv'],options['max_benefit_npv'],options['min_cost_npv'],options['max_cost_npv'],options['min_adapt_npv'],options['max_adapt_npv'],options['min_bc_ratio'],options['max_bc_ratio']])
						scenarios_list.append({**dict(row),**options})
						# print (options['strategy'])



				new_cols = ['adapt_strategy','min_initial_cost','max_initial_cost','min_benefit_npv','max_benefit_npv','min_cost_npv',
							'max_cost_npv','min_adapt_npv','max_adapt_npv','min_bc_ratio','max_bc_ratio']
				scenarios_list_df = pd.DataFrame(scenarios_list,columns = index_cols + new_cols)
				# df_path = os.path.join(output_path,'hazard_scenarios','roads_hazard_intersections_{}_risks_2.csv'.format(province_name))
				# scenarios_list_df.to_csv(df_path,index = False)
				scenarios_list_df.to_excel(excl_wrtr,grth[1] + '_' + str(duration_max[dur]),index = False)
				excl_wrtr.save()

				print ('Done with national {0} in {1} growth scenario {2} days'.format(modes[m],grth[1],duration_max[dur]))


if __name__ == "__main__":
	main()