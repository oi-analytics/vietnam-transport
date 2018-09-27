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
import matplotlib.pyplot as plt
import matplotlib as mpl


from vtra.utils import load_config
from vtra.transport_network_creation import province_shapefile_to_network, add_igraph_generalised_costs_province_roads, province_shapefile_to_dataframe


def main():
	mpl.style.use('ggplot')
	data_path,calc_path,output_path = load_config()['paths']['data'],load_config()['paths']['calc'],load_config()['paths']['output']

	truck_unit_wt = [5.0,20.0]
	duration_max = [10,15,20,25,30]
	# provinces to consider
	province_list = ['Lao Cai','Binh Dinh','Thanh Hoa']
	province_terrian = ['mountain','flat','flat']

	growth_scenarios = [(5,'low'),(6.5,'forecast'),(10,'high')]
	base_year = 2016

	index_cols = ['edge_id','hazard_type','model','climate_scenario','year','road_cond','asset_type','width','road_length']
	new_cols = ['min_band','max_band','min_height','max_height','min_exposure_percent','max_exposure_percent',
				'min_duration','max_duration','min_exposure_length','max_exposure_length','min_daily_loss','max_daily_loss',
				'min_npv_nooption','max_npv_nooption','adapt_strategy','min_npv','max_npv','min_bc_ratio','max_bc_ratio']

	index_cols = index_cols + new_cols
	'''
	Path OD flow disruptions
	'''
	# for prn in range(len(province_list)):
	for prn in range(0,1):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		province_name = province.replace(' ','').lower()
		for tr_wt in truck_unit_wt:
			adapt_output_excel = os.path.join(output_path,'failure_results','single_edge_failures_scenarios_{0}_{1}_tons_adapt_options.xlsx'.format(province_name,int(tr_wt)))
			for grth in growth_scenarios:
				for dur in duration_max:
					scenarios_df = pd.read_excel(adapt_output_excel,sheet_name = grth[1] + '_' + str(dur))
					scenarios_df = scenarios_df[scenarios_df['max_daily_loss'] > 0]
					scenarios_df['min_npv'] = 1.0e-6*scenarios_df['min_npv']
					scenarios_df['max_npv'] = 1.0e-6*scenarios_df['max_npv']

					# scenarios_df = scenarios_df.set_index(['model','climate_scenario','year','adapt_strategy'])
					# scenarios = list(set(scenarios_df.index.values.tolist()))
					# scenarios_df = scenarios_df.reset_index()

					adapt_strategy = list(set(scenarios_df['adapt_strategy'].values.tolist()))
					for adapt in adapt_strategy:
						plt.figure(figsize=(8,4))
						adapt_df = scenarios_df[scenarios_df['adapt_strategy'] == adapt]
						adapt_df = adapt_df.set_index(['model','climate_scenario','year','hazard_type'])
						scenarios = list(set(adapt_df.index.values.tolist()))
						adapt_df = adapt_df.reset_index()
						adapt_df['edge_nos'] = adapt_df.index.values.tolist()

						for sc in scenarios:
							model = sc[0]
							cl = sc[1]
							year = sc[2]
							if sc[1] == 'none':
								label = '{0} {1} {2}'.format(sc[3],sc[0],sc[2])
							else:
								label = '{0} {1} {2} {3}'.format(sc[3],sc[0],sc[1],sc[2])

							sc_df = adapt_df[(adapt_df['model'] == sc[0]) & (adapt_df['climate_scenario'] == sc[1]) & (adapt_df['year'] == sc[2]) & (adapt_df['hazard_type'] == sc[3])]
							# sc_df.plot('edge_nos','min_npv',linestyle = '-',Linewidth = 1,label = label)
							# sc_df.plot('edge_nos','max_npv',linestyle = '-',Linewidth = 1,label = label)
							# sc_df = sc_df.sort_values(by = ['max_npv'])
							# xs = sc_df.edge_nos
							# ys = sc_df.min_npv
							# plt.plot(xs, ys, label=label)
							ys = np.sort(sc_df.max_npv)
							xs = np.arange(0,len(ys),1)
							plt.plot(xs ,ys, label=label)

						plt.xlabel('Number of road assests', fontweight='bold')
						plt.ylabel('Net Present Value (USD million)', fontweight='bold')
						plt.title('Adaptation option: {}'.format(adapt.title()))
						plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)

						fig_path = os.path.join(output_path,'failure_results','Figures','{0}_{1}_tons_adapt_options_{2}_growth_{3}_days_{4}.png'.format(province_name,int(tr_wt),grth[1],dur,adapt))
						plt.tight_layout()
						plt.savefig(fig_path,dpi=500)
						plt.close()












if __name__ == "__main__":
	main()
