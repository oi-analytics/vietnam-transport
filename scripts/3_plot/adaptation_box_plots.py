# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 15:42:31 2018

@author: cenv0574
"""
import os
import sys

import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches


from vtra.utils import *


mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 10.
mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 10.
mpl.rcParams['xtick.labelsize'] = 9.
mpl.rcParams['ytick.labelsize'] = 9.

def plot_boxplots(input_data,input_labels,input_colors,x_label,y_label,plot_title,plot_file_path):
	fig, ax = plt.subplots(figsize=(8,4))
	bp = ax.boxplot(input_data, meanline=True,labels=input_labels, showfliers=True,patch_artist = True)

	for element in ['boxes','whiskers', 'means', 'medians', 'caps']:
		for e, color in zip(bp[element], input_colors):
			plt.setp(e,color = '#252525',linewidth = 0.5)

	for patch, color in zip(bp['boxes'], input_colors):
		patch.set_facecolor(color)

	for flier,color in zip(bp['fliers'], input_colors):
		flier.set(marker='o', markersize=2, markerfacecolor=color,linestyle='none',markeredgecolor=color)


	ax.tick_params(axis='x', rotation=45)
	plt.xlabel(x_label, fontweight='bold')
	plt.ylabel(y_label, fontweight='bold')
	plt.title(plot_title)
	# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)

	# print (bp['boxes'][0])
	plt.tight_layout()
	plt.savefig(plot_file_path,dpi=500)
	plt.close()

def plot_boxplots_subplots(scenario_groups,input_data,x_label,y_label,plot_title,plot_file_path):
	adapt_options = list(set([x for (x,y,z) in input_data]))
	fig, axs = plt.subplots(nrows = 1, ncols = len(adapt_options),sharey = True,figsize=(8,4))

	max_input_val = np.max(np.array([np.max(z) for (x,y,z) in input_data]))
	min_input_val = np.min(np.array([np.min(z) for (x,y,z) in input_data]))
	print (min_input_val)
	if 0 <= min_input_val < 10:
		if max_input_val < 1000:
			min_input_val = -0.1*max_input_val
		else:
			min_input_val = -0.02*max_input_val
	else:
		min_input_val = 1.1*min_input_val

	if abs(max_input_val) < 0.4*abs(min_input_val):
		max_input_val = 0.4*abs(min_input_val)

	for a in range(len(axs)):
		adapt = adapt_options[a]
		gr_vals = []
		gr_colors = []
		gr_labels = []
		for gr in scenario_groups:
			vals = np.array([z for (x,y,z) in input_data if x == adapt and y == gr[0]])
			if len(vals) > 0:
				vals = np.concatenate(vals).ravel()
				# print ('values',vals)
				gr_vals.append(np.unique(vals))
				gr_labels.append('s{}'.format(gr[0]))
				gr_colors.append(gr[1])


		bp = axs[a].boxplot(gr_vals, meanline=True,labels=gr_labels, showfliers=True,patch_artist = True)

		for element in ['boxes','whiskers', 'means', 'medians', 'caps']:
			for e, color in zip(bp[element], gr_colors):
				plt.setp(e,color = '#252525',linewidth = 0.5)

		for patch, color in zip(bp['boxes'], gr_colors):
			patch.set_facecolor(color)

		for flier,color in zip(bp['fliers'], gr_colors):
			flier.set(marker='o', markersize=2, markerfacecolor=color,linestyle='none',markeredgecolor=color)


		axs[a].text(0.7*len(gr_vals),1.04*max_input_val,adapt.title(),ha="center", va="center",fontweight='bold')
		axs[a].set_ylim([min_input_val,1.2*max_input_val])
		axs[a].tick_params(axis='x', rotation=45)

	fig.text(0.5,0.04, 'Hazard scenarios', ha="center", va="center",fontweight='bold')
	fig.text(0.015,0.5, y_label, ha="center", va="center", rotation=90,fontweight='bold')

	fig.text(0.5,0.98, plot_title, ha="center", va="center",fontweight='bold')
	# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)
	fig.subplots_adjust(hspace=0)
	# print (bp['boxes'][0])
	plt.tight_layout()
	plt.savefig(plot_file_path,dpi=500)
	plt.close()

def main():
	data_path, output_path,figure_path = load_config()['paths']['data'],load_config()['paths']['output'],load_config()['paths']['figures']

	# provinces to consider
	# province_list = ['Lao Cai','Binh Dinh','Thanh Hoa','National Road']
	province_list = ['National Road']
	types = ['min','max']
	base_year = 2016
	growth_scenarios = [(5,'low'),(6.5,'forecast'),(10,'high')]
	duration_max = [10,15,20,25,30]
	index_cols = ['hazard_type','model','climate_scenario','year']
	val_cols = ['adapt_strategy','min_daily_loss_2016','max_daily_loss_2016','min_initial_cost','max_initial_cost','min_benefit_npv','max_benefit_npv','min_cost_npv','max_cost_npv','min_adapt_npv','max_adapt_npv','min_bc_ratio','max_bc_ratio']
	var_cols = ['daily_loss_2016','initial_cost','benefit_npv','cost_npv','adapt_npv','bc_ratio']

	var_dict = [
		{
			'variable':'daily_loss_2016',
			'x_label':'Hazard models',
			'y_label':'Daily economic impacts (USD million/day)',
			'title':'Hazard economic impacts',
			'divisor':1000000

		},
		{
			'variable':'initial_cost',
			'x_label':'Hazard models',
			'y_label':'Initial investment (USD million)',
			'title':'Initial cost of adaptation',
			'divisor':1000000

		},
		{
			'variable':'benefit_npv',
			'x_label':'Hazard models',
			'y_label':'Benefit NPV (USD million)',
			'title':'NPV of benefits over time',
			'divisor':1000000

		},
		{
			'variable':'cost_npv',
			'x_label':'Hazard models',
			'y_label':'Cost NPV (USD million)',
			'title':'NPV of costs over time',
			'divisor':1000000

		},
		{
			'variable':'adapt_npv',
			'x_label':'Hazard models',
			'y_label':'Adaptation NPV (USD million)',
			'title':'NPV of adaptation over time',
			'divisor':1000000

		},
		{
			'variable':'bc_ratio',
			'x_label':'Hazard models',
			'y_label':'BCR',
			'title':'BCR over time',
			'divisor':1

		},

	]

	hazard_df = pd.read_excel(os.path.join(data_path,'Hazard_data','hazard_summary.xlsx'),sheet_name = 'model_nos')
	# hazard_df = hazard_df[hazard_df['model_no'].isin(['m{}'.format(x) for x in range(16,28)])]
	hazard_df = hazard_df[hazard_df['model_no'].isin(['m1'] + ['m{}'.format(x) for x in range(3,12)])]
	scenarios = list(h for h in hazard_df.itertuples(index = False))
	scenario_groups = list(set(list(zip(hazard_df['group_no'].values.tolist(),hazard_df['color'].values.tolist()))))
	scenario_groups = [(x,y) for (x,y) in sorted(scenario_groups, key=lambda pair: pair[0])]
	print (scenario_groups)

	# for prn in range(len(province_list)):
	for prn in range(0,1):
		province = province_list[prn]
		# set all paths for all input files we are going to use
		# province_name = province.replace(' ','').lower()
		province_name = province.replace(' ','_').lower()
		# adapt_output_excel = os.path.join(output_path,'failure_results','single_edge_failures_commune_access_scenarios_{0}_5_tons_adapt_options.xlsx'.format(province_name))
		adapt_output_excel = os.path.join(output_path,'failure_results','single_edge_failures_scenarios_{0}_adapt_options.xlsx'.format(province_name))

		'''
		Box plots for the hazard losses
		'''
		# for g in range(len(growth_scenarios)):
		for g in range(1,2):
			grth = growth_scenarios[g]
			# for dur in range(len(duration_max)):
			for dur in range(0,1):
				adapt_df = pd.read_excel(adapt_output_excel,sheet_name = grth[1] + '_' + str(duration_max[dur]))
				print ('Number of options',len(adapt_df.index))
				loss_df = adapt_df[index_cols + val_cols]
				del adapt_df
				loss_df = loss_df[loss_df['max_daily_loss_2016'] > 0]
				adapt_strategy = list(set(loss_df['adapt_strategy'].values.tolist()))
				# loss_df = loss_df.set_index(index_cols)
				# scenarios = list(set(loss_df.index.values.tolist()))
				# loss_df = loss_df.reset_index()
				loss_vals = []
				adapt_dict = dict((var,[]) for var in var_cols)
				for adapt in adapt_strategy:
					sc_dict = dict((var,[]) for var in var_cols)
					label_list = []
					label_color = []
					for s in range(len(scenarios)):
						sc = scenarios[s]
						sc_df = loss_df[(loss_df.hazard_type == sc.hazard_type) & (loss_df.model == sc.model) & (loss_df.climate_scenario == sc.climate_scenario) & (loss_df.year == sc.year) & (loss_df.adapt_strategy == adapt)]
						# vals_list = []
						if len(sc_df.index) > 0:
							for var in var_dict:
								div = var['divisor']
								v = var['variable']
								vals = sc_df['min_{}'.format(v)].append(sc_df['max_{}'.format(v)])
								vals = vals.unique()
								sc_dict[v].append(1.0*vals/div)
								adapt_dict[v].append((adapt,sc.group_no,1.0*vals/div))

							label_list.append(sc.model_no)
							label_color.append(sc.color)


					# for var in var_dict:
					# 	# fig_path = os.path.join(figure_path,'{0}_5_tons_hazard_{1}_{2}.png'.format(province_name,var['variable'],adapt))
					# 	fig_path = os.path.join(figure_path,'{0}_hazard_{1}_{2}.png'.format(province_name,var['variable'],adapt))
					# 	plot_boxplots(sc_dict[var['variable']],label_list,label_color,var['x_label'],var['y_label'],province + ': ' + var['title']+ ' for {}'.format(adapt),fig_path)

				for var in var_dict:
					# fig_path = os.path.join(figure_path,'{0}_5_tons_hazard_{1}_all_options_landslise_flashflood.png'.format(province_name,var['variable']))
					fig_path = os.path.join(figure_path,'{0}_hazard_{1}_all_options_landslide_flashflood.png'.format(province_name,var['variable']))
					plot_boxplots_subplots(scenario_groups,adapt_dict[var['variable']],var['x_label'],var['y_label'],province + ': ' + var['title'],fig_path)



				# label_list = []
				# label_color = []
				# vals_list = []
				# for s in range(len(scenarios)):
				# 	sc = scenarios[s]
				# 	sc_df = loss_df[(loss_df.hazard_type == sc.hazard_type) & (loss_df.model == sc.model) & (loss_df.climate_scenario == sc.climate_scenario) & (loss_df.year == sc.year)]
				# 	if len(sc_df.index) > 0:
				# 		vals = sc_df['min_daily_loss_2016'].append(sc_df['max_daily_loss_2016'])
				# 		vals = vals.unique()
				# 		vals_list.append(1.0*vals/1000000)

				# 		label_list.append(sc.model_no)
				# 		label_color.append(sc.color)

				# fig, ax = plt.subplots(figsize=(8,4))
				# bp = ax.boxplot(vals_list, meanline=True,labels=label_list, showfliers=True,patch_artist = True)

				# for element in ['boxes','whiskers', 'means', 'medians', 'caps']:
				# 	for e, color in zip(bp[element], label_color):
				# 		plt.setp(e,color = '#252525',linewidth = 0.5)

				# for patch, color in zip(bp['boxes'], label_color):
				# 	patch.set_facecolor(color)

				# for flier,color in zip(bp['fliers'], label_color):
				# 	flier.set(marker='o', markersize=2, markerfacecolor=color,linestyle='none',markeredgecolor=color)


				# ax.tick_params(axis='x', rotation=45)
				# plt.xlabel('Hazard models', fontweight='bold')
				# plt.ylabel('Daily economic impacts (USD million/day)', fontweight='bold')
				# plt.title('{}: Hazard economic impacts'.format(province))
				# # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 8)

				# # print (bp['boxes'][0])
				# # fig_path = os.path.join(figure_path,'{0}_5_tons_hazard_economic_impacts.png'.format(province_name))
				# fig_path = os.path.join(figure_path,'{0}_hazard_economic_impacts.png'.format(province_name))
				# plt.tight_layout()
				# plt.savefig(fig_path,dpi=500)
				# plt.close()



if __name__ == "__main__":
	main()
