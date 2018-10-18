"""
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

import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from vtra.utils import *
from vtra.adaptation.adaptation_options_functions import *



def main():
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    start_year = 2016
    end_year = 2050
    discount_rate = 12.0
    percentage = [100.0]
    duration_max = [10, 15, 20, 25, 30]
    length_thr = 500.0
    single_edge = True

    cols = ['band_num', 'climate_scenario',
            'edge_id', 'hazard_type', 'max_val', 'min_val', 'model', 'probability',
            'year', 'length']
    growth_scenarios = [(5, 'low'), (6.5, 'forecast'), (10, 'high')]
    base_year = 2016
    types = ['min', 'max']
    index_cols = ['edge_id', 'hazard_type', 'model', 'climate_scenario', 'year', 
        'terrain','surface','road_class','road_cond','asset_type','width','road_length']
    

    """Give the paths to the input data files
    """
    network_loss_data = os.path.join(data_path,'failure_results')
    fail_scenarios_data = os.path.join(
        output_path, 'hazard_scenarios')

    adaptation_data_path = os.path.join(data_path, 'Adaptation_options', 'adaptation_costs_road_types.xlsx')
    
    """Specify the output files and paths to be created
    """
    adapt_output_path = os.path.join(output_path, 'adaptation_results')
    if os.path.exists(adapt_output_path) == False:
        os.mkdir(adapt_output_path)
    modes = ['road']
    for m in range(len(modes)):
        scenarios_df = pd.read_csv(os.path.join(fail_scenarios_data,
            'national_{}_hazard_intersections_risks.csv'.format(modes[m])))

        for perct in percentage:
            """Load failure impact results
            """
            impacts_df = pd.read_csv(os.path.join(network_loss_data, 
                'single_edge_failures_minmax_national_{0}_{1}_percent_disrupt.csv'.format(mode,int(perct))))

            impacts_df = impacts_df[['edge_id','min_econ_impact','max_econ_impact']]
            impacts_df.rename(columns={'min_econ_impact': 'min_daily_loss_{}'.format(start_year),
                'max_econ_impact':'max_daily_loss_{}'.format(start_year)}, inplace=True)

            scenarios_df = pd.merge(scenarios_df,impacts_df,on=['edge_id'], how='left').fillna(0)

            print('done with scenarios creation')

            scenarios_df = scenarios_df[scenarios_df['max_daily_loss_{}'.format(start_year)] > 0]

            if single_edge == True:
                    file_name = 'single_edge_failures_scenarios_national_{0}_{1}_percent_disrupt_adapt_options.xlsx'.format(modes[m],int(perct))
                else:
                    file_name = 'multiple_edge_failures_scenarios_national_{0}_{1}_percent_disrupt_adapt_options.xlsx'.format(modes[m],int(perct))

            adapt_output_excel = os.path.join(adapt_output_path,file_name)
            excl_wrtr = pd.ExcelWriter(adapt_output_excel)

            for gr in range(len(growth_scenarios)):
                grth = growth_scenarios[gr]
                for dur in range(len(duration_max)):
                    scenarios_list = []
                    scenarios_df['min_duration'] = duration_max[dur]*scenarios_df['min_duration_wt']
                    scenarios_df['max_duration'] = duration_max[dur]*scenarios_df['max_duration_wt']
                    index_cols = scenarios_df.columns.values.tolist()
                    for iter_, row in scenarios_df.iterrows():
                        edge_options = []
                        min_eael = row['risk_wt']*row['min_duration']*row['min_daily_loss_{}'.format(start_year)]
                        max_eael = row['risk_wt']*row['max_duration']*row['max_daily_loss_{}'.format(start_year)]
                        road_cond = row['road_cond']
                        road_terrian = row['terrain']
                        road_level = row['class']
                        width = row['width']
                        max_exposure_len = row['max_exposure_length']

                        if max_height > 4:
                            max_height = 6
                        if max_height > 0:
                            edge_options = net_present_value(prob_wt, duration_max[dur], adaptation_options, 'height_m', [max_height], min_loss, max_loss, edge_options,
                                                             start_year, end_year, grth[0], discount_rate, edge_width=1.0, edge_length=0.001*max_exposure_len)

                        if road_cond == 'unpaved':
                            edge_options = net_present_value(prob_wt, duration_max[dur], adaptation_options, 'asset_cond', ['unpaved', 'unpaved-paved'], min_loss, max_loss,
                                                             edge_options, start_year, end_year, grth[0], discount_rate, edge_width=width, edge_length=max_exposure_len)

                        if road_cond == 'paved':
                            edge_options = net_present_value(prob_wt, duration_max[dur], adaptation_options, 'asset_cond', ['paved'], min_loss, max_loss, edge_options, start_year,
                                                             end_year, grth[0], discount_rate, edge_width=width, edge_length=max_exposure_len)

                        for options in edge_options:
                            scenarios_list.append({**dict(row), **options})

                    new_cols = ['adapt_strategy', 'min_initial_cost', 'max_initial_cost', 'min_benefit_npv', 'max_benefit_npv', 'min_cost_npv',
                                'max_cost_npv', 'min_adapt_npv', 'max_adapt_npv', 'min_bc_ratio', 'max_bc_ratio']
                    scenarios_list_df = pd.DataFrame(scenarios_list, columns=index_cols + new_cols)
                    scenarios_list_df.to_excel(
                        excl_wrtr, grth[1] + '_' + str(duration_max[dur]), index=False)
                    excl_wrtr.save()

                    print('Done with national {0} in {1} growth scenario {2} days'.format(
                        modes[m], grth[1], duration_max[dur]))


if __name__ == "__main__":
    main()
