# -*- coding: utf-8 -*-
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
from vtra.transport_flow_and_failure_functions import *


def net_present_value(probability_wt, duration, adaptation_options_dataframe,
                      strategy_parameter, parameter_value_list, min_economic_loss,
                      max_economic_loss, edge_options_dictionary, start_year, end_year,
                      growth_rate, discount_rate, edge_width=1.0, edge_length=1.0):
   
    total_discount_ratio = []
    for year in range(start_year, end_year):
        # total_discount_ratio += 1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year)
        total_discount_ratio.append(
            1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))

    # total_discount_ratio = sum(total_discount_ratio_list)

    df_path = os.path.join(data_path, 'Adaptation_options', 'adaptation_options.xlsx')
    adaptation_df = pd.read_excel(df_path, sheet_name='adapt_options')

    adaptation_df['total_discount_ratio'] = sum(total_discount_ratio)

    min_maintain_discount_ratio_list = []
    max_maintain_discount_ratio_list = []

    for iter_, row in adaptation_df.iterrows():
        min_maintain_schedule = row['maintenance_times_min']
        max_maintain_schedule = row['maintenance_times_max']

        min_maintain_discount_ratio = 0
        max_maintain_discount_ratio = 0

        max_maintain_discount_years = np.arange(start_year, end_year, min_maintain_schedule)
        min_maintain_discount_years = np.arange(start_year, end_year, max_maintain_schedule)
        for year in max_maintain_discount_years[1:]:
            max_maintain_discount_ratio += 1.0 / \
                math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year)

        for year in min_maintain_discount_years[1:]:
            min_maintain_discount_ratio += 1.0 / \
                math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year)

        min_maintain_discount_ratio_list.append(min_maintain_discount_ratio)
        max_maintain_discount_ratio_list.append(max_maintain_discount_ratio)
    
    for param_val in parameter_value_list:

        st = adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                              == param_val, 'strategy'].values[0]

        cost_unit = adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                                     == param_val, 'cost_unit'].values[0]
        if cost_unit == 'million USD':
            cost_fact = 1000000.0
        else:
            cost_fact = 1.0

        adapt_cost_min = cost_fact*edge_width*edge_length * \
            adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                             == param_val, 'min_initial_cost'].values[0]
        adapt_cost_max = cost_fact*edge_width*edge_length * \
            adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                             == param_val, 'max_initial_cost'].values[0]

        st_min_rehab_benefit = cost_fact*probability_wt*edge_width*edge_length * \
            adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                             == param_val, 'min_rehab_benefit'].sum()
        st_max_rehab_benefit = cost_fact*probability_wt*edge_width*edge_length * \
            adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                             == param_val, 'max_rehab_benefit'].sum()

        st_min_cost = cost_fact*edge_width*edge_length * \
            adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                             == param_val, 'min_cost'].sum()
        st_max_cost = cost_fact*edge_width*edge_length * \
            adaptation_options_dataframe.loc[adaptation_options_dataframe[strategy_parameter]
                                             == param_val, 'max_cost'].sum()

        total_discount_ratio = []
        for year in range(start_year, end_year):
            # total_discount_ratio += 1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year)
            total_discount_ratio.append(1.0*math.pow(1.0 + 1.0*growth_rate/100.0, year -
                                                     start_year)/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))

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

        edge_options_dictionary.append({'adapt_strategy': st, 'min_initial_cost': adapt_cost_min, 'max_initial_cost': adapt_cost_max,
                                        'min_benefit_npv': st_min_benefit, 'max_benefit_npv': st_max_benefit,
                                        'min_cost_npv': st_min_cost, 'max_cost_npv': st_max_cost,
                                        'min_adapt_npv': min_npv, 'max_adapt_npv': max_npv, 'min_bc_ratio': min_bc_ratio, 'max_bc_ratio': max_bc_ratio})

    return edge_options_dictionary
