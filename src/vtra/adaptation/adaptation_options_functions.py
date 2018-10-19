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


def net_present_value(adaptation_options_file,
                      min_economic_loss,
                      max_economic_loss, start_year, end_year,
                      growth_rate, discount_rate,road_class,road_width,road_terrain='flat',edge_width=1.0, edge_length=1.0):
   
    """Find min-max economic losses
    """
    total_discount_ratio = []
    for year in range(start_year, end_year):
        total_discount_ratio.append(1.0*math.pow(1.0 + 1.0*growth_rate/100.0, year - \
            start_year)/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))

    min_eal = sum(min_economic_loss*np.array(total_discount_ratio))
    max_eal = sum(max_economic_loss*np.array(total_discount_ratio))

    """Find the damage costs
    """
    total_discount_ratio = []
    for year in range(start_year, end_year):
        total_discount_ratio.append(
            1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year))


    rehab_costs_df = pd.read_excel(adaptation_options_file,sheet_name='baseline_costs',thousands=',').fillna(0)
    rehab_costs = list(rehab_costs_df.itertuples(index=False))
    del rehab_costs_df
    for rc in rehab_costs:
        if road_terrain in rc.terrain and rc.road_class == road_class:
            min_ead = sum((1.0*road_width/rc.design_width)*float(rc.basic_cost)*np.array(total_discount_ratio))
            max_ead = min_ead
            break

    del rehab_costs

    """Get the rehab costs
    """
    rehab_costs = list(rehab_costs_df.itertuples(index=False))
    del rehab_costs_df
    for rc in rehab_costs:
        if road_terrain in rc.terrain and rc.road_class == road_class:
            min_ead = sum((1.0*road_width/rc.design_width)*float(rc.basic_cost)*np.array(total_discount_ratio))
            max_ead = min_ead
            break

    del rehab_costs

    """Create options lists
    """
    options_df = pd.read_excel(adaptation_options_file,sheet_name='options',thousands=',').fillna(0)
    options_list = []
    """Look at pavement options and their fractions
    """
    options_codes = options_df.loc[options_df['adaption_group'] == 'Pavement','option_code'].values.tolist()
    if road_cond == 'paved':
        options_frac = [0] + list(np.random.dirichlet(np.ones(3)))
    else:
        options_frac = list(np.random.dirichlet(np.ones(4)))

    options_list += list(zip(options_codes,options_frac))
    
    """Look at pavement drainage options
    """
    options_codes = options_df.loc[options_df['adaption_group'] == 'Pavement Drain','option_code'].values.tolist()
    DL = np.random.uniform(low=0,high=1)
    DR = np.random.uniform(low=0,high=1)
    DT = DL + DR
    options_frac = [DT,DR,DT]
    options_list += list(zip(options_codes,options_frac))

    """Look at Earthwork
    """
    options_codes = options_df.loc[options_df['adaption_group'] == 'Earthwork','option_code'].values.tolist()
    EW1 = np.random.uniform(low=0,high=1)
    EW2 = np.random.uniform(low=0,high=1)
    options_frac = [EW1,EW2]
    options_list += list(zip(options_codes,options_frac))

    """Look at Slope Protection
    """
    options_codes = options_df.loc[options_df['adaption_group'] == 'Slope Protection','option_code'].values.tolist()
    for oc in options_codes:
        if oc == 'SS5':
            options_list.append((oc,np.random.uniform(low=0,high=3)))
        else:
            options_list.append((oc,np.random.uniform(low=0,high=1)))

    """Look at culverts
    """
    options_codes = options_df.loc[options_df['adaption_group'] == 'Culverts','option_code'].values.tolist()
    for oc in options_codes:
        if oc == 'CV1':
            options_list.append((oc,np.ceil(1.0*edge_length/200)))
        if oc == 'CV2':
            options_list.append((oc,np.ceil(1.0*edge_length/1000)))

    options_df = pd.DataFrame(options_list,columns=['option_code','frac'])
    """Maintenance schedule
    """
    maintain_times_df = pd.read_excel(adaptation_options_file,sheet_name='maintenance_schedule',thousands=',').fillna(0)
    maintain_times_df = pd.merge(maintain_times_df,options_df,on=['option_code'], how='left').fillna(0)
    maintain_times_df['min_maintain_time'] = np.ceil(maintain_times_df['frac'])*maintain_times_df['min_maintain_time']
    maintain_times_df['max_maintain_time'] = np.ceil(maintain_times_df['frac'])*maintain_times_df['max_maintain_time']

    
    maintain_times = list(maintain_times_df.itertuples(index=False))
    del maintain_times_df
    for mt in maintain_times:
        for 
        if road_terrain in rc.terrain and rc.road_class == road_class:
            min_ead = sum((1.0*road_width/rc.design_width)*float(rc.basic_cost)*np.array(total_discount_ratio))
            max_ead = min_ead
            break

    del rehab_costs

    min_maintain_schedule = 8
    min_maintain_schedule = 8

    min_maintain_discount_ratio = 0
    max_maintain_discount_ratio = 0

    max_maintain_discount_years = np.arange(start_year, end_year, min_maintain_schedule)
    min_maintain_discount_years = np.arange(start_year, end_year, min_maintain_schedule)
    for year in max_maintain_discount_years[1:]:
        max_maintain_discount_ratio += 1.0 / \
            math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year)

    for year in min_maintain_discount_years[1:]:
        min_maintain_discount_ratio += 1.0 / \
            math.pow(1.0 + 1.0*discount_rate/100.0, year - start_year)


    """Cost of inital investment
    """
    adaptation_df = pd.read_excel(df_path, sheet_name='adapt_options')    
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


def combine_hazards_and_network_attributes_and_impacts(hazard_dataframe,network_dataframe):
    hazard_dataframe.loc[hazard_dataframe['probability']
                                == 'none', 'probability'] = 1.0
    hazard_dataframe['probability'] = pd.to_numeric(
        hazard_dataframe['probability'])
    hazard_dataframe.rename(columns={'length': 'exposure_length','min_val':'min_flood_depth','max_val':'max_flood_depth'}, inplace=True)

    network_dataframe.rename(columns={'length': 'road_length'}, inplace=True)
    network_dataframe['road_length'] = 1000.0*network_dataframe['road_length']

    all_edge_fail_scenarios = pd.merge(hazard_dataframe, network_dataframe, on=[
                                           'edge_id'], how='left').fillna(0)

    all_edge_fail_scenarios['percent_exposure'] = 100.0 * \
            all_edge_fail_scenarios['exposure_length']/all_edge_fail_scenarios['road_length']

    all_edge_fail_scenarios.to_csv('test.csv',index = False)

    del hazard_dataframe, network_dataframe

    return all_edge_fail_scenarios

def create_hazard_scenarios_for_adaptation(all_edge_fail_scenarios,index_cols,length_thr):
    all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(index_cols)
    scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
    print('Number of failure scenarios', len(scenarios))
    scenarios_list = []
    for sc in scenarios:
        min_height = max(all_edge_fail_scenarios.loc[[sc], 'min_flood_depth'].values.tolist())
        max_height = max(all_edge_fail_scenarios.loc[[sc], 'max_flood_depth'].values.tolist())
        min_band_num = min(all_edge_fail_scenarios.loc[[sc], 'band_num'].values.tolist())
        max_band_num = max(all_edge_fail_scenarios.loc[[sc], 'band_num'].values.tolist())
        prob = all_edge_fail_scenarios.loc[[sc], 'probability'].values
        if len(list(set(prob))) > 1:
            exposure_len = all_edge_fail_scenarios.loc[[sc], 'exposure_length'].values
            per = all_edge_fail_scenarios.loc[[sc], 'percent_exposure'].values

            prob_tup = list(zip(prob, exposure_len, per))
            u_pr = sorted(list(set(prob.tolist())))
            exposure_len = []
            per = []
            r_wt = []
            for pr in u_pr:
                per_exp = sum([z for (x, y, z) in prob_tup if x == pr])
                exp_len = sum([y for (x, y, z) in prob_tup if x == pr])
                if per_exp > 100.0:
                    exposure_len.append(100.0*exp_len/per_exp)
                    per.append(100.0)
                    r_wt.append(1.0)
                else:
                    exposure_len.append(exp_len)
                    per.append(per_exp)
                    if exp_len < length_thr:
                        r_wt.append(0.01*per_exp)
                    else:
                        r_wt.append(1.0)

            max_exposure_len = max(exposure_len)
            min_exposure_len = min(exposure_len)

            min_per = min(per)
            max_per = max(per)
            min_dur = 0.01*min_per
            max_dur = 0.01*max_per
            risk_wt = 0
            dam_wt = 0
            for p in range(len(u_pr)-1):
                risk_wt += 0.5*(u_pr[p+1]-u_pr[p])*(r_wt[p+1]+r_wt[p])
                dam_wt +=  0.5*(u_pr[p+1]-u_pr[p])*(exposure_len[p+1]+exposure_len[p])

        else:
            prob_wt = prob[0]
            min_exposure_len = sum(
                all_edge_fail_scenarios.loc[[sc], 'exposure_length'].values.tolist())
            min_per = sum(all_edge_fail_scenarios.loc[[
                          sc], 'percent_exposure'].values.tolist())
            if min_per > 100.0:
                min_exposure_len = 100.0*min_exposure_len/min_per
                min_per = 100.0

            max_per = min_per
            max_exposure_len = min_exposure_len
            dam_wt = max_exposure_len
            min_dur = 0.01*min_per
            if max_exposure_len < length_thr:
                max_dur = 0.01*max_per
                risk_wt = 0.01*max_per*prob_wt 
            else:
                max_dur = 1.0
                risk_wt = prob_wt

        scenarios_list.append(list(sc) + [min_band_num, max_band_num, min_height, max_height,
                                          min_per, max_per, min_dur, max_dur, min_exposure_len, max_exposure_len, risk_wt, dam_wt])

    new_cols = ['min_band', 'max_band', 'min_height', 'max_height', 'min_exposure_percent', 'max_exposure_percent',
                'min_duration_wt', 'max_duration_wt', 'min_exposure_length', 'max_exposure_length', 'risk_wt','dam_wt']
    scenarios_df = pd.DataFrame(scenarios_list, columns=index_cols + new_cols)
    
    del all_edge_fail_scenarios, scenarios_list
    return scenarios_df
