# -*- coding: utf-8 -*-
"""Estimate costs and benefits, either under fixed parameters or under a sensitivity analysis,
varying the cost components.
"""
import os
import sys
import pandas as pd
import numpy as np
import math
from SALib.sample import morris
from tqdm import tqdm

from dask import dataframe as dd
from dask.multiprocessing import get
from multiprocessing import cpu_count
from vtra.utils import *

nCores = cpu_count()


def calculate_discounting_arrays(discount_rate=12,growth_rate=6):
    """Set discount rates for yearly and period maintenance costs

    Parameters
    ----------
    discount_rate
        yearly discount rate
    growth_rate
        yearly growth rate

    Returns
    -------
    discount_rate_norm
        discount rates to be used for the costs
    discount_rate_growth
        discount rates to be used for the losses
    min_main_dr
        discount rates for 4-year periodic maintenance
    max_main_dr
        discount rates for 8-year periodic maintenance

    """
    discount_rate_norm = []
    discount_rate_growth = []

    for year in range(2016, 2050):
        discount_rate_norm.append(
            1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016))

        discount_rate_growth.append(
            1.0*math.pow(1.0 + 1.0*growth_rate/100.0, year -
                                                         2016)/math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016))

    min_maintain_discount_years = np.arange(2016, 2050, 4)
    maintain_discount_ratio_list = []
    for year in min_maintain_discount_years[1:]:
        maintain_discount_ratio = 1.0 / math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016)
        maintain_discount_ratio_list.append(maintain_discount_ratio)

    max_main_dr = np.array(maintain_discount_ratio_list)

    maintain_discount_ratio_list = []
    max_maintain_discount_years = np.arange(2016, 2050, 8)
    for year in max_maintain_discount_years[1:]:
        maintain_discount_ratio = 1.0 / math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016)
        maintain_discount_ratio_list.append(maintain_discount_ratio)

    min_main_dr = np.array(maintain_discount_ratio_list)

    return np.array(discount_rate_norm), np.array(discount_rate_growth), min_main_dr, max_main_dr

def sum_tuples(l):
    return list(sum(x) for x in zip(*l))

def average_tuples(l):
    return list(np.mean(x) for x in zip(*l))

def max_tuples(l):
    return list(np.max(x) for x in zip(*l))


def calc_costs(x, param_values, mnt_dis_cost, mnt_nat_cost, cst_dis_cost, cst_nat_cost,
               brdg_cost, pavement, mnt_main_cost, cst_main_cost, discount_rates,
               discount_growth_rates, rehab_costs, min_main_dr, max_main_dr, duration_max=10,
               min_exp=True, national=False, min_loss=True):
    """Estimate the total cost and benefits for a road segment. This function is used within a
    pandas apply

    Parameters
    ----------
    x
        a row from the road segment dataframe that we are considering
    param_values
        numpy array with a set of parameter combinations
    mnt_dis_cost
        adaptation costs for a district road in the mountains
    mnt_nat_cost
        adaptation costs for a national road in the mountains
    cst_dis_cost
        adaptation costs for a district road on flat terrain
    cst_nat_cost
        adaptation costs for a national road on flat terrain
    pavement
        set of paving combinations. This corresponds with the cost table and the param_values
    mnt_main_cost
        maintenance costs for roads in the mountains
    cst_main_cost
        maintenance costs for roads on flat terrain
    discount_rates
        discount rates to be used for the costs
    discount_growth_rates
        discount rates to be used for the losses
    rehab_costs
        rehabilitation costs after a disaster
    min_main_dr
        discount rates for 4-year periodic maintenance
    max_main_dr
        discount rates for 8-year periodic maintenance
    min_exp : bool, optional
        Specify whether we want to use the minimum or maximum exposure length. The default value is set to **True**
    national : bool, optional
        Specify whether we are looking at national roads. The default value is set to **False**
    min_loss : bool, optional
        Specify whether we want to use the minimum or maximum economic losses.  The default value is set to **True**

    Returns
    -------
    uncer_output : list
        outcomes for the initial adaptation costs of this road segment
    tot_uncer_output : list
        outcomes for the total adaptation costs of this road segment
    rel_share : list
        relative share of each factor in the initial adaptation cost of this road segment
    tot_rel_share : list
        relative share of each factor in the total adaptation cost of this road segment
    bc_ratio : list
        benefit cost ratios for this road segment

    """
    # Identify terrain type of the road
    if x.terrain.lower().strip() == 'mountain' or x.asset_type == 'Bridge':
        main_cost = mnt_main_cost
        ter_type = 'mountain'
    elif x.terrain.lower().strip() == 'flat':
        main_cost = cst_main_cost
        ter_type = 'flat'


    # Set which exposure length to use
    if min_exp == True:
        exp_length = x.min_exposure_length
    else:
        exp_length = x.max_exposure_length

    # Set which loss to use
    if min_loss == True:
        loss = x.min_econ_impact
        duration = duration_max*x.min_duration_wt
    else:
        loss = x.max_econ_impact
        duration = duration_max*x.max_duration_wt

    # Identify asset type, which is the main driver of the costs
    if (x.asset_type == 'Expressway') | ((national == True) & (x.road_class == 1)):
        rehab_cost = rehab_costs.loc[('Expressway',ter_type),'rate_m']
        rehab_corr = rehab_costs.loc[('Expressway',ter_type),'design_width']
        general_road_class = 'national'
    elif (x.asset_type == 'National roads') | ((national == True) & (x.road_class == 2)):
        rehab_cost = rehab_costs.loc[('National  2x Carriageway',ter_type),'rate_m']
        rehab_corr = rehab_costs.loc[('National  2x Carriageway',ter_type),'design_width']
        general_road_class = 'national'
    elif (x.asset_type == 'National roads') | ((national == True) & (x.road_class == 3)):
        rehab_cost = rehab_costs.loc[('National  1x Carriageway',ter_type),'rate_m']
        rehab_corr = rehab_costs.loc[('National  1x Carriageway',ter_type),'design_width']
        general_road_class = 'national'
    elif (x.asset_type == 'Provincial roads') | ((national == True) & (x.road_class == 4)):
        rehab_cost = rehab_costs.loc[('Provincial',ter_type),'rate_m']
        rehab_corr = rehab_costs.loc[('Provincial',ter_type),'design_width']
        general_road_class = 'district'
    elif ((x.asset_type == 'Urban roads/Named roads') | (x.asset_type == 'Boulevard')) | ((national == True) & (x.road_class == 5)):
        rehab_cost = rehab_costs.loc[('District',ter_type),'rate_m']
        rehab_corr = rehab_costs.loc[('District',ter_type),'design_width']
        general_road_class = 'district'
    elif (x.asset_type == 'Other roads') | ((national == True) & (x.road_class == 6)):
        rehab_cost = rehab_costs.loc[('Commune',ter_type),'rate_m']
        rehab_corr = rehab_costs.loc[('Commune',ter_type),'design_width']
        general_road_class = 'district'
    elif x.asset_type == 'Bridge':
        rehab_cost = rehab_costs.rate_m.max()
        rehab_corr = rehab_costs.design_width.max()
        general_road_class = 'bridge'
    else:
        rehab_cost = rehab_costs.rate_m.min()
        rehab_corr = rehab_costs.design_width.min()
        general_road_class = 'district'

    rehab_cost = (rehab_cost*x.width)/rehab_corr

    if (x.asset_type in ['Urban roads/Named roads','Other roads',
                         'Boulevard','Provincial roads']) & (x.terrain == 'mountain'):
        costs = mnt_dis_cost
        width_corr = 4.5
    elif (x.asset_type in ['Urban roads/Named roads','Other roads',
                           'Boulevard','Provincial roads']) & (x.terrain == 'flat'):
        costs = cst_dis_cost
        width_corr = 4.5
    elif ((x.asset_type in ['National roads', 'Expressway']) | (national == True)) & (x.terrain == 'mountain'):
        costs = mnt_nat_cost
        width_corr = 15
    elif ((x.asset_type in ['National roads', 'Expressway']) | (national == True)) & (x.terrain == 'flat'):
        costs = cst_nat_cost
        width_corr = 15
    elif x.asset_type == 'Bridge':
        # 'bridges'
        costs = brdg_cost
        width_corr = 4.5
    else:
        costs = mnt_dis_cost
        width_corr = 4.5


    # Identify costs for paved roads
    pav_rates = costs.loc['Pavement','rate_m']

    # Identify costs for drainage
    if x.terrain == 'mountain' or x.asset_type == 'Bridge':
        drain_rates = costs.loc[('Pavement Drain','DT'),'rate_m']
    else:
        drain_rates = 0

    # Identify all other costs, mainly dependent on whether it is a mountain or a flat road
    if x.asset_type == 'Bridge':
        EW1_rates = costs.loc[('Earthwork','EW1'),'rate_m']*2
        EW2_rates = costs.loc[('Earthwork','EW2'),'rate_m']*0
        SS1_rates = costs.loc[('Slope Protection','SS1'),'rate_m']*2
        SS2_rates = costs.loc[('Slope Protection','SS2'),'rate_m']*0
        concr_rates = costs.loc[('Slope Protection','SS3'),'rate_m']*0
        riverb_rates = costs.loc[('Slope Protection','SS4'),'rate_m']*2

    elif x.asset_type == 'Culvert':
        EW1_rates = 0
        EW2_rates = 0
        SS1_rates = 0
        SS2_rates = 0
        concr_rates = 0
        riverb_rates = 0

    elif x.terrain == 'mountain':
        EW1_rates = costs.loc[('Earthwork','EW1'),'rate_m']/2
        EW2_rates = costs.loc[('Earthwork','EW2'),'rate_m']/2
        SS1_rates = costs.loc[('Slope Protection','SS1'),'rate_m']*2
        SS2_rates = costs.loc[('Slope Protection','SS2'),'rate_m']*2
        concr_rates = costs.loc[('Slope Protection','SS3'),'rate_m']
        riverb_rates = costs.loc[('Slope Protection','SS4'),'rate_m']/10

    elif x.terrain == 'flat':
        EW1_rates = costs.loc[('Earthwork','EW1'),'rate_m']*2
        EW2_rates = costs.loc[('Earthwork','EW2'),'rate_m']*0
        SS1_rates = costs.loc[('Slope Protection','SS1'),'rate_m']*0
        SS2_rates = costs.loc[('Slope Protection','SS2'),'rate_m']*2
        concr_rates = costs.loc[('Slope Protection','SS3'),'rate_m']*2
        riverb_rates = costs.loc[('Slope Protection','SS4'),'rate_m']/5


    if x.asset_type == 'Culvert':
        face_dr_rates = 0
        bioeng_rates = 0
    else:
        face_dr_rates = costs.loc[('Slope Protection','S55'),'estimated _amount_ fraction']*costs.loc[('Slope Protection','S55'),'rate_m']
        bioeng_rates = costs.loc[('Slope Protection','SS6'),'estimated _amount_ fraction']*costs.loc[('Slope Protection','SS6'),'rate_m']

    costs_culvert = (np.ceil(exp_length/200)*costs.loc[('Culverts','CV1'),'rate_m']+
                     np.ceil(exp_length/1000)*costs.loc[('Culverts','CV2'),'rate_m'])


    drain_recur = main_cost.loc[('Pavement Drain','DT'),'rec_rate_m']
    pav_recur = main_cost.loc[('Pavement'),'rec_rate_m']
    EW1_recur = main_cost.loc[('Earthwork','EW1'),'rec_rate_m']
    EW2_recur = main_cost.loc[('Earthwork','EW2'),'rec_rate_m']
    SS1_recur = main_cost.loc[('Slope Protection','SS1'),'rec_rate_m']
    SS2_recur = main_cost.loc[('Slope Protection','SS2'),'rec_rate_m']
    concr_recur = main_cost.loc[('Slope Protection','SS3'),'rec_rate_m']
    riverb_recur = main_cost.loc[('Slope Protection','SS4'),'rec_rate_m']
    face_dr_recur = main_cost.loc[('Slope Protection','S55'),'rec_rate_m']
    bioeng_recur = main_cost.loc[('Slope Protection','SS6'),'rec_rate_m']

    drain_period = main_cost.loc[('Pavement Drain','DT'),'periodic_cost']
    pav_period = main_cost.loc[('Pavement'),'periodic_cost']
    EW1_period = main_cost.loc[('Earthwork','EW1'),'periodic_cost']
    EW2_period = main_cost.loc[('Earthwork','EW2'),'periodic_cost']
    SS1_period = main_cost.loc[('Slope Protection','SS1'),'periodic_cost']
    SS2_period = main_cost.loc[('Slope Protection','SS2'),'periodic_cost']
    concr_period = main_cost.loc[('Slope Protection','SS3'),'periodic_cost']
    riverb_period = main_cost.loc[('Slope Protection','SS4'),'periodic_cost']
    face_dr_period = main_cost.loc[('Slope Protection','S55'),'periodic_cost']
    bioeng_period = main_cost.loc[('Slope Protection','SS6'),'periodic_cost']


    # Estimate benefit
    benefit = (sum(loss*discount_growth_rates*x.risk_wt*duration)+sum(discount_rates*x.dam_wt*rehab_cost))

    # Prepare lists for output
    uncer_output = []
    tot_uncer_output = []
    rel_share = []
    tot_rel_share = []
    bc_ratio = []
    bc_diff = []

    # preset parameter values, use if not reading from file
    if isinstance(param_values, dict):
        assert x.terrain in ('mountain', 'flat'), 'Expected terrain to be either "mountain" or "flat", got "%s"' % x.terrain
        param_set = "{}_{}".format(general_road_class, x.terrain)
        param_values = [
            param_values[param_set]
        ]


    # Loop through all parameter combinations for the uncertainty and sensitivity analysis
    for param in param_values:
        if  x.road_cond.lower() == 'paved':
            cost_pav = sum(pav_rates*np.array([0]+list(pavement[int(param[7]-1)][1:]))*exp_length)
            tot_cost_pav = (sum(pav_rates*np.array([0]+list(pavement[int(param[7]-1)][1:]))*exp_length)
                        + sum(discount_rates*exp_length*pav_recur[0])
                        + sum(min_main_dr*exp_length*x.width*pav_period[0])
                        + sum(sum([discount_rates*z*exp_length for z in pav_recur[1:]]))
                        + sum(sum([max_main_dr*z*exp_length*x.width for z in pav_period[1:]])))
        else:
            cost_pav = sum(pav_rates*pavement[int(param[7]-1)]*exp_length)
            tot_cost_pav = (sum(pav_rates*pavement[int(param[7]-1)]*exp_length)
                        + sum(discount_rates*exp_length*pav_recur[0]*pavement[int(param[7]-1)][0])
                        + sum(min_main_dr*exp_length*x.width*pav_period[0]*pavement[int(param[7]-1)][0])
                        + sum(sum([discount_rates*z*exp_length for z in pav_recur[1:]]))
                        + sum(sum([max_main_dr*z*exp_length*x.width for z in pav_period[1:]])))

        # Calculate the initial adaptation cost
        cost_drain = (drain_rates*param[0]*exp_length)
        cost_EW1 = (EW1_rates*param[1]*exp_length)
        cost_EW2 = (EW2_rates*param[2]*exp_length)
        cost_SS1 = (SS1_rates*param[3]*exp_length)
        cost_SS2 = (SS2_rates*param[4]*exp_length)
        cost_concr = (concr_rates*param[5]*exp_length)
        cost_riverb = (riverb_rates*param[6]*exp_length)
        cost_face_drain = (face_dr_rates*exp_length)
        cost_bioeng = (bioeng_rates*exp_length)

        # Calculate the total adaptation cost (initial + recurring + periodic costs)
        tot_cost_drain = (drain_rates*param[0]*exp_length+ sum(discount_rates*(drain_recur*exp_length))
                    + sum(min_main_dr*(drain_period*exp_length)))
        tot_cost_EW1 = (EW1_rates*param[1]*exp_length + sum(discount_rates*(EW1_recur*exp_length))
                    + sum(min_main_dr*(EW1_period*exp_length)))
        tot_cost_EW2 = (EW2_rates*param[2]*exp_length + sum(discount_rates*(EW2_recur*exp_length))
                    + sum(min_main_dr*(EW2_period*exp_length)))
        tot_cost_SS1 = (SS1_rates*param[3]*exp_length + sum(discount_rates*(SS1_recur*exp_length))
                    + sum(min_main_dr*(SS1_period*exp_length)))
        tot_cost_SS2 = (SS2_rates*param[4]*exp_length + sum(discount_rates*(SS2_recur*exp_length))
                    + sum(min_main_dr*(SS2_period*exp_length)))
        tot_cost_concr = (concr_rates*param[5]*exp_length + sum(discount_rates*(concr_recur*exp_length))
                    + sum(min_main_dr*(concr_period*exp_length)))
        tot_cost_riverb = (riverb_rates*param[6]*exp_length + sum(discount_rates*(riverb_recur*exp_length))
                    + sum(min_main_dr*(riverb_period*exp_length)))
        tot_cost_face_drain = (face_dr_rates*exp_length + sum(discount_rates*(face_dr_recur*exp_length))
                    + sum(min_main_dr*(face_dr_period*exp_length)))
        tot_cost_bioeng = (bioeng_rates*exp_length + sum(discount_rates*(bioeng_recur*exp_length))
                    + sum(min_main_dr*(bioeng_period*exp_length)))

        # Sum everything
        tot_cost_sum = (tot_cost_pav+tot_cost_drain+tot_cost_EW1+tot_cost_EW2
                            +tot_cost_SS1+tot_cost_SS2+tot_cost_concr+tot_cost_riverb+tot_cost_face_drain+tot_cost_bioeng+costs_culvert)/width_corr*x.width
        tot_uncer_output.append(tot_cost_sum)
        cost_sum = (costs_culvert+cost_pav+cost_drain+cost_EW1+cost_EW2
                            +cost_SS1+cost_SS2+cost_concr+cost_riverb+cost_face_drain+cost_bioeng)/width_corr*x.width
        uncer_output.append(cost_sum)
        rel_share.append(list([cost_pav,cost_drain,cost_EW1,cost_EW2,
                            cost_SS1,cost_SS2,cost_concr,cost_riverb,cost_face_drain,cost_bioeng,costs_culvert]))
        tot_rel_share.append(list([tot_cost_pav,tot_cost_drain,tot_cost_EW1,tot_cost_EW2,
                            tot_cost_SS1,tot_cost_SS2,tot_cost_concr,tot_cost_riverb,tot_cost_face_drain,tot_cost_bioeng,costs_culvert]))

        # Calculate the benefit cost ratio
        bc_ratio.append(benefit/tot_cost_sum)
        bc_diff.append(benefit-tot_cost_sum)

    rel_share = np.array(sum_tuples(rel_share))/sum(np.array(sum_tuples(rel_share)))*100
    tot_rel_share = np.array(sum_tuples(tot_rel_share))/sum(np.array(sum_tuples(tot_rel_share)))*100

    return benefit, uncer_output, tot_uncer_output, rel_share, tot_rel_share, bc_ratio, bc_diff

def run_file(file_id,data_path,output_path,duration_max=10,discount_rate=12,growth_rate=6, read_from_file=False):

    print('{} started!'.format(file_id))

    # load cost file
    adapt_path = os.path.join(data_path,'Adaptation_options')
    mnt_dis_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='costs_district_mountain',
                                 usecols=7,index_col=[0,1]).fillna(0)
    mnt_dis_cost.columns = map(str.lower, mnt_dis_cost.columns)
    mnt_dis_cost['rate_m'] = mnt_dis_cost.factor*mnt_dis_cost.rate

    mnt_nat_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='costs_national_mountain',
                                 usecols=7,index_col=[0,1]).fillna(0)
    mnt_nat_cost.columns = map(str.lower, mnt_nat_cost.columns)
    mnt_nat_cost['rate_m'] = mnt_nat_cost.factor*mnt_nat_cost.rate

    cst_dis_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='costs_district_flat',
                                 usecols=7,index_col=[0,1]).fillna(0)
    cst_dis_cost.columns = map(str.lower, cst_dis_cost.columns)
    cst_dis_cost['rate_m'] = cst_dis_cost.factor*cst_dis_cost.rate

    cst_nat_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='costs_national_flat',
                                 usecols=10,index_col=[0,1]).fillna(0)
    cst_nat_cost.columns = map(str.lower, cst_nat_cost.columns)
    cst_nat_cost['rate_m'] = cst_nat_cost.factor*cst_nat_cost.rate

    brdg_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='bridges',
                                 usecols=10,index_col=[0,1]).fillna(0)
    brdg_cost.columns = map(str.lower, brdg_cost.columns)
    brdg_cost['rate_m'] = brdg_cost.factor*brdg_cost.rate


    # load maintenance costs
    mnt_main_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='maintenance_mountain',
                                 usecols=10,index_col=[0,1]).fillna(0)
    mnt_main_cost['rec_rate_m'] = mnt_main_cost.recurrent_cost*mnt_main_cost.recurrent_factor
    cst_main_cost = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='maintenance_flat',
                                 usecols=10,index_col=[0,1]).fillna(0)
    cst_main_cost['rec_rate_m'] = cst_main_cost.recurrent_cost*cst_main_cost.recurrent_factor

    # load rehab costs
    rehab_costs = pd.read_excel(os.path.join(adapt_path,'adaptation_costs_road_types.xlsx'),sheet_name='rehabilitation_costs',
                                 usecols=4,index_col=[0,1]).fillna(0)
    rehab_costs['rate_m'] = rehab_costs.basic_cost*0.001

    # load param values
    if read_from_file:
        param_values = [np.fromfile(os.path.join(adapt_path,'param_values.pkl'))[x:x+8] for x in range(0, len(np.fromfile(os.path.join(adapt_path,'param_values.pkl'))), 8)]
    else:
        param_values = {
            'district_mountain': [
                2, # drainage
                0.05, # ew1
                0.25, # ew2
                1.5, # ss1
                0.5, # ss2
                1, # ss3
                0.05, # ss4
                2 # pavement option

            ],
            'district_flat': [
                0, # drainage
                2, # ew1
                0, # ew2
                0, # ss1
                2, # ss2
                2, # ss3
                0.1, # ss4
                3 # pavement option
            ],
            'national_mountain': [
                2, # drainage
                0.2, # ew1
                1, # ew2
                1.5, # ss1
                0.5, # ss2
                1, # ss3
                0.05, # ss4
                1 # pavement option
            ],
            'national_flat': [
                0, # drainage
                2, # ew1
                0, # ew2
                0, # ss1
                2, # ss2
                2, # ss3
                0.1, # ss4
                4 # pavement option
            ],
            'bridge_mountain': [
                2, # drainage
                2, # ew1
                0, # ew2
                0, # ss1
                2, # ss2
                0, # ss3
                2, # ss4
                2 # pavement option
            ],
            'bridge_flat': [
                2, # drainage
                2, # ew1
                0, # ew2
                0, # ss1
                2, # ss2
                0, # ss3
                2, # ss4
                2 # pavement option
            ]
        }

    pavement = np.array([
        [1,0,0,0], #
        [0,0.75,0.2,0.05],
        [0,0.9,0,0.1],
        [0,2,0,0]
    ])

    # load provinces
    if file_id == 'national_road':
        nat = True
        prov_roads = pd.read_csv(os.path.join(output_path,'hazard_scenarios','{}_hazard_intersections_risks.csv'.format(file_id)))
        loss_roads = pd.read_csv(os.path.join(output_path,'failure_results','minmax_combined_scenarios','single_edge_failures_minmax_{}_100_percent_disrupt.csv'.format(file_id)))[['edge_id','min_econ_impact','max_econ_impact']]
    else:
        prov_roads = pd.read_csv(os.path.join(output_path,'hazard_scenarios','roads_hazard_intersections_{}_risks.csv'.format(file_id)))
        loss_roads = pd.read_csv(os.path.join(output_path,'failure_results','minmax_combined_scenarios','single_edge_failures_minmax_{}_5_tons_100_percent_disrupt.csv'.format(file_id)))[['edge_id','min_econ_impact','max_econ_impact']]
        nat = False

    prov_roads = prov_roads.merge(loss_roads)

    dr_norm,dr_growth,min_main_dr,max_main_dr = calculate_discounting_arrays(discount_rate,growth_rate)
    print (sum(dr_norm),sum(dr_growth),sum(min_main_dr),sum(max_main_dr))

    tqdm.pandas()

    prov_roads['min_benefit'],prov_roads['min_ini_adap_cost'],prov_roads['min_tot_adap_cost'], \
    prov_roads['min_ini_rel_share'],prov_roads['min_tot_rel_share'],prov_roads['min_bc_ratio'],prov_roads['min_bc_diff'] = zip(*prov_roads.progress_apply(lambda x: calc_costs(x,param_values,
                                mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,brdg_cost,pavement,
                                mnt_main_cost,cst_main_cost,dr_norm,dr_growth,rehab_costs,min_main_dr,max_main_dr,
                                duration_max=duration_max,min_exp=True,national=nat,min_loss=True),axis=1))


    prov_roads['max_benefit'],prov_roads['max_ini_adap_cost'],prov_roads['max_tot_adap_cost'], \
    prov_roads['max_ini_rel_share'],prov_roads['max_tot_rel_share'],prov_roads['max_bc_ratio'],prov_roads['max_bc_diff'] = zip(*prov_roads.progress_apply(lambda x: calc_costs(x,param_values,
                                mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,brdg_cost,pavement,
                                mnt_main_cost,cst_main_cost,dr_norm,dr_growth,rehab_costs,min_main_dr,max_main_dr,
                                duration_max=duration_max,min_exp=False,national=nat,min_loss=False),axis=1))

    if not read_from_file:
        prov_roads['max_tot_adap_cost'] = prov_roads.max_tot_adap_cost.apply(lambda item: max(item))
        prov_roads['min_tot_adap_cost'] = prov_roads.min_tot_adap_cost.apply(lambda item: max(item))
        prov_roads['min_bc_ratio'] = prov_roads['min_benefit']/prov_roads['max_tot_adap_cost']
        prov_roads['max_bc_ratio'] = prov_roads['max_benefit']/prov_roads['min_tot_adap_cost']

    if read_from_file:
        prov_roads.to_csv(os.path.join(output_path,'adaptation_results','output_adaptation_{}_{}_days_max_disruption.csv'.format(file_id,duration_max)))
    else:
        prov_roads.to_csv(os.path.join(output_path,'adaptation_results','output_adaptation_{}_{}_days_max_disruption_fixed_parameters.csv'.format(file_id,duration_max)))

if __name__ == '__main__':
    read_from_file = False
    if len(sys.argv) == 2:
        if sys.argv[1] == '--read_from_file':
            read_from_file = True

    if read_from_file:
        print('Reading param values from file')
    else:
        print('Running with fixed param values, --read_from_file to use file')

    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    duration_list = [10]
    discount_rate = 12
    growth_rate = 6.5
    regions = ['laocai','binhdinh','thanhhoa','national_road']
    for dur in duration_list:
        for file_id in regions:
            run_file(file_id,data_path,output_path,duration_max=dur,discount_rate=discount_rate,growth_rate=growth_rate, read_from_file=read_from_file)
