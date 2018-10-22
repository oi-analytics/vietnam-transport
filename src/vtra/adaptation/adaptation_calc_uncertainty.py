# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 11:28:18 2018

@author: cenv0574
"""


import pandas as pd
import numpy as np
from SALib.sample import morris
from tqdm import tqdm
import math 

import matplotlib.pyplot as plt
import SALib.analyze.morris


def sum_tuples(l):
    return list(sum(x) for x in zip(*l))


def calc_costs(x,param_values,mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,
               pavement,mnt_main_cost,cst_main_cost,discount_rates,min_main_dr,max_main_dr,min_exp=True,national=False):
    
    if x.terrain == 'mountain':
        main_cost = mnt_main_cost
    elif x.terrain == 'flat':
        main_cost = cst_main_cost
    
    if (x.asset_type in ['Urban roads/Named roads','Other roads',
                         'Boulevard','Provincial roads']) & (x.terrain == 'mountain'):
        costs = mnt_dis_cost
    elif (x.asset_type in ['Urban roads/Named roads','Other roads',
                           'Boulevard','Provincial roads']) & (x.terrain == 'flat'):
        costs = cst_dis_cost
    elif ((x.asset_type in ['National roads', 'Expressway']) | (national == True)) & (x.terrain == 'mountain'):
        costs = mnt_nat_cost
    elif ((x.asset_type in ['National roads', 'Expressway']) | (national == True)) & (x.terrain == 'flat'):
        costs = cst_nat_cost
    else:
        return [0]*len(param_values)
    
    if min_exp == True:
        exp_length = x.min_exposure_length
    else:
        exp_length = x.max_exposure_length
    
    road_cond = x.road_cond
    
    pav_rates = costs.loc['Pavement','rate_m']
    if x.terrain == 'mountain':
        drain_rates = costs.loc[('Pavement Drain','DT'),'rate_m']
    else:
        drain_rates = 0
        
    if x.terrain == 'mountain':
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

    face_dr_rates = costs.loc[('Slope Protection','S55'),'estimated _amount_ fraction']*costs.loc[('Slope Protection','S55'),'rate_m']
    bioeng_rates = costs.loc[('Slope Protection','SS6'),'estimated _amount_ fraction']*costs.loc[('Slope Protection','SS6'),'rate_m']
     
    costs_culvert = (np.ceil(exp_length/200)*costs.loc[('Culverts','CV1'),'rate_m']+
                     np.ceil(exp_length/1000)*costs.loc[('Culverts','CV2'),'rate_m'])

    pav_recur = main_cost.loc[('Pavement'),'rec_rate_m']
    EW1_recur = main_cost.loc[('Earthwork','EW1'),'rec_rate_m']
    EW2_recur = main_cost.loc[('Earthwork','EW1'),'rec_rate_m']
    SS1_recur = main_cost.loc[('Slope Protection','SS1'),'rec_rate_m']
    SS2_recur = main_cost.loc[('Slope Protection','SS2'),'rec_rate_m']
    concr_recur = main_cost.loc[('Slope Protection','SS3'),'rec_rate_m']
    riverb_recur = main_cost.loc[('Slope Protection','SS4'),'rec_rate_m']
    face_dr_recur = main_cost.loc[('Slope Protection','S55'),'rec_rate_m']
    bioeng_recur = main_cost.loc[('Slope Protection','SS6'),'rec_rate_m']

    pav_period = main_cost.loc[('Pavement'),'periodic_cost']
    EW1_period = main_cost.loc[('Earthwork','EW1'),'periodic_cost']
    EW2_period = main_cost.loc[('Earthwork','EW1'),'periodic_cost']
    SS1_period = main_cost.loc[('Slope Protection','SS1'),'periodic_cost']
    SS2_period = main_cost.loc[('Slope Protection','SS2'),'periodic_cost']
    concr_period = main_cost.loc[('Slope Protection','SS3'),'periodic_cost']
    riverb_period = main_cost.loc[('Slope Protection','SS4'),'periodic_cost']
    face_dr_period = main_cost.loc[('Slope Protection','S55'),'periodic_cost']
    bioeng_period = main_cost.loc[('Slope Protection','SS6'),'periodic_cost']
    
    uncer_output = []
    for param in param_values:
        # pavement costs, come up with adjusting if param value is first one
        if road_cond == 'Paved':
            cost_pav = (sum(pav_rates*np.array([0]+list(pavement[int(param[7]-1)][1:]))*exp_length)
                        + sum(min_main_dr*exp_length*x.width*pav_recur[0])
                        + sum(min_main_dr*exp_length*x.width*pav_period[0])
                        + sum(sum([max_main_dr*z*exp_length*x.width for z in pav_recur[1:]]))
                        + sum(sum([max_main_dr*z*exp_length*x.width for z in pav_period[1:]])))
        else:
            cost_pav = (sum(pav_rates*pavement[int(param[7]-1)]*exp_length)
                        + sum(min_main_dr*exp_length*x.width*pav_recur[0]*pavement[int(param[7]-1)][0])
                        + sum(min_main_dr*exp_length*x.width*pav_period[0]*pavement[int(param[7]-1)][0])
                        + sum(sum([max_main_dr*z*exp_length*x.width for z in pav_recur[1:]]))
                        + sum(sum([max_main_dr*z*exp_length*x.width for z in pav_period[1:]])))
                        
        cost_drain = drain_rates*param[0]*exp_length
        cost_EW1 = (EW1_rates*param[1]*exp_length + sum(discount_rates*(EW1_recur*exp_length))
                    + sum(min_main_dr*(EW1_period*exp_length)))
        cost_EW2 = (EW2_rates*param[2]*exp_length + sum(discount_rates*(EW2_recur*exp_length))
                    + sum(min_main_dr*(EW2_period*exp_length)))
        cost_SS1 = (SS1_rates*param[3]*exp_length + sum(discount_rates*(SS1_recur*exp_length))
                    + sum(min_main_dr*(SS1_period*exp_length)))
        cost_SS2 = (SS2_rates*param[4]*exp_length + sum(discount_rates*(SS2_recur*exp_length))
                    + sum(min_main_dr*(SS2_period*exp_length)))
        cost_concr = (concr_rates*param[5]*exp_length + sum(discount_rates*(concr_recur*exp_length))
                    + sum(min_main_dr*(concr_period*exp_length)))
        cost_riverb = (riverb_rates*param[6]*exp_length + sum(discount_rates*(riverb_recur*exp_length))
                    + sum(min_main_dr*(riverb_period*exp_length)))
        cost_face_drain = (face_dr_rates*exp_length + sum(discount_rates*(face_dr_recur*exp_length))
                    + sum(min_main_dr*(face_dr_period*exp_length)))
        cost_bioeng = (bioeng_rates*exp_length + sum(discount_rates*(bioeng_recur*exp_length))
                    + sum(min_main_dr*(bioeng_period*exp_length)))
        
        uncer_output.append((costs_culvert+cost_pav+cost_drain+cost_EW1+cost_EW2
                            +cost_SS1+cost_SS2+cost_concr+cost_riverb+cost_face_drain+cost_bioeng)/4.5*x.width)

    return uncer_output

if __name__ == '__main__':

    # load provinces
    prov_roads = pd.read_csv('roads_hazard_intersections_laocai_risks.csv')
    
    # load cost file
    costs = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='district_mountain',usecols=7,index_col=[0,1])
    costs['rate_m'] = costs.factor*costs.rate

    # load cost file
    mnt_dis_cost = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='costs_district_mountain',
                                 usecols=7,index_col=[0,1]).fillna(0)
    mnt_dis_cost['rate_m'] = mnt_dis_cost.factor*mnt_dis_cost.rate
    mnt_nat_cost = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='costs_national_mountain',
                                 usecols=7,index_col=[0,1]).fillna(0)
    mnt_nat_cost['rate_m'] = mnt_nat_cost.Factor*mnt_nat_cost.Rate
    cst_dis_cost = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='costs_district_flat',
                                 usecols=7,index_col=[0,1]).fillna(0)
    cst_dis_cost['rate_m'] = cst_dis_cost.Factor*cst_dis_cost.rate
    cst_nat_cost = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='costs_national_flat',
                                 usecols=10,index_col=[0,1]).fillna(0)
    cst_nat_cost['rate_m'] = cst_nat_cost.Factor*cst_nat_cost.Rate
    
    # load maintenance costs
    mnt_main_cost = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='periodic_maintenance_mountain',
                                 usecols=10,index_col=[0,1]).fillna(0)
    mnt_main_cost['rec_rate_m'] = mnt_main_cost.recurrent_cost*mnt_main_cost.recurrent_factor
    cst_main_cost = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='periodic_maintenance_flat',
                                 usecols=10,index_col=[0,1]).fillna(0)
    cst_main_cost['rec_rate_m'] = cst_main_cost.recurrent_cost*cst_main_cost.recurrent_factor

    # set discount rates for yearly and period maintenance costs
    discount_rate_norm = []
    discount_rate_growth = []
    
    discount_rate = 12
    growth_rate = 6
    for year in range(2016, 2050):
        discount_rate_norm.append(
            1.0/math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016))
        
        discount_rate_growth.append(
            1.0*math.pow(1.0 + 1.0*growth_rate/100.0, year -
                                                         2016)/math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016))
        
    dr_norm = np.array(discount_rate_norm)
    dr_growth = np.array(discount_rate_growth)
    
    min_maintain_discount_years = np.arange(2016, 2050, 4)
    maintain_discount_ratio = 0
    maintain_discount_ratio_list = []
    for year in min_maintain_discount_years[1:]:
        maintain_discount_ratio += 1.0 / math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016)
        maintain_discount_ratio_list.append(maintain_discount_ratio)
        
    min_main_dr = np.array(maintain_discount_ratio_list)
    
    max_maintain_discount_years = np.arange(2016, 2050, 8)
    maintain_discount_ratio = 0
    maintain_discount_ratio_list = []
    for year in max_maintain_discount_years[1:]:
        maintain_discount_ratio += 1.0 / math.pow(1.0 + 1.0*discount_rate/100.0, year - 2016)
        maintain_discount_ratio_list.append(maintain_discount_ratio)
        
    max_main_dr = np.array(maintain_discount_ratio_list)    

    # set up sensitivity
    problem = {
          'num_vars': 8,
          'names': ['drain', 'EW1','EW2','SS1','SS2','concrete','riverb','pavement'],
          'bounds': [[0,2],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[1,4]]}
    
    param_values = morris.sample(problem, 10, num_levels=4, grid_jump=2,local_optimization =True)
    
    pavement = np.array([[1,0,0,0],[0,0.75,0.2,0.05],[0,0.9,0,0.1],[0,2,0,0]])

    tqdm.pandas(desc='Minimum Cost')
    prov_roads['min_adap_cost'] = prov_roads.progress_apply(lambda x: calc_costs(x,param_values,
                                    mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,pavement,
                                    mnt_main_cost,cst_main_cost,dr_norm,min_main_dr,max_main_dr,min_exp=True,national=True),axis=1)
    tqdm.pandas(desc='Maximum Cost')
    prov_roads['max_adap_cost'] = prov_roads.progress_apply(lambda x: calc_costs(x,param_values,
                                    mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,pavement,
                                    mnt_main_cost,cst_main_cost,dr_norm,min_main_dr,max_main_dr,min_exp=False,national=True),axis=1)    # analyze uncertainty and sensitivity of results

    Si = SALib.analyze.morris.analyze(problem, np.array(param_values), 
                                              np.array(sum_tuples(list(prov_roads['max_adap_cost']))),
                                         print_to_console=False, grid_jump=2, num_levels=4)
    
    risk_sens = pd.DataFrame.from_dict(Si)
    risk_sens['rel'] = risk_sens['mu']/risk_sens['mu'].sum()*100
    risk_sens = risk_sens.groupby('names').sum()
    risk_sens = risk_sens.T
    
    stats=risk_sens.loc['rel',np.array(risk_sens.columns)].values

    # plot results    
    fig, ax = plt.subplots(1, 1,figsize=(6,6),subplot_kw=dict(projection='polar')) #
    
    angles=np.linspace(0, 2*np.pi, len(np.array(risk_sens.columns)), endpoint=False)
    # close the plot
    stats=np.concatenate((stats,[stats[0]]))
    angles=np.concatenate((angles,[angles[0]]))
    
    fig=plt.figure()
    ax.plot(angles, stats, 'o-', linewidth=2)
    ax.set_ylim([0, 50])   
    ax.set_yticklabels([])
    ax.fill(angles, stats, alpha=0.25)
    ax.set_thetagrids(angles * 180/np.pi, np.array(risk_sens.columns))
    ax.tick_params(axis='x',labelsize=10,labelcolor='black',color='black') # pad=12
    
    ax.set_title('Lao Cai Sensitivity Analysis', y=1.1,fontweight='bold')
    ax.grid(True)
    
    
