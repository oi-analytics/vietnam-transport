# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 11:28:18 2018

@author: cenv0574
"""


import pandas as pd
import numpy as np
from SALib.sample import morris
from tqdm import tqdm
    
import matplotlib.pyplot as plt
import SALib.analyze.morris

from vtra.utils import *

def sum_tuples(l):
    return list(sum(x) for x in zip(*l))


def calc_costs(x,costs,param_values,mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,min_exp=True):
    
    if (x.asset_type in ['Urban roads/Named roads','Other roads',
                         'Boulevard','Provincial roads']) & (x.terrain == 'mountain'):
        costs = mnt_dis_cost
    elif (x.asset_type in ['Urban roads/Named roads','Other roads',
                           'Boulevard','Provincial roads']) & (x.terrain == 'flat'):
        costs = cst_dis_cost
    elif (x.asset_type in ['National roads', 'Expressway']) & (x.terrain == 'mountain'):
        costs = mnt_nat_cost
    elif (x.asset_type in ['National roads', 'Expressway']) & (x.terrain == 'flat'):
        costs = cst_nat_cost
    
    if min_exp == True:
        exp_length = x.min_exposure_length
    else:
        exp_length = x.max_exposure_length
    
    road_cond = x.road_cond
    
    pav_rates = costs.loc['Pavement','rate_m']
    drain_rates = costs.loc[('Pavement Drain','DT'),'rate_m']
    embank_rates = costs.loc[('Earthwork','EW1'),'rate_m']+costs.loc[('Slope Protection','SS2'),'rate_m']
    cut_slope_rates = costs.loc[('Earthwork','EW1'),'rate_m']+costs.loc[('Slope Protection','SS2'),'rate_m']
    concr_rates = costs.loc[('Slope Protection','SS3'),'rate_m']
    riverb_rates = costs.loc[('Slope Protection','SS4'),'rate_m']
    face_dr_rates = costs.loc[('Slope Protection','S55'),'rate_m']
    bioeng_rates = costs.loc[('Slope Protection','SS6'),'rate_m']
     
    costs_culvert = (np.ceil(exp_length/200)*costs.loc[('Culverts','CV1'),'rate_m']+
                     np.ceil(exp_length/1000)*costs.loc[('Culverts','CV2'),'rate_m'])

    uncer_output = []
    for param in param_values:
        # pavement costs
        if road_cond == 'Paved':
            cost_pav = sum(pav_rates*costs.loc['Pavement','estimated _amount_ fraction']*exp_length)
        else:
            cost_pav = costs.loc[('Pavement','OO'),'rate_m']*exp_length
        
        cost_drain = drain_rates*param[0]*exp_length
        cost_embank = embank_rates*param[1]*exp_length
        cost_slope = cut_slope_rates*param[2]*exp_length
        cost_concr = concr_rates*param[3]*exp_length
        cost_riverb = riverb_rates*param[4]*exp_length
        cost_face_drain = face_dr_rates*param[5]*exp_length
        cost_bioeng = bioeng_rates*param[6]*exp_length
        
        uncer_output.append((cost_pav+cost_drain+cost_embank+cost_slope
                            +costs_culvert+cost_concr+cost_riverb+cost_face_drain+cost_bioeng)/4.5*x.width)

    return uncer_output

if __name__ == '__main__':

    # load provinces
    prov_roads = pd.read_csv('roads_hazard_intersections_laocai_risks.csv')
    
    # load cost file
    costs = pd.read_excel('adaptation_costs_road_types.xlsx',sheet_name='district_mountain',usecols=7,index_col=[0,1])
    costs['rate_m'] = costs.factor*costs.rate

    # define all different costs (now all the same)    
    mnt_dis_cost = costs 
    mnt_nat_cost = costs
    cst_dis_cost = costs
    cst_nat_cost = costs
    
    # set up sensitivity
    problem = {
          'num_vars': 7,
          'names': ['drain', 'embank','cut_slope','concrete','riverb','face_drain','bioeng'],
          'bounds': [[0,2],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]}
    
    param_values = morris.sample(problem, 10, num_levels=4, grid_jump=2,local_optimization =True)

    # run minimum and maximum adaptation cost calculation
    tqdm.pandas(desc='Minimum Cost')
    prov_roads['min_adap_cost'] = prov_roads.progress_apply(lambda x: calc_costs(x,costs,param_values,
                                                                                 mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,min_exp=True),axis=1)
    tqdm.pandas(desc='Maximum Cost')
    prov_roads['max_adap_cost'] = prov_roads.progress_apply(lambda x: calc_costs(x,costs,param_values,
                                                                                 mnt_dis_cost,mnt_nat_cost,cst_dis_cost,cst_nat_cost,min_exp=False),axis=1)
        
    # analyze uncertainty and sensitivity of results
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
    
    
