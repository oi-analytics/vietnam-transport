# -*- coding: utf-8 -*-
"""
Created on Tue May  1 08:15:32 2018

@author: cenv0574
"""

import os
import sys
import json

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config

def create_disruption(input_file):
    
    data_path = load_config()['paths']['data']

    # read commodity descriptions
    comm_des = pd.read_excel(os.path.join(data_path,'Results','commodity_descriptions_and_totals.xlsx'),sheet_name='Sheet1')
    
    # read vietnam mpof domestic results
    vnm_mpof_dom = pd.read_csv(input_file)
    
    # get industry daily tonnage totals
    ind_des = comm_des.copy()
    ind_des['industry'] = ind_des.commodity_field.apply(map_comm_ind)
    ind_des['ind'] = ind_des.industry.apply(map_ind)
    ind_des = ind_des.drop(['commodity_name','commodity_field'],axis=1)
    ind_des = ind_des.groupby('ind').sum()
    
    # convert commodities to industries for EORA impact modelling
    vnm_mpof_dom['industry'] = vnm_mpof_dom.industry.apply(map_comm_ind)
    vnm_mpof_dom['ind'] = vnm_mpof_dom.industry.apply(map_ind)
    
    # create scenarios per region
    regional_scenarios = vnm_mpof_dom.groupby(['region_flooded','typhoon_level']).sum()
      
    # get scenarios and industries for typhoons and regions
    typhoon_region_unique = vnm_mpof_dom.groupby(['region_flooded','typhoon_level','ind']).sum()
    
    #  create disruption events for typhoon per region
    typhoon_region_events = create_events(list(regional_scenarios.index),typhoon_region_unique,ind_des)
    
    event_dict = {}
    for event,edata in typhoon_region_events.groupby(axis=0,level=1):
        edata.index = edata.index.droplevel(1)
        edata.index = [x.replace(' ','_') for x in list(edata.index)]
        event_dict[event] = edata
        
    return event_dict

def create_events(scenario_list,flow_data,ind_des):

    events = {}
   
    for event in scenario_list:
        get_event  = flow_data.loc[event]
        get_event = get_event.merge(ind_des, left_index=True,right_index=True)
        get_event['tonnage'] = get_event['tonnage'].div(get_event['total_daily_tonnage'])
        events[event] = get_event['tonnage']     
        
    return pd.DataFrame(events).fillna(0).T

def map_comm_ind(x):

    comm_ind_map = { 'sugar' : 'Agriculture',
    'wood' : 'Agriculture',
    'steel' : 'Processing',
    'constructi' : 'Construction',
    'cement' : 'Mining',
    'fertilizer' : 'Processing',
    'coal' : 'Mining',
    'petroluem' : 'Processing',
    'manufactur' : 'Manufacturing',
    'fishery' : 'Agriculture',
    'meat' : 'Processing',
    'rice' : 'Agriculture',
    'cash' : 'Agriculture',
    'cass' : 'Agriculture',
    'teas' : 'Processing',
    'maiz' : 'Agriculture',
    'rubb' : 'Manufacturing',
    'swpo' : 'Agriculture',
    'acof' : 'Processing',
    'rcof' : 'Processing',
    'pepp' : 'Agriculture' }

    return comm_ind_map[x]

def map_ind(x):
    
    ind_map = {'secA':'Agriculture',
           'secB':'Mining',
           'secC':'Processing',
           'secD':'Textile and Garment',
           'secE':'Wood and Paper',
           'secF':'Manufacturing',
           'secG':'Construction',
           'secH':'Trade',
           'secI':'Services'}
    
    return  {v: k for k, v in ind_map.items()}[x]


if __name__ == "__main__":
    # Define current directory and data directory
    data_path = load_config()['paths']['data']
   
    input_file = os.path.join(data_path,'Results','vnm_road_rail_edge_multi_failure.csv')

