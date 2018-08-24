# -*- coding: utf-8 -*-
"""
Created on Tue May  1 08:15:32 2018

@author: cenv0574
"""

import os
import sys

import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config

def create_disruption(input_file,output_dir,min_rice=True,single_point=True):
    
     # Define current directory and data directory
    data_path = load_config()['paths']['data']
   
     # read national IO
    comm_des = df_com_to_ind(pd.read_excel(os.path.join(data_path,'OD_data','national_scale_od_matrix.xlsx'),sheet_name='total'),min_rice=min_rice)

   # read vietnam failure results
    vnm_failure_all_failure = pd.read_csv(input_file,index_col=[0])

    if single_point == False:
        mapper = dict(zip(vnm_failure_all_failure.index.unique(),['multi_{}'.format(n) for n in np.arange(1,len(vnm_failure_all_failure.index.unique())+1,1)]))
        vnm_failure_all_failure.index = vnm_failure_all_failure.index.map(lambda x: mapper[x])
        pd.DataFrame.from_dict(mapper,orient='index').to_csv(os.path.join(output_dir,'mapper.csv'))

    event_dict = {}
    for id_,scenario in vnm_failure_all_failure.groupby(level=0):
        if id_ == 'multi_102':
            break
        disruption = df_com_to_ind(scenario,min_rice=min_rice).div(comm_des).dropna(axis=0, how='all').fillna(0)
        disruption = pd.DataFrame(disruption.stack(0))
        disruption.columns = ['value']
        disruption = disruption.loc[disruption.value > 0]
        disruption['value'] = 1- disruption['value']
        event_dict[id_] = disruption['value'].to_dict()
    
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

    comm_ind_map = { 
'acof'	 : 'Agriculture',
'cash': 'Agriculture',
'cass'	: 'Agriculture',
'cement'	: 'Processing',
'coal'	: 'Processing',
'constructi'	:'Construction',
'fertilizer'	: 'Processing',
'fishery'	: 'Agriculture',
'maiz'	: 'Agriculture',
'manufactur'	:'Manufacturing',
'meat'	: 'Agriculture',
'min_rice'	: 'Agriculture',
'max_rice'	: 'Agriculture',
'pepp'	: 'Agriculture',
'petroluem'	: 'Processing',
'rcof'	: 'Agriculture',
'rubb': 'Processing',
'steel': 'Processing',
'sugar': 'Agriculture',
'swpo'	:'Processing',
'teas'	: 'Agriculture',
'wood':'Wood and Paper' }
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

def df_com_to_ind(comm_des,min_rice=True):
    comm_des.o_region = comm_des.o_region.apply(lambda x: x.replace(' ','_').replace('-','_'))
    comm_des.d_region = comm_des.d_region.apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    comm_des = comm_des.groupby([comm_des.o_region,comm_des.d_region]).sum()
    df_od = comm_des.stack(0).reset_index()
    df_od.drop(['o_region'],inplace=True,axis=1)
    df_od.columns = ['d_region','good','value']
    
    if min_rice == True:
        df_od = df_od.loc[~(df_od.good.isin(['max_rice','min_tons','max_tons']))]
    else:
         df_od = df_od.loc[~(df_od.good.isin(['min_rice','min_tons','max_tons']))]
    
    df_od['good'] = df_od.good.apply(lambda x: map_comm_ind(x))
    df_od['good'] = df_od.good.apply(lambda x: map_ind(x))
    
    df_od_ind = df_od.groupby(['d_region','good']).sum()
    df_od_ind = df_od_ind.unstack(1)
    df_od_ind.columns = df_od_ind.columns.get_level_values(1)
    
    return df_od_ind
    
if __name__ == "__main__":

    data_path = load_config()['paths']['data']

    input_file = os.path.join(data_path,'Results','Failure_results','single_edge_failures_totals_national_road_min.csv')

    event_dict = create_disruption(input_file)
