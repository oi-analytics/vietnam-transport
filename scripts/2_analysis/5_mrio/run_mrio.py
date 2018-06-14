# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 15:57:42 2018

@author: cenv0574
"""

import os
import pandas as pd

import subprocess

from prepare_table import load_config,load_sectors,load_table,load_provincial_stats,estimate_gva
from read_table import load_output
from ras_method import ras_method


if __name__ == "__main__":

    data_path = load_config()['paths']['data']
   
    vnm_IO = load_table(data_path)
    vnm_IO_rowcol = load_sectors(data_path)

    # load provincial shapefile
    provinces = load_provincial_stats(data_path)
    provinces.name_eng = provinces.name_eng.apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    # get list of names and concat
    list_names = list(provinces.name_eng)
    list_string = " , ".join(list_names)

    # estimate gross value added
    provinces['raw_gva'] = estimate_gva(provinces,in_million=True)

    # run mrio_disaggregate
    mrio_tool_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','mrio_disaggregate')
    settings_tool_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','settings_trade.yml')

    p = subprocess.Popen(['mrio_disaggregate','settings_trade.yml'], cwd=os.path.join(data_path,'IO_analysis','MRIO_TABLE'))
    p.wait()
    # get reordered mrio with new region classification
    Xin = load_output(data_path,provinces)

    X0 = Xin.as_matrix() #[:9]
   
    u = X0.sum(axis=1)
    v = X0.sum(axis=0)
    v[:(len(u)-3)] = u[:-3]
    
    X1 = ras_method(X0,u,v,eps=1e-6)
#    
    u_new = X1.sum(axis=1)
    v_new = X1.sum(axis=0)
    
    Xin.iloc[:,:] = X1
    
    TotIn = Xin.sum(axis='columns')
    TotOut = Xin.sum(axis='index')
    
    Xout = Xin.iloc[:-3,:]
    Xout.index = pd.MultiIndex.from_tuples(Xout.index, names=['region', 'sector'])
    Region_sum = Xout.groupby(Xout.columns.get_level_values(0),axis='columns').sum().groupby(Xout.index.get_level_values(0),axis='index').sum()
    
    Region_sum.to_csv(os.path.join(data_path,'IO_analysis','MRIO_TABLE','region_trade.csv'))
    