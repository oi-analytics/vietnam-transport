# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 13:52:00 2017

@author: cenv0574
"""

# Import
import pandas as pd 
import geopandas as gp
import numpy as np
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config
from table import io_basic
from model import MRIA_IO as MRIA
from ratmarg import ratmarg_IO


if __name__ == '__main__':

    data_path = load_config()['paths']['data']
   
    ''' Specify file path '''
    filepath =  os.path.join(data_path,'input_data','IO_VIETNAM.xlsx')

    '''Specify which countries should be included in the subset'''
  
    '''Create data input'''
    DATA = io_basic('Vietnam',filepath,2010)
    DATA.prep_data()
    
        
    '''Create model '''    
    MRIA_model = MRIA(DATA.name, DATA.countries,DATA.sectors,DATA.FD_cat)
    MRIA_model.create_sets(FD_SET=['FinDem'])
    MRIA_model.create_alias()
#
#
    '''Run model and create some output'''
    output = pd.DataFrame()
# 
    '''Specify disruption'''
    disr_dict_sup = {} #{(k,r): v for r, kv in disr.iterrows() for k,v in kv.to_dict().items()}
    disr_dict_fd = {} #{(disrupted_org[0],k,r): v for r, kv in disr.iterrows() for k,v in kv.to_dict().items()}

    '''Create model'''
    MRIA_RUN = MRIA(DATA.name,DATA.countries,DATA.sectors,EORA=False)
    
    '''Define sets and alias'''
    # CREATE SETS
    MRIA_RUN.create_sets()
    
    # CREATE ALIAS
    MRIA_RUN.create_alias()
    
    ''' Define tables and parameters'''
    Regmaxcap = 0.98
    
    ratmarg_IO(DATA,EORA=False)
#    
#    MRIA_RUN.baseline_data(DATA,disr_dict_sup,disr_dict_fd)
#    MRIA_RUN.impact_data(DATA,disr_dict_sup,disr_dict_fd)
#
#    output['x_in'] = pd.Series(MRIA_RUN.X.get_values())
#   
#    MRIA_RUN.run_impactmodel()
#  
#    output['x_out'] = pd.Series(MRIA_RUN.X.get_values())
#    output['loss'] = output['x_out'] - output['x_in']
