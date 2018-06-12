# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 13:29:58 2018

@author: cenv0574
"""

import os
import json

import pandas as pd
import geopandas as gpd

from prepare_table import load_config,load_provincial_stats,load_output,estimate_gva


def load_od(data_path):
    
    od_path = os.path.join(data_path,'OD_data','OD_Transport_Data_Set.xlsx')

    od_table = pd.read_excel(od_path,sheet_name='64 x 64_2008_13 goods')
    
    od_table = od_table.drop(['O','D'],axis='columns')
    od_table = od_table.dropna(subset=['Name O','name D'],axis='index')
    od_table['Name O'] = od_table['Name O'].apply(lambda x: x.replace(' ','_').replace('-','_'))
    od_table['name D'] = od_table['name D'].apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    od_table = od_table.rename(columns = {'Name O' : 'Origin','name D':'Destination'})
    
    return od_table

if __name__ == "__main__":
    
    # load path
    data_path = load_config()['paths']['data']

    # load provincial shapefile
    provinces = load_provincial_stats(data_path)
    provinces.name_eng = provinces.name_eng.apply(lambda x: x.replace(' ','_').replace('-','_'))
   
    # estimate gross value added
    provinces['raw_gva'] = estimate_gva(provinces,in_million=True)
    
    # get reordered mrio with new region classification
    mrio_vnm = load_output(data_path,provinces)

    # load od matrix    
    od_vnm = load_od(data_path)
