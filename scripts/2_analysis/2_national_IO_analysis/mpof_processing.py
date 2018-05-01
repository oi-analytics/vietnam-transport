# -*- coding: utf-8 -*-
"""
Created on Tue May  1 08:15:32 2018

@author: cenv0574
"""

import os
import json

import pandas as pd

#def main():

def map_comm_ind(x):

    comm_ind_map = { 'sugar' : 'Agriculture',
    'wood' : 'Agriculture',
    'steel' : 'Metal Products',
    'constructi' : 'Construction',
    'cement' : 'Mining and Quarrying',
    'fertilizer' : 'Petroleum, Chemical and Non-Metallic Mineral Products',
    'coal' : 'Mining and Quarrying',
    'petroluem' : 'Petroleum, Chemical and Non-Metallic Mineral Products',
    'manufactur' : 'Electrical and Machinery',
    'fishery' : 'Fishing',
    'meat' : 'Food & Beverages',
    'rice' : 'Agriculture',
    'cash' : 'Agriculture',
    'cass' : 'Agriculture',
    'teas' : 'Food & Beverages',
    'maiz' : 'Agriculture',
    'rubb' : 'Other Manufacturing',
    'swpo' : 'Agriculture',
    'acof' : 'Food & Beverages',
    'rcof' : 'Food & Beverages',
    'pepp' : 'Agriculture' }

    return comm_ind_map[x]

def map_ind(x):
    
    ind_map = {'i1':'Agriculture',
           'i2':'Fishing',
           'i3':'Mining and Quarrying',
           'i4':'Food & Beverages',
           'i5':'Textiles and Wearing Apparel',
           'i6':'Wood and Paper',
           'i7':'Petroleum, Chemical and Non-Metallic Mineral Products',
           'i8':'Metal Products',
           'i9':'Electrical and Machinery',
           'i10':'Transport Equipment',
           'i11':'Other Manufacturing',
           'i12':'Recycling',
           'i13':'Electricity, Gas and Water',
           'i14':'Construction',
           'i15':'Maintenance and Repair',
           'i16':'Wholesale Trade',
           'i17':'Retail Trade',
           'i18':'Hotels and Restraurants',
           'i19':'Transport',
           'i20':'Post and Telecommunications',
           'i21':'Finacial Intermediation and Business Activities','i22':'Public Administration','i23':'Education, Health and Other Services','i24':'Private Households','i25':'Others','i26':'Re-export & Re-import'}

    return  {v: k for k, v in ind_map.items()}[x]


if __name__ == "__main__":

    # Define current directory and data directory
    config_path = os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', '..', '..', 'config.json')
    )
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    data_path = config['paths']['data']
    output_path = config['paths']['output']
    figure_path = config['paths']['figures']

    # read commodity descriptions
    comm_des = pd.read_excel(os.path.join(data_path,'Results','commodity_descriptions_and_totals.xlsx'),sheet_name='Sheet1')
    
    # read vietnam mpof domestic results
    vnm_mpof_dom = pd.read_csv(os.path.join(data_path,'Results','vnm_road_rail_edge_multi_failure.csv'))
    
    # get industry daily tonnage totals
    ind_des = comm_des.copy()
    ind_des['industry'] = ind_des.commodity_field.apply(map_comm_ind)
    ind_des['ind'] = ind_des.industry.apply(map_ind)
    ind_des = ind_des.drop(['commodity_name','commodity_field'],axis=1)
    
    # convert commodities to industries for EORA impact modelling
    vnm_mpof_dom['industry'] = vnm_mpof_dom.industry.apply(map_comm_ind)
    vnm_mpof_dom['ind'] = vnm_mpof_dom.industry.apply(map_ind)
    
    # create scenarios per region
    regional_scenarios = vnm_mpof_dom.groupby(['region_flooded','typhoon_level']).sum()
    
    # create typhoon scenarios
    typhoon_scenarios = vnm_mpof_dom.groupby(['typhoon_level']).sum()
    
    # get scenarios and industries for typhoons total
    typhoon_unique = vnm_mpof_dom.groupby(['typhoon_level','ind']).sum()
    
    # get scenarios and industries for typhoons and regions
    typhoon_region_unique = vnm_mpof_dom.groupby(['region_flooded','typhoon_level','ind']).sum()
    
    # =============================================================================
    #     # create disruption events
    # =============================================================================
    
    typhoon_events = {}
    
    
    