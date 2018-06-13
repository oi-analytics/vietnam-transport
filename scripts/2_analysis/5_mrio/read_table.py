# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 08:47:38 2018

@author: cenv0574
"""

import os

import pandas as pd
#import geopandas as gpd

from prepare_table import load_sectors,load_config,load_provincial_stats,estimate_gva

def load_output(data_path,provinces):
    
    # prepare index and cols
    region_names = list(provinces.name_eng)
    rowcol_names = list(load_sectors(data_path)['mapped'].unique())

    rows = [x for x in rowcol_names if (x.startswith('sec') | x.startswith('row'))]*len(region_names)
    cols = [x for x in rowcol_names if (x.startswith('sec') | x.startswith('col'))]*len(region_names)
    
    region_names_list = [item for sublist in [[x]*12 for x in region_names] for item in sublist]
    
    index_mi = pd.MultiIndex.from_arrays([region_names_list,rows], names=('region', 'row'))
    column_mi = pd.MultiIndex.from_arrays([region_names_list,cols], names=('region', 'col'))
    
    # read output
    output_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','output.csv')
    output_df = pd.read_csv(output_path,header=None)
    output_df.index = index_mi
    output_df.columns = column_mi
    
    row_sum = output_df.loc[[x.find('sec')==0 for x in list(output_df.index.get_level_values(1))]].sum(axis=1)
    col_sum = output_df[[x for x in list(output_df.columns) if x[1].find('sec')==0]].sum(axis=0)
    
    diff =  (row_sum - col_sum).unstack(1)
    balanced = (row_sum - col_sum).sum() < 1
    
    # create predefined index and col, which is easier to read
    sector_only = [x for x in rowcol_names if x.startswith('sec')]*len(region_names)
    row_only =  [x for x in rowcol_names if x.startswith('row')]*len(region_names) 
    col_only  =  [x for x in rowcol_names if x.startswith('col')]*len(region_names) 

    region_row = [item for sublist in [[x]*9 for x in region_names] for item in sublist] + [item for sublist in [[x]*3 for x in region_names] for item in sublist] 
    region_col = [item for sublist in [[x]*9 for x in region_names] for item in sublist] + [item for sublist in [[x]*3 for x in region_names] for item in sublist] 

#    index_mi_reorder = pd.MultiIndex.from_arrays([region_row,sector_only+row_only], names=('region', 'row'))
    column_mi_reorder = pd.MultiIndex.from_arrays([region_col,sector_only+col_only], names=('region', 'col'))
 
    #sum va and imports
    tax_sub = output_df.loc[output_df.index.get_level_values(1)=='row1'].sum(axis='index')
    import_ = output_df.loc[output_df.index.get_level_values(1)=='row2'].sum(axis='index')
    valueA = output_df.loc[output_df.index.get_level_values(1)=='row3'].sum(axis='index')
    
    output_new = pd.concat([output_df.loc[~output_df.index.get_level_values(1).isin(['row1','row2','row3'])],pd.DataFrame(tax_sub).T,
                             pd.DataFrame(import_).T,pd.DataFrame(valueA).T])
    
#    output_new = output_new.reindex(index_mi_reorder,axis='index')
    output_new = output_new.reindex(column_mi_reorder,axis='columns')

    # write to new csv    
    output_path_new = os.path.join(data_path,'IO_analysis','MRIO_TABLE','output_reordered.csv')
    output_new.to_csv(output_path_new)

    return output_new

if __name__ == "__main__":

    # load path
    data_path = load_config()['paths']['data']

    # load provincial shapefile
    provinces = load_provincial_stats(data_path)
    provinces.name_eng = provinces.name_eng.apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    # get list of names and concat
    list_names = list(provinces.name_eng)
    list_string = " , ".join(list_names)
    
    # estimate gross value added
    provinces['raw_gva'] = estimate_gva(provinces,in_million=True)
    
    # get reordered mrio with new region classification
    mrio_vnm = load_output(data_path,provinces)
    
    # get sums per region 
    Xcol = mrio_vnm.sum(axis='columns')
    Xrow = mrio_vnm.sum(axis='index')
    
#    
#    X = X.unstack(1)
#    X = X[['secA','secB','secC','secD','secE','secF','secG','secH','secI']]
#    X = X.T    
    
    
    
