# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 08:47:38 2018

@author: cenv0574
"""

import os
import json

import pandas as pd
import geopandas as gpd

def load_config():
    # Define current directory and data directory
    config_path = os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', '..', '..', 'config.json')
    )
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    
    return config

def load_table(data_path):
    vnm_IO_path = os.path.join(data_path,"INPUT-OUTPUT TABLE 2012","IO Table 2012 English.xlsx")
    return pd.read_excel(vnm_IO_path,sheet_name='IO_clean',index_col=0)
    
def load_sectors(data_path):
    vnm_IO_path = os.path.join(data_path,"INPUT-OUTPUT TABLE 2012","IO Table 2012 English.xlsx")
    vnmIO_rowcol = pd.read_excel(vnm_IO_path,sheet_name='SectorName')

    return vnmIO_rowcol    

def map_sectors(vnm_IO_rowcol):
    
    row_only = vnm_IO_rowcol[vnm_IO_rowcol['mapped'].str.contains("row") | vnm_IO_rowcol['mapped'].str.contains("sec") ]
    col_only = vnm_IO_rowcol[vnm_IO_rowcol['mapped'].str.contains("col") | vnm_IO_rowcol['mapped'].str.contains("sec") ]
    
    return dict(zip(row_only.code,row_only.mapped)),dict(zip(col_only.code,col_only.mapped))

def aggregate_table(vnm_IO,vnm_IO_rowcol,in_million=True):
    
    #aggregate table
    mapper_row,mapper_col = map_sectors(vnm_IO_rowcol)
    vnm_IO.index = vnm_IO.index.map(mapper_row.get)
    vnm_IO.columns = vnm_IO.columns.to_series().map(mapper_col)
    
    aggregated =  vnm_IO.groupby(vnm_IO.index,axis=0).sum().groupby(vnm_IO.columns, axis=1).sum()
    
    if in_million == True:
        return aggregated/1000000
    else:
        return aggregated

def is_balanced(io_table):
        
    row = io_table.sum(axis=0)
    col = io_table.sum(axis=1)
    
    if ((row-col).sum() < 1):
        print('Table is balanced')

def load_provincial_stats(data_path):
    
    prov_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','province_level_stats.shp')
    
    return gpd.read_file(prov_path)

def estimate_gva(regions,in_million=True):
 
    if in_million == True:
        return list(((regions.pro_nfirm*regions.laborcost)+(regions.pro_nfirm*regions.capital))/1000000)
    else:
        return list(((regions.pro_nfirm*regions.laborcost)+(regions.pro_nfirm*regions.capital))/1000000)

def create_regional_proxy(regions,write_to_csv=True):
    
    subset = regions.loc[:,['name_eng','raw_gva']]
    subset['year'] = 2010
    subset['raw_gva'] = subset.raw_gva.apply(int)
    subset = subset[['year','name_eng','raw_gva']]
    subset.columns = ['year','id','gdp']
    
    if write_to_csv == True:
        csv_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','proxy_reg_vnm.csv')
        subset.to_csv(csv_path,index=False)

def create_sector_proxies(regions,write_to_csv=True):
    
    #list of sectors
    sector_list = ['secA','secB','secC','secD','secE','secF','secG','secH','secI']
    
    #get own sector classification for region file
    map_dict = map_sect_vnm_to_eng()
    regions=regions.rename(columns = map_dict)
    
    # get sectoral gva based on proportion of firms in the region
    sector_shares = regions[sector_list].multiply(regions['raw_gva'],axis='index')
    sector_shares.index = regions.name_eng
    
    for sector in sector_list:

        subset = pd.DataFrame(sector_shares.loc[:,sector])
        subset.reset_index(inplace=True,drop=False)
        subset['year'] = 2010
        subset['sector'] = sector+str(1)
        subset[sector] = subset[sector].apply(lambda x: round(x,2))
        subset = subset[['year','sector','name_eng',sector]]
        subset.columns = ['year','sector','region','gdp']
        
        if write_to_csv == True:
            csv_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','proxy_{}.csv'.format(sector))
            subset.to_csv(csv_path,index=False)
        
    

def map_sect_vnm_to_eng():
    
    map_dict = { 'nongnghiep' : 'secA',
                'khaikhoang': 'secB',
                'chebien': 'secC',
                'detmay': 'secD',
                'gogiay': 'secE',
                'sanxuat': 'secF',
                'xaydung': 'secG',
                'thuongmai': 'secH',
                'dichvu': 'secI' }

    return map_dict

def load_output(data_path,provinces):
    
    # prepare index and cols
    region_names = list(provinces.name_eng)
    rowcol_names = list(vnm_IO_rowcol['mapped'].unique())

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
    
    balanced = (row_sum - col_sum).sum() < 1
    
    # create predefined index and col, which is easier to read
    sector_only = [x for x in rowcol_names if x.startswith('sec')]*len(region_names)
    row_only =  [x for x in rowcol_names if x.startswith('row')]*len(region_names) 
    col_only  =  [x for x in rowcol_names if x.startswith('col')]*len(region_names) 

    region_row = [item for sublist in [[x]*9 for x in region_names] for item in sublist] + [item for sublist in [[x]*3 for x in region_names] for item in sublist] 
    region_col = [item for sublist in [[x]*9 for x in region_names] for item in sublist] + [item for sublist in [[x]*3 for x in region_names] for item in sublist] 

    index_mi_reorder = pd.MultiIndex.from_arrays([region_row,sector_only+row_only], names=('region', 'row'))
    column_mi_reorder = pd.MultiIndex.from_arrays([region_col,sector_only+col_only], names=('region', 'col'))
 
    output_df = output_df.reindex(index_mi_reorder,axis='index')
    output_df = output_df.reindex(column_mi_reorder,axis='columns')

    # write to new csv    
    output_path_new = os.path.join(data_path,'IO_analysis','MRIO_TABLE','output_reordered.csv')
    output_df.to_csv(output_path_new)

    return output_df

if __name__ == "__main__":

    # load path
    data_path = load_config()['paths']['data']

    # load IO table
    vnm_IO = load_table(data_path)

    # load row and col table
    vnm_IO_rowcol = load_sectors(data_path)
    
    # aggregate table
    vnm_IO_agg = aggregate_table(vnm_IO,vnm_IO_rowcol)
    
    # load provincial shapefile
    provinces = load_provincial_stats(data_path)
    provinces.name_eng = provinces.name_eng.apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    # get list of names and concat
    list_names = list(provinces.name_eng)
    list_string = " , ".join(list_names)
    
    # estimate gross value added
    provinces['raw_gva'] = estimate_gva(provinces,in_million=True)
    
    # create regional proxy file
    create_regional_proxy(provinces,write_to_csv=True)
    
    #create sectoral proxies
    create_sector_proxies(provinces,write_to_csv=True)

    # get reordered mrio with new region classification
    mrio_vnm = load_output(data_path,provinces)
    
    # get sums per region 
    X = mrio_vnm.sum(axis='columns')
    X = X.unstack(1)
    X = X[['secA','secB','secC','secD','secE','secF','secG','secH','secI']]
    X = X.T    
    
    
    
