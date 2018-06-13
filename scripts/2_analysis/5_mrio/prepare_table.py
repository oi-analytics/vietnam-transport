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

def get_final_sector_classification():
    return ['secA','secB','secC','secD','secE','secF','secG','secH','secI']

def map_sectors(vnm_IO_rowcol):
    
    row_only = vnm_IO_rowcol[vnm_IO_rowcol['mapped'].str.contains("row") | vnm_IO_rowcol['mapped'].str.contains("sec") ]
    col_only = vnm_IO_rowcol[vnm_IO_rowcol['mapped'].str.contains("col") | vnm_IO_rowcol['mapped'].str.contains("sec") ]
    
    return dict(zip(row_only.code,row_only.mapped)),dict(zip(col_only.code,col_only.mapped))

def aggregate_table(vnm_IO,vnm_IO_rowcol,in_million=True):
    
    sectors = get_final_sector_classification()
    
    #aggregate table
    mapper_row,mapper_col = map_sectors(vnm_IO_rowcol)
    vnm_IO.index = vnm_IO.index.map(mapper_row.get)
    vnm_IO.columns = vnm_IO.columns.to_series().map(mapper_col)
    
    aggregated =  vnm_IO.groupby(vnm_IO.index,axis=0).sum().groupby(vnm_IO.columns, axis=1).sum()

    aggregated = aggregated.reindex(sectors+['col1','col2','col3'],axis='columns')
    aggregated = aggregated.reindex(sectors+['row1','row2','row3'],axis='index')

    
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
    subset['raw_gva'] = subset.raw_gva.apply(int)/(subset['raw_gva'].sum(axis='index'))
    subset = subset[['year','name_eng','raw_gva']]
    subset.columns = ['year','id','gdp']
    
    if write_to_csv == True:
        csv_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','proxy_reg_vnm.csv')
        subset.to_csv(csv_path,index=False)

def create_indices(data_path,provinces,write_to_csv=True):

    # prepare index and cols
    region_names = list(provinces.name_eng)
    rowcol_names = list(load_sectors(data_path)['mapped'].unique())

    rows = [x for x in rowcol_names if (x.startswith('sec') | x.startswith('row'))]*len(region_names)
    
    region_names_list = [item for sublist in [[x]*12 for x in region_names] for item in sublist]
    
    indices = pd.DataFrame([region_names_list,rows]).T
    indices.columns = ['region','sector']
    indices['sector'] = indices['sector'].apply(lambda x: x.replace('row','other'))

    if write_to_csv == True:
        csv_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','indices_mrio.csv')
        indices.to_csv(csv_path,index=False)

def create_sector_proxies(regions,write_to_csv=True):
    
    #list of sectors
    sector_list = get_final_sector_classification()
    
    #get own sector classification for region file
    map_dict = map_sect_vnm_to_eng()
    regions=regions.rename(columns = map_dict)
    
    # get sectoral gva based on proportion of firms in the region
    sector_shares = regions[sector_list].multiply(regions['raw_gva'],axis='index')
    sector_shares.index = regions.name_eng
    
    for sector in sector_list+['other1','other2','other3']:
        if sector in ['other1','other2','other3']:
            subset = pd.DataFrame(sector_shares.sum(axis='columns')).divide(pd.DataFrame(sector_shares.sum(axis='columns')).sum(axis='index'))
            subset.columns  = [sector]
        else:
            subset = pd.DataFrame(sector_shares.loc[:,sector]).divide(pd.DataFrame(sector_shares.loc[:,sector]).sum(axis='index'))
        subset.reset_index(inplace=True,drop=False)
        subset['year'] = 2010
        subset['sector'] = sector+str(1)
        subset[sector] = subset[sector].apply(lambda x: round(x,7))
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


    
    
