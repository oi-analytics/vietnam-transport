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

def create_regional_proxy(data_path,regions,write_to_csv=True):
    
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

def create_sector_proxies(data_path,regions,write_to_csv=True):
    
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

def load_od(data_path):
    
    od_path = os.path.join(data_path,'OD_data','OD_Transport_Data_Set.xlsx')

    od_table = pd.read_excel(od_path,sheet_name='64 x 64_2008_13 goods')
    
    od_table = od_table.drop(['O','D'],axis='columns')
    od_table = od_table.dropna(subset=['Name O','name D'],axis='index')
    od_table['Name O'] = od_table['Name O'].apply(lambda x: x.replace(' ','_').replace('-','_'))
    od_table['name D'] = od_table['name D'].apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    od_table = od_table.rename(columns = {'Name O' : 'Origin','name D':'Destination'})
    
    return od_table

def map_sectors_to_od(od_table):

    goods = [x for x in od_table.columns if not (x.startswith('Origin') | x.startswith('Destination'))]

    sectors_conn = ['secA','secC','secE','secF','secG','secF','secF','secB','secF','secF','secF','secA','secC']
    
    return dict(zip(goods,sectors_conn))

def map_regions():
    
    return {'An_Giang':'An_Giang',
 'Ba_Ria_Vung_Tau': 'Ba_Ria_Vung_Tau',
 'Bac_Giang': 'Bac_Giang',
 'Bac_Kan': 'Bac_Kan',
 'Bac_Lieu': 'Bac_Lieu',
 'Bac_Ninh': 'Bac_Ninh',
 'Ben_Tre': 'Ben_Tre',
 'Binh_Dinh': 'Binh_Dinh',
 'Binh_Duong': 'Binh_Duong',
 'Binh_Phuoc': 'Binh_Phuoc',
 'Binh_Thuan': 'Binh_Thuan',
 'Ca_Mau': 'Ca_Mau',
 'Can_Tho': 'Can_Tho',
 'Cao_Bang': 'Cao_Bang',
 'Da_Nang': 'Da_Nang',
 'Dak_Lak': 'Dak_Lak',
 'Dak_Nong': 'Dak_Nong',
 'Dien_Bien_Phu': 'Dien_Bien',
 'Dong_Nai': 'Dong_Nai',
 'Dong_Thap': 'Dong_Thap',
 'Gia_Lai': 'Gia_Lai',
 'Ha_Giang': 'Ha_Giang',
 'Ha_Nam': 'Ha_Nam',
 'Ha_Tay': 'Ha_Noi',
 'Ha_Noi': 'Ha_Noi',
 'Ha_Tinh': 'Ha_Tinh',
 'Hai_Duong': 'Hai_Duong',
 'Hai_Phong': 'Hai_Phong',
 'Hau_Giang': 'Hau_Giang',
 'HCM': 'Ho_Chi_Minh',
 'Hoa_Binh': 'Hoa_Binh',
 'Hung_Yen': 'Hung_Yen',
 'Khanh_Hoa': 'Khanh_Hoa',
 'Kien_Giang': 'Kien_Giang',
 'Kon_Tum': 'Kon_Tum',
 'Lai_Chau': 'Lai_Chau',
 'Lam_Dong': 'Lam_Dong',
 'Lang_Son': 'Lang_Son',
 'Lao_Cai': 'Lao_Cai',
 'Long_An': 'Long_An',
 'Nam_Dinh': 'Nam_Dinh',
 'Nghe_An': 'Nghe_An',
 'Ninh_Binh': 'Ninh_Binh',
 'Ninh_Thuan': 'Ninh_Thuan',
 'Phu_Tho': 'Phu_Tho',
 'Phu_Yen': 'Phu_Yen',
 'Quang_Binh': 'Quang_Binh',
 'Quang_Nam': 'Quang_Nam',
 'Quang_Ngai': 'Quang_Ngai',
 'Quang_Ninh': 'Quang_Ninh',
 'Quang_Tri': 'Quang_Tri',
 'Soc_Trang': 'Soc_Trang',
 'Son_La': 'Son_La',
 'Tay_Ninh': 'Tay_Ninh',
 'Thai_Binh': 'Thai_Binh',
 'Thai_Nguyen': 'Thai_Nguyen',
 'Thanh_Hoa': 'Thanh_Hoa',
 'Thua_Thien_Hue': 'Thua_Thien_Hue',
 'Tien_Giang': 'Tien_Giang',
 'Tra_Vinh': 'Tra_Vinh',
 'Tuyen_Quang': 'Tuyen_Quang',
 'Vinh_Long': 'Vinh_Long',
 'Vinh_Phuc': 'Vinh_Phuc',
 'Yen_Bai': 'Yen_Bai'}

def create_zero_proxies(data_path,od_table,write_to_csv=True):
 
    # get sector list
    sector_list = get_final_sector_classification()+['other1','other2','other3']
    sector_list = [x+str(1) for x in sector_list]

    #map sectors to be the same
    mapper = map_regions()
    od_table['Destination'] = od_table['Destination'].apply(lambda x: mapper[x])
    od_table['Origin'] = od_table['Origin'].apply(lambda x: mapper[x])
    
    od_table = od_table.loc[od_table['Destination'] != od_table['Origin']]

    od_sum = pd.DataFrame(od_table.groupby(['Destination','Origin']).sum().sum(axis=1))
    od_sum.reset_index(inplace=True)
    od_sum.columns = ['Destination','Origin','gdp']
    
    for sector in sector_list:
        if sector in ['other1','other2','other3']:
            subset = od_sum.copy()
            subset['year'] = 2010
            subset['sector'] = sector
            subset['gdp'] = 0
            combine = []
            for sector2 in sector_list:
                sub_subset = subset.copy()
                sub_subset['subsector'] = sector2
                combine.append(sub_subset)
        else:
            subset = od_sum.copy()
            subset = subset.loc[od_sum.gdp == 0]
            subset['year'] = 2010
            subset['sector'] = sector
            subset['gdp'] = 0 #subset['gdp'].apply(lambda x: round(x,2))
            combine = []
            for sector2 in sector_list:
                sub_subset = subset.copy()
                sub_subset['subsector'] = sector2
                combine.append(sub_subset)
    
        all_ = pd.concat(combine)
        final_sub = all_[['year','sector','Origin','subsector','Destination','gdp']]
        final_sub.columns = ['year','sector','region','sector','region','gdp']
        
        if write_to_csv == True:
            csv_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','proxy_trade_{}.csv'.format(sector[:-1]))
            final_sub.to_csv(csv_path,index=False)

def create_proxies(data_path):
    
    provinces = load_provincial_stats(data_path)
    od_table = load_od(data_path)
    
    create_indices(data_path,provinces,write_to_csv=True)
    create_regional_proxy(data_path,provinces,write_to_csv=True)
    create_sector_proxies(data_path,provinces,write_to_csv=True)
    create_zero_proxies(data_path,od_table,write_to_csv=True)
    
    
    
