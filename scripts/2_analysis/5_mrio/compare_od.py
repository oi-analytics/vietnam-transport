# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 13:29:58 2018

@author: cenv0574
"""

import os
import json

import pandas as pd
import geopandas as gpd

from prepare_table import load_config,load_provincial_stats,load_output,estimate_gva,get_final_sector_classification


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

def create_zero_proxies(od_table,write_to_csv=True):
 
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
            subset['gdp'] = subset['gdp'].apply(lambda x: round(x,2))
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
    od_table = load_od(data_path)
    
    #rename goods to sectors
    mapper = map_sectors_to_od(od_table)
    
    od_table = od_table.rename(columns = mapper)
    
    od_table = od_table.groupby(od_table.columns, axis=1).sum()
    
    create_zero_proxies(od_table,write_to_csv=True)
    
#    od_sum = od_table.groupby(['Destination','Origin']).sum().sum(axis=1).unstack(1)
    
#    mrio_sum = mrio_vnm.groupby(mrio_vnm.columns.get_level_values(0),axis=1).sum().groupby(mrio_vnm.index.get_level_values(0),axis=0).sum()
    