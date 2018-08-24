# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 08:47:38 2018

@author: cenv0574
"""

import os
import json

import pandas as pd
import geopandas as gpd

from itertools import product

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
        return list(((regions.pro_nfirm*regions.laborcost)+(regions.pro_nfirm*regions.capital)))

def create_proxies(data_path,notrade=False,own_production_ratio=0.9):
    
    provinces = load_provincial_stats(data_path)
    provinces.name_eng = provinces.name_eng.apply(lambda x: x.replace(' ','_').replace('-','_'))
    od_table = load_od(data_path)
    
    create_indices(data_path,provinces,write_to_csv=True)
    create_regional_proxy(data_path,provinces,write_to_csv=True)
    create_sector_proxies(data_path,provinces,write_to_csv=True)
    create_zero_proxies(data_path,od_table,notrade=notrade,write_to_csv=True)
    if notrade == False:
        create_level14_proxies(data_path,od_table,own_production_ratio,write_to_csv=True)

def create_regional_proxy(data_path,regions,write_to_csv=True):
    
    regions['raw_gva'] = estimate_gva(regions) #regions['pro_nfirm']*regions['laborcost'] + regions['pro_nfirm']*regions['capital']
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

def get_trade_value(x,sum_use,sector,own_production_ratio=0.9):
    if x.Destination == x.Origin:
        try:
            return list(sum_use.loc[(sum_use['region'] == x.Destination) & (sum_use['sector'] == sector)]['value'])[0]*own_production_ratio
        except:
            return 1
    elif x.gdp == 0:
        return 0
    else:
        try:
            return list(sum_use.loc[(sum_use['region'] == x.Destination) & (sum_use['sector'] == sector)]['value'])[0]*(1-own_production_ratio)*x.ratio
        except:
            return 0
        
def create_level14_proxies(data_path,od_table,own_production_ratio=0.9,write_to_csv=True):

    # get sector list
    sector_list_ini = get_final_sector_classification()+['other1','other2','other3']
    sector_list = [x+str(1) for x in sector_list_ini]

    #map sectors to be the same
#    mapper = map_regions()
#    od_table['Destination'] = od_table['Destination'].apply(lambda x: mapper[x])
#    od_table['Origin'] = od_table['Origin'].apply(lambda x: mapper[x])
     
    od_table.loc[od_table['Destination'] == od_table['Origin'],'gdp'] = 10

    od_sum = pd.DataFrame(od_table.groupby(['Destination','Origin']).sum().sum(axis=1))

    od_sum['ratio'] = od_sum.groupby(level=0).apply(lambda x:
                                                     x / float(x.sum()))
    od_sum.reset_index(inplace=True)
    od_sum.columns = ['Destination','Origin','gdp','ratio']
    
    df_pretable = pd.read_csv(os.path.join(data_path,'IO_analysis','MRIO_TABLE','notrade_trade.csv'),index_col=[0,1],header=[0,1])
    df_pretable = df_pretable.iloc[:,:567]
    sum_use = df_pretable.sum(axis=1)
    sum_use = pd.DataFrame(sum_use*0.1)
    sum_use.reset_index(inplace=True)
    sum_use.columns = ['region','sector','value']
    
    combine = []

    for sector in sector_list:
        if sector[:-1] in ['other1','other2','other3']:
            subset = od_sum.copy()
            subset['year'] = 2010
            subset['sector'] = sector
            subset['gdp'] = 0
            subset.drop('ratio',axis=1,inplace=True)
            combine.append(subset)
        else:
            subset = od_sum.copy()
            subset = subset.loc[od_sum.gdp != 0]
            subset['year'] = 2010
            subset['sector'] = sector
            subset['gdp'] =  subset.apply(lambda x: get_trade_value(x,sum_use,sector[:-1],own_production_ratio),axis=1) #subset['gdp'].apply(lambda x: round(x,2))
            subset.drop('ratio',axis=1,inplace=True)
            combine.append(subset)

    
        all_ = pd.concat(combine)
        final_sub = all_[['year','sector','Origin','Destination','gdp']]
        final_sub.columns = ['year','sector','region','region','gdp']
        
        if write_to_csv == True:
            csv_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','proxy_trade14_{}.csv'.format(sector[:-1]))
            final_sub.to_csv(csv_path,index=False)

def create_zero_proxies(data_path,od_table,notrade=False,write_to_csv=True):
 
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
    
    if notrade == True:
        od_sum['gdp'] = 0
    
    for sector in sector_list:
        if sector[:-1] in ['other1','other2','other3']:
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
            if notrade == False:
                subset = subset.loc[od_sum.gdp == 0]
            subset['year'] = 2010
            subset['sector'] = sector
            subset['gdp'] = 0 
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

def load_output(data_path,provinces,notrade=True):
    
    # prepare index and cols
    region_names = list(provinces.name_eng)
    rowcol_names = list(load_sectors(data_path)['mapped'].unique())

    rows = [x for x in rowcol_names if (x.startswith('sec') | x.startswith('row'))]*len(region_names)
    cols = [x for x in rowcol_names if (x.startswith('sec') | x.startswith('col'))]*len(region_names)
    
    region_names_list = [item for sublist in [[x]*12 for x in region_names] for item in sublist]
    
    index_mi = pd.MultiIndex.from_arrays([region_names_list,rows], names=('region', 'row'))
    column_mi = pd.MultiIndex.from_arrays([region_names_list,cols], names=('region', 'col'))
    
    # read output
    if notrade == True:
        output_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','output_notrade.csv')
    else:
         output_path = os.path.join(data_path,'IO_analysis','MRIO_TABLE','output.csv')
       
    output_df = pd.read_csv(output_path,header=None)
    output_df.index = index_mi
    output_df.columns = column_mi
    
    # create predefined index and col, which is easier to read
    sector_only = [x for x in rowcol_names if x.startswith('sec')]*len(region_names)
#    row_only =  [x for x in rowcol_names if x.startswith('row')]*len(region_names) 
    col_only  =  [x for x in rowcol_names if x.startswith('col')]*len(region_names) 

#    region_row = [item for sublist in [[x]*9 for x in region_names] for item in sublist] + [item for sublist in [[x]*3 for x in region_names] for item in sublist] 
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
    
    od_path = os.path.join(data_path,'OD_data','national_scale_od_matrix.xlsx')
    od_table = pd.read_excel(od_path,sheet_name='total')
    od_table.drop(['max_rice','min_tons','max_tons'],inplace=True,axis=1)
    
    od_table = od_table.dropna(subset=['o_region','d_region'],axis='index')
    od_table['o_region'] = od_table['d_region'].apply(lambda x: x.replace(' ','_').replace('-','_'))
    od_table['d_region'] = od_table['d_region'].apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    od_table = od_table.rename(columns = {'o_region' : 'Origin','d_region':'Destination'})
    
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
 'Dien_Bien': 'Dien_Bien',
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
 'Ho_Chi_Minh': 'Ho_Chi_Minh',
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