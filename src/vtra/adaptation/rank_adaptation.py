# -*- coding: utf-8 -*-
"""Rank adaptation results
"""
import os
import sys

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from vtra.utils import *


def main():
    """Process results
    """
    config = load_config()
    data_path = config['paths']['data']
    output_path = config['paths']['output']
    figure_path = config['paths']['figures']

    province_list = ['National Road','Lao Cai','Binh Dinh','Thanh Hoa']
    base_year = 2016
    growth_scenarios = [(6.5,'forecast')] # could include low, high
    duration_max = [10] # could extend to consider other durations

    val_cols = ['hazard_type','model','climate_scenario','year','adapt_strategy','road_cond',
            'road_length','width','min_exposure_length','max_exposure_length',
            'min_daily_loss_2016','max_daily_loss_2016','min_initial_cost','max_initial_cost',
            'min_benefit_npv','max_benefit_npv','min_cost_npv','max_cost_npv','min_adapt_npv',
            'max_adapt_npv','min_bc_ratio','max_bc_ratio']
    var_cols = ['min_exposure_length','max_exposure_length',
            'min_daily_loss_2016','max_daily_loss_2016','min_initial_cost','max_initial_cost',
            'min_benefit_npv','max_benefit_npv','min_cost_npv','max_cost_npv','min_adapt_npv',
            'max_adapt_npv','min_bc_ratio','max_bc_ratio']
    comb_cols = ['Hazard exposure (km)','Loss (USD million/day)','Initial Cost (USD million)',
            'NPV benefit (USD million)','NPV Cost (USD million)','NPV adaptation (USD million)','BCR']
    divisor = [1000] +[1000000]*5 + [1]
    adapt_rank_excel = os.path.join(data_path,'Results','Adaptation_results','asset_ranks_adapt_options_2.xlsx')
    excel_writer = pd.ExcelWriter(adapt_rank_excel)
    for prn in range(len(province_list)):
        province = province_list[prn]
        # set all paths for all input files we are going to use
        if province == 'National Road':
            province_name = province.replace(' ','_').lower()
            adapt_output_excel = os.path.join(data_path,'Results','Adaptation_results','single_edge_failures_scenarios_{0}_adapt_options.xlsx'.format(province_name))
            hazard_excel = os.path.join(output_path,'hazard_scenarios','national_scale_hazard_intersections.xlsx')
            hazard_sheet = 'road'
        else:
            province_name = province.replace(' ','').lower()
            adapt_output_excel = os.path.join(data_path,'Results','Adaptation_results','single_edge_failures_commune_access_scenarios_{0}_5_tons_adapt_options.xlsx'.format(province_name))
            hazard_excel = os.path.join(output_path,'hazard_scenarios','province_roads_hazard_intersections.xlsx')
            hazard_sheet = province_name
            val_cols = ['asset_type'] + val_cols

        hazard_df = pd.read_excel(hazard_excel, sheet_name = hazard_sheet)
        hazard_df = hazard_df[['edge_id','commune_name','district_name','province_name']]
        hazard_df = hazard_df.drop_duplicates(subset=['edge_id','commune_name','district_name','province_name'], keep='first')
        hazard_df = hazard_df.set_index('edge_id')

        edge_list = []
        for eid in hazard_df.index:
            commune_names = ", ".join(list(set(hazard_df.loc[[eid],'commune_name'].values.tolist())))
            district_names = ", ".join(list(set(hazard_df.loc[[eid],'district_name'].values.tolist())))
            province_names = ", ".join(list(set(hazard_df.loc[[eid],'province_name'].values.tolist())))
            edge_list.append((eid, commune_names, district_names, province_names))

        del hazard_df

        edge_df = pd.DataFrame(edge_list, columns = ['edge_id','commune_name','district_name','province_name'])

        for g in range(len(growth_scenarios)):
            grth = growth_scenarios[g]
            for dur in range(len(duration_max)):
                adapt_df = pd.read_excel(adapt_output_excel, sheet_name = grth[1] + '_' + str(duration_max[dur]))
                adapt_df = pd.merge(adapt_df, edge_df, on = ['edge_id'], how = 'left').fillna(0)
                adapt_df = adapt_df[adapt_df['max_daily_loss_2016'] > 0]
                adapt_df = adapt_df[['edge_id','commune_name','district_name','province_name'] + val_cols]
                adapt_df = adapt_df.drop_duplicates(
                    subset=['edge_id','commune_name','district_name','province_name']+val_cols,
                    keep='first')
                adapt_df = adapt_df.sort_values(
                    'max_bc_ratio', ascending=False
                ).groupby(
                    ['edge_id','commune_name','district_name','province_name'], as_index=False
                ).first()

                adapt_df = adapt_df.sort_values(by='max_daily_loss_2016', ascending=False)
                sel_df = adapt_df.head(50)
                sel_df['road_length'] = sel_df.apply(
                    lambda x:'{:.2f}'.format(0.001*x['road_length']),
                    axis = 1)
                cols_rename = {
                    'commune_name': 'Road communes',
                    'district_name': 'Road districts',
                    'province_name': 'Road provinces',
                    'road_length': 'Length (km)',
                    'width':'Width',
                    'adapt_strategy':'Adaptation strategy'
                }
                sel_df.rename(columns=cols_rename, inplace=True)
                for c in range(len(comb_cols)):
                    sel_df[comb_cols[c]] = sel_df.apply(lambda x: '{:.2f} - {:.2f}'.format(
                        x[var_cols[2*c]]/divisor[c],
                        x[var_cols[2*c+1]]/divisor[c]
                    ), axis =1)
                    sel_df.drop([var_cols[2*c], var_cols[2*c+1]], axis=1, inplace=True)

                sel_df.to_excel(excel_writer, sheet_name = province_name + '_maxflow', index = False)
                excel_writer.save()

                print ('Done with', province)




if __name__ == "__main__":
    main()
