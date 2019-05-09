"""Provincial road network risks and adaptation maps
"""
import os
import sys
from collections import OrderedDict

import ast
import numpy as np
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from vtra.utils import *


def main():
    config = load_config()

    regions = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
    hazard_cols = ['hazard_type','climate_scenario','year']
    duration = 10

    hazard_set = [
        {
            'hazard': 'landslide',
            'name': 'Landslide'
        },
        {
            'hazard': 'flashflood',
            'name':'Flashflood'

        },
        {
            'hazard': 'flooding',
            'name': 'Fluvial flooding'
        },
        {
            'hazard': 'typhoon flooding',
            'name': 'Typhoon flooding'
        }
    ]

    adapt_cols = ['min_benefit','min_ini_adap_cost','min_ini_adap_cost_perkm','min_tot_adap_cost','min_tot_adap_cost_perkm','min_bc_ratio',\
                'max_benefit','max_ini_adap_cost','max_ini_adap_cost_perkm','max_tot_adap_cost','max_tot_adap_cost_perkm','max_bc_ratio']
    for region in regions:
        region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                                   '{}_roads_edges.shp'.format(region.lower().replace(' ', '')))

        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                                   'single_edge_failures_minmax_{}_5_tons_100_percent_disrupt.csv'.format(region.lower().replace(' ', '')))

        region_file = gpd.read_file(region_file_path,encoding='utf-8')
        region_file.drop('geometry',axis=1,inplace=True)
        flow_file = pd.read_csv(flow_file_path)
        region_file = pd.merge(region_file,flow_file,how='left', on=['edge_id']).fillna(0)
        del flow_file

        asset_locations = pd.read_excel(os.path.join(config['paths']['output'],
            'hazard_scenarios','province_roads_hazard_intersections.xlsx'),sheet_name=region.replace(' ','').lower())
        asset_locations = asset_locations[['edge_id','commune_name','district_name']]
        asset_locations = asset_locations.drop_duplicates(subset=['edge_id','commune_name','district_name'], keep='first')

        region_file = pd.merge(region_file,asset_locations,how='left', on=['edge_id']).fillna(0)
        region_file = region_file[region_file['max_econ_impact'] > 0]

        region_file[['edge_id','asset_type','commune_name','district_name','min_econ_impact','max_econ_impact']].to_csv(os.path.join(config['paths']['output'],
            'network_stats',
            '{}_criticality.csv'.format(region.replace(' ','').lower())
            ), index=False
        )

        flow_file_path = os.path.join(config['paths']['output'], 'adaptation_results',
                                   'output_adaptation_{}_10_days_max_disruption_fixed_parameters.csv'.format(region.lower().replace(' ', '')))

        fail_scenarios = pd.read_csv(flow_file_path)
        fail_scenarios = pd.merge(fail_scenarios,region_file[['edge_id','commune_name','district_name']],how='left', on=['edge_id']).fillna(0)
        fail_scenarios['min_eael'] = duration*fail_scenarios['min_duration_wt']*fail_scenarios['risk_wt']*fail_scenarios['min_econ_impact']
        fail_scenarios['max_eael'] = duration*fail_scenarios['max_duration_wt']*fail_scenarios['risk_wt']*fail_scenarios['max_econ_impact']
        all_edge_fail_scenarios = fail_scenarios[hazard_cols + ['edge_id','asset_type','commune_name','district_name','min_eael','max_eael']]
        all_edge_fail_scenarios = all_edge_fail_scenarios.groupby(hazard_cols + ['edge_id','asset_type','commune_name','district_name'])['min_eael','max_eael'].max().reset_index()

        all_edge_fail_scenarios.to_csv(os.path.join(config['paths']['output'],
            'network_stats',
            '{}_eael.csv'.format(region.replace(' ','').lower())
            ), index=False
        )


        # Climate change effects
        all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(['hazard_type','edge_id','asset_type','commune_name','district_name'])
        scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
        change_tup = []
        for sc in scenarios:
            eael = all_edge_fail_scenarios.loc[[sc], 'max_eael'].values.tolist()
            yrs = all_edge_fail_scenarios.loc[[sc], 'year'].values.tolist()
            cl = all_edge_fail_scenarios.loc[[sc], 'climate_scenario'].values.tolist()
            if 2016 not in yrs:
                change_tup += list(zip([sc[0]]*len(cl),[sc[1]]*len(cl),[sc[2]]*len(cl),[sc[3]]*len(cl),[sc[4]]*len(cl),cl,yrs,[0]*len(cl),eael,[1e9]*len(cl)))
            elif len(yrs) > 1:
                vals = list(zip(cl,eael,yrs))
                vals = sorted(vals, key=lambda pair: pair[-1])
                change = 100.0*(np.array([p for (c,p,y) in vals[1:]]) - vals[0][1])/vals[0][1]
                cl = [c for (c,p,y) in vals[1:]]
                yrs = [y for (c,p,y) in vals[1:]]
                fut = [p for (c,p,y) in vals[1:]]
                change_tup += list(zip([sc[0]]*len(cl),[sc[1]]*len(cl),[sc[2]]*len(cl),[sc[3]]*len(cl),[sc[4]]*len(cl),cl,yrs,[vals[0][1]]*len(cl),fut,change))

        change_df = pd.DataFrame(change_tup,columns=['hazard_type','edge_id','asset_type','commune_name','district_name','climate_scenario','year','current','future','change(%)'])
        change_df.to_csv(os.path.join(config['paths']['output'],
            'network_stats',
            '{}_eael_climate_change.csv'.format(region.replace(' ','').lower())
            ), index=False
        )


if __name__ == '__main__':
    main()
