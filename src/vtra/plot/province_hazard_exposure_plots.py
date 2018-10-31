"""Road network flows
"""
import os
import sys
from collections import OrderedDict

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
    plot_set = [
        {
            'hazard': 'landslide',
            'color': ['#fdd0a2', '#fdae6b', '#fd8d3c', '#e6550d', '#a63603','#d9d9d9'],
            'name': 'Landslide'
        },
        {
            'hazard': 'flashflood',
            'color': ['#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f','#d9d9d9'],
            'name':'Flashflood'

        },
        {
            'hazard': 'flooding',
            'color': ['#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c','#d9d9d9'],
            'name': 'Fluvial flooding'
        },
        {
            'hazard': 'typhoon flooding',
            'color': ['#c7e9c0', '#a1d99b', '#74c476', '#31a354', '#006d2c','#d9d9d9'],
            'name': 'Typhoon flooding'
        }
    ]
    province_pth = os.path.join(config['paths']['output'], 
            'network_stats',
            'province_roads_hazards_stats.xlsx')
    commune_shp = os.path.join(config['paths']['data'], 'Vietnam_boundaries',
                                'who_boundaries', 'who_communes.shp')

    labels = ['0 to 10', '10 to 20', '20 to 30', '30 to 40', '40 to 100', 'No value']
    change_colors = ['#1a9850','#66bd63','#a6d96a','#d9ef8b','#fee08b','#fdae61','#f46d43','#d73027','#d9d9d9']
    change_labels = ['< -40','-40 to -20','-20 to -10','-10 to 0','0 to 10','10 to 20','20 to 40','> 40','No change/value']
    change_ranges = [(-1e10,-40),(-40,-20),(-20,-10),(-10,0),(0.001,10),(10,20),(20,40),(40,1e10)]

    for region in regions:
        region_file = gpd.read_file(commune_shp,encoding='utf-8')
        region_file.columns = map(str.lower, region_file.columns)
        region_file = region_file[region_file['pro_name_e'] == region]
        all_edge_fail_scenarios = pd.read_excel(province_pth,sheet_name=region.lower().replace(' ', ''))
        all_edge_fail_scenarios = all_edge_fail_scenarios[hazard_cols + ['commune_id','percentage']]
        all_edge_fail_scenarios = all_edge_fail_scenarios.groupby(hazard_cols + ['commune_id'])['percentage'].max().reset_index()

        """Climate change effects
        """
        all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(['hazard_type','commune_id'])
        scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
        change_tup = []
        for sc in scenarios:
            perc = all_edge_fail_scenarios.loc[[sc], 'percentage'].values.tolist()
            yrs = all_edge_fail_scenarios.loc[[sc], 'year'].values.tolist()
            cl = all_edge_fail_scenarios.loc[[sc], 'climate_scenario'].values.tolist()
            if 2016 not in yrs:
                change_tup += list(zip([sc[0]]*len(cl),[sc[1]]*len(cl),cl,yrs,[1e9]*len(cl)))
            elif len(cl) > 1:
                vals = list(zip(cl,perc,yrs))
                vals = sorted(vals, key=lambda pair: pair[-1])
                change = np.array([p for (c,p,y) in vals[1:]]) - vals[0][1]
                cl = [c for (c,p,y) in vals[1:]]
                yrs = [y for (c,p,y) in vals[1:]]
                change_tup += list(zip([sc[0]]*len(cl),[sc[1]]*len(cl),cl,yrs,change))

        change_df = pd.DataFrame(change_tup,columns=['hazard_type','commune_id','climate_scenario','year','change'])
        change_df.to_csv(os.path.join(config['paths']['output'],
            'network_stats',
            '{}_climate_change.csv'.format(region.replace(' ','').lower())
            ), index=False
        )

        """Change effects
        """
        change_df = change_df.set_index(hazard_cols)
        scenarios = list(set(change_df.index.values.tolist()))
        for sc in scenarios:
            hazard_type = sc[0]
            climate_scenario = sc[1]
            year = sc[2]
            percentage = change_df.loc[[sc], 'change'].values.tolist()
            communes = change_df.loc[[sc], 'commune_id'].values.tolist()
            communes_df = pd.DataFrame(list(zip(communes,percentage)),columns=['commune_id','change'])
            commune_vals = pd.merge(region_file,communes_df,how='left',on=['commune_id']).fillna(0)
            del percentage,communes,communes_df

            plot_settings = get_region_plot_settings(region)

            ax = get_axes(plot_settings['bbox'], figsize=plot_settings['figure_size'])

            if region == 'Binh Dinh':
                plot_basemap(ax, config['paths']['data'], country_border='none',
                             plot_states=False, plot_districts=True, highlight_region=region)
            else:
                plot_basemap(ax, config['paths']['data'], country_border='none',
                             plot_states=True, plot_districts=True, highlight_region=region)

            scale_bar(ax, location=(0.8, 0.05), length=plot_settings['scale_legend'])
            proj = ccrs.PlateCarree()

            name = [c['name'] for c in plot_set if c['hazard'] == hazard_type][0]
            for iter_,record in commune_vals.iterrows():
                geom = record.geometry
                region_val = record.change
                if region_val:
                    cl = [c for c in range(len((change_ranges))) if region_val >= change_ranges[c][0] and region_val < change_ranges[c][1]]
                    if cl:
                        c = cl[0]
                        ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                        facecolor=change_colors[c], label=change_labels[c],zorder=2)
                else:
                    ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                    facecolor=change_colors[-1], label=change_labels[-1],zorder=2)

            # Legend
            legend_handles = []
            for c in range(len(change_colors)):
                legend_handles.append(mpatches.Patch(color=change_colors[c], label=change_labels[c]))

            ax.legend(
                handles=legend_handles,
                title='Percentage change in exposure',
                loc='lower left'
            )
            if climate_scenario == 'none':
                climate_scenario = 'current'
            else:
                climate_scenario = climate_scenario.upper()

            # district labels
            plot_district_labels(ax, config['paths']['data'], highlight_region=region)

            plt.title('Percentage change in exposure of road kms for {} {} {}'.format(name,climate_scenario,year), fontsize=14)
            output_file = os.path.join(config['paths']['figures'],
                                       '{}-{}-{}-{}-exposure-change-percentage.png'.format(region.replace(' ',''),name,climate_scenario.replace('.',''),year))
            save_fig(output_file)
            plt.close()

        """Absolute effects
        """
        all_edge_fail_scenarios = all_edge_fail_scenarios.reset_index()
        all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(hazard_cols)
        scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
        for sc in scenarios:
            hazard_type = sc[0]
            climate_scenario = sc[1]
            year = sc[2]
            percentage = all_edge_fail_scenarios.loc[[sc], 'percentage'].values.tolist()
            communes = all_edge_fail_scenarios.loc[[sc], 'commune_id'].values.tolist()
            communes_df = pd.DataFrame(list(zip(communes,percentage)),columns=['commune_id','percentage'])
            commune_vals = pd.merge(region_file,communes_df,how='left',on=['commune_id']).fillna(0)
            del percentage,communes,communes_df

            plot_settings = get_region_plot_settings(region)

            ax = get_axes(plot_settings['bbox'], figsize=plot_settings['figure_size'])

            if region == 'Binh Dinh':
                plot_basemap(ax, config['paths']['data'], country_border='none',
                             plot_states=False, plot_districts=True, highlight_region=region)
            else:
                plot_basemap(ax, config['paths']['data'], country_border='none',
                             plot_states=True, plot_districts=True, highlight_region=region)

            scale_bar(ax, location=(0.8, 0.05), length=plot_settings['scale_legend'])
            proj = ccrs.PlateCarree()

            colors = [c['color'] for c in plot_set if c['hazard'] == hazard_type][0]
            name = [c['name'] for c in plot_set if c['hazard'] == hazard_type][0]

            for iter_,record in commune_vals.iterrows():
                geom = record.geometry
                region_val = record.percentage
                if region_val:
                    if region_val > 0 and region_val <= 10:
                        ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                          facecolor=colors[0], label=labels[0])
                    elif region_val > 10 and region_val <= 20:
                        ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                          facecolor=colors[1], label=labels[1])
                    if region_val > 20 and region_val <= 30:
                        ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                          facecolor=colors[2], label=labels[2])
                    elif region_val > 30 and region_val <= 40:
                        ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                          facecolor=colors[3], label=labels[3])
                    elif region_val > 40 and region_val <= 100:
                        ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                          facecolor=colors[4], label=labels[4])

                else:
                    ax.add_geometries([geom], crs=proj, edgecolor='#ffffff',
                                      facecolor=colors[5], label=labels[5])

            # Legend
            legend_handles = []
            for c in range(len(colors)):
                legend_handles.append(mpatches.Patch(color=colors[c], label=labels[c]))

            ax.legend(
                handles=legend_handles,
                title='Percentage exposure',
                loc='lower left'
            )
            if climate_scenario == 'none':
                climate_scenario = 'current'
            else:
                climate_scenario = climate_scenario.upper()

            # district labels
            plot_district_labels(ax, config['paths']['data'], highlight_region=region)
            
            plt.title('Percentage exposure of road kms for {} {} {}'.format(name,climate_scenario,year), fontsize=14)
            output_file = os.path.join(config['paths']['figures'],
                                       '{}-roads-{}-{}-{}-exposure-percentage.png'.format(region.replace(' ',''),name,climate_scenario.replace('.',''),year))
            save_fig(output_file)
            plt.close()


if __name__ == '__main__':
    main()
