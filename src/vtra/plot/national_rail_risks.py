"""Rail hazard exposure maps
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
    change_colors = ['#1a9850','#66bd63','#a6d96a','#d9ef8b','#fee08b','#fdae61','#f46d43','#d73027','#969696']
    change_labels = ['< -40','-40 to -20','-20 to -10','-10 to 0','0 to 10','10 to 20','20 to 40',' > 40','No change/value']
    change_ranges = [(-1e10,-40),(-40,-20),(-20,-10),(-10,0),(0.001,10),(10,20),(20,40),(40,1e10)]

    eael_set = [
        {
            'column': 'min_eael',
            'title': 'Min EAEL',
            'legend_label': "Expected Annual losses (million USD)",
            'divisor': 1000000,
            'significance': 0
        },
        {
            'column': 'max_eael',
            'title': 'Max EAEL',
            'legend_label': "Expected Annual losses (million USD)",
            'divisor': 1000000,
            'significance': 0
        }
    ]

    adapt_set = [
        {
            'column': 'min_eael',
            'title': 'Min EAEL',
            'legend_label': "Expected Annual losses (million USD)",
            'divisor': 1000000,
            'significance': 0
        },
        {
            'column': 'max_eael',
            'title': 'Max EAEL',
            'legend_label': "Expected Annual losses (million USD)",
            'divisor': 1000000,
            'significance': 0
        },
    ]

    region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                               'rail_edges.shp')

    flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                               'single_edge_failures_minmax_national_rail_100_percent_disrupt.csv')

    region_file = gpd.read_file(region_file_path,encoding='utf-8')
    flow_file = pd.read_csv(flow_file_path)
    region_file = pd.merge(region_file,flow_file,how='left', on=['edge_id']).fillna(0)
    del flow_file

    flow_file_path = os.path.join(config['paths']['output'], 'hazard_scenarios',
                               'national_rail_hazard_intersections_risks.csv')

    fail_scenarios = pd.read_csv(flow_file_path)
    fail_scenarios = pd.merge(fail_scenarios,region_file[['edge_id','min_econ_impact','max_econ_impact']],how='left', on=['edge_id']).fillna(0)
    fail_scenarios['min_eael'] = duration*fail_scenarios['min_duration_wt']*fail_scenarios['risk_wt']*fail_scenarios['min_econ_impact']
    fail_scenarios['max_eael'] = duration*fail_scenarios['max_duration_wt']*fail_scenarios['risk_wt']*fail_scenarios['max_econ_impact']
    all_edge_fail_scenarios = fail_scenarios[hazard_cols + ['edge_id','min_eael','max_eael']]
    all_edge_fail_scenarios = all_edge_fail_scenarios.groupby(hazard_cols + ['edge_id'])['min_eael','max_eael'].max().reset_index()

    # Climate change effects
    all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(['hazard_type','edge_id'])
    scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
    change_tup = []
    for sc in scenarios:
        eael = all_edge_fail_scenarios.loc[[sc], 'max_eael'].values.tolist()
        yrs = all_edge_fail_scenarios.loc[[sc], 'year'].values.tolist()
        cl = all_edge_fail_scenarios.loc[[sc], 'climate_scenario'].values.tolist()
        if 2016 not in yrs:
            change_tup += list(zip([sc[0]]*len(cl),[sc[1]]*len(cl),cl,yrs,[0]*len(cl),eael,[1e9]*len(cl)))
        elif len(yrs) > 1:
            vals = list(zip(cl,eael,yrs))
            vals = sorted(vals, key=lambda pair: pair[-1])
            change = 100.0*(np.array([p for (c,p,y) in vals[1:]]) - vals[0][1])/vals[0][1]
            cl = [c for (c,p,y) in vals[1:]]
            yrs = [y for (c,p,y) in vals[1:]]
            fut = [p for (c,p,y) in vals[1:]]
            change_tup += list(zip([sc[0]]*len(cl),[sc[1]]*len(cl),cl,yrs,[vals[0][1]]*len(cl),fut,change))

    change_df = pd.DataFrame(change_tup,columns=['hazard_type','edge_id','climate_scenario','year','current','future','change'])
    change_df.to_csv(os.path.join(config['paths']['output'],
        'network_stats',
        'national_rails_eael_climate_change.csv'
        ), index=False
    )

    # # Change effects
    # change_df = change_df.set_index(hazard_cols)
    # scenarios = list(set(change_df.index.values.tolist()))
    # for sc in scenarios:
    #     hazard_type = sc[0]
    #     climate_scenario = sc[1]
    #     year = sc[2]
    #     percentage = change_df.loc[[sc], 'change'].values.tolist()
    #     edges = change_df.loc[[sc], 'edge_id'].values.tolist()
    #     edges_df = pd.DataFrame(list(zip(edges,percentage)),columns=['edge_id','change'])
    #     edges_vals = pd.merge(region_file,edges_df,how='left',on=['edge_id']).fillna(0)
    #     del percentage,edges,edges_df

    #     ax = get_axes()
    #     plot_basemap(ax, config['paths']['data'], highlight_region=[])
    #     scale_bar(ax, location=(0.8, 0.05))
    #     plot_basemap_labels(ax, config['paths']['data'],plot_international_left=False)
    #     proj = ccrs.PlateCarree()

    #     name = [c['name'] for c in hazard_set if c['hazard'] == hazard_type][0]
    #     for iter_,record in edges_vals.iterrows():
    #         geom = record.geometry
    #         region_val = record.change
    #         if region_val:
    #             cl = [c for c in range(len((change_ranges))) if region_val >= change_ranges[c][0] and region_val < change_ranges[c][1]]
    #             if cl:
    #                 c = cl[0]
    #                 ax.add_geometries([geom],crs=proj,linewidth=1.5,edgecolor=change_colors[c],facecolor='none',zorder=2)
    #         else:
    #             ax.add_geometries([geom], crs=proj, linewidth=1.5,edgecolor=change_colors[-1],facecolor='none',zorder=1)


    #     # Legend
    #     legend_handles = []
    #     for c in range(len(change_colors)):
    #         legend_handles.append(mpatches.Patch(color=change_colors[c], label=change_labels[c]))

    #     ax.legend(
    #         handles=legend_handles,
    #         title='Percentage change in EAEL',
    #         loc='center left'
    #     )
    #     if climate_scenario == 'none':
    #         climate_scenario = 'current'
    #     else:
    #         climate_scenario = climate_scenario.upper()

    #     title = 'Percentage change in EAEL for {} {} {}'.format(name,climate_scenario,year)
    #     print(" * Plotting {}".format(title))

    #     plt.title(title, fontsize=14)
    #     output_file = os.path.join(config['paths']['figures'],
    #                                'national-rail-{}-{}-{}-risks-change-percentage.png'.format(name,climate_scenario.replace('.',''),year))
    #     save_fig(output_file)
    #     plt.close()

    # # Absolute effects
    # all_edge_fail_scenarios = all_edge_fail_scenarios.reset_index()
    # all_edge_fail_scenarios = all_edge_fail_scenarios.set_index(hazard_cols)
    # scenarios = list(set(all_edge_fail_scenarios.index.values.tolist()))
    # for sc in scenarios:
    #     hazard_type = sc[0]
    #     climate_scenario = sc[1]
    #     if climate_scenario == 'none':
    #         climate_scenario = 'current'
    #     else:
    #         climate_scenario = climate_scenario.upper()
    #     year = sc[2]
    #     min_eael = all_edge_fail_scenarios.loc[[sc], 'min_eael'].values.tolist()
    #     max_eael = all_edge_fail_scenarios.loc[[sc], 'max_eael'].values.tolist()
    #     edges = all_edge_fail_scenarios.loc[[sc], 'edge_id'].values.tolist()
    #     edges_df = pd.DataFrame(list(zip(edges,min_eael,max_eael)),columns=['edge_id','min_eael','max_eael'])
    #     edges_vals = pd.merge(region_file,edges_df,how='left',on=['edge_id']).fillna(0)
    #     del edges_df

    #     for c in range(len(eael_set)):
    #         ax = get_axes()
    #         plot_basemap(ax, config['paths']['data'], highlight_region=[])
    #         scale_bar(ax, location=(0.8, 0.05))
    #         plot_basemap_labels(ax, config['paths']['data'],plot_international_left=False)
    #         proj_lat_lon = ccrs.PlateCarree()

    #         # generate weight bins
    #         column = eael_set[c]['column']
    #         weights = [record[column] for iter_, record in edges_vals.iterrows()]

    #         max_weight = max(weights)
    #         width_by_range = generate_weight_bins(weights)

    #         rail_geoms_by_category = {
    #             '1': [],
    #             '2': []
    #         }

    #         for iter_,record in edges_vals.iterrows():
    #             geom = record.geometry
    #             val = record[column]
    #             if val == 0:
    #                 cat = '2'
    #             else:
    #                 cat = '1'

    #             buffered_geom = None
    #             for (nmin, nmax), width in width_by_range.items():
    #                 if nmin <= val and val < nmax:
    #                     buffered_geom = geom.buffer(width)

    #             if buffered_geom is not None:
    #                 rail_geoms_by_category[cat].append(buffered_geom)
    #             else:
    #                 print("Feature was outside range to plot", iter_)

    #         styles = OrderedDict([
    #             ('1',  Style(color='#006d2c', zindex=9, label='Hazard failure effect')),  # green
    #             ('2', Style(color='#969696', zindex=7, label='No hazard exposure/effect'))
    #         ])

    #         for cat, geoms in rail_geoms_by_category.items():
    #             cat_style = styles[cat]
    #             ax.add_geometries(
    #                 geoms,
    #                 crs=proj_lat_lon,
    #                 linewidth=0,
    #                 facecolor=cat_style.color,
    #                 edgecolor='none',
    #                 zorder=cat_style.zindex
    #             )
    #         name = [h['name'] for h in hazard_set if h['hazard'] == hazard_type][0]

    #         x_l = 102.3
    #         x_r = x_l + 0.4
    #         base_y = 14
    #         y_step = 0.4
    #         y_text_nudge = 0.1
    #         x_text_nudge = 0.1

    #         ax.text(
    #             x_l,
    #             base_y + y_step - y_text_nudge,
    #             eael_set[c]['legend_label'],
    #             horizontalalignment='left',
    #             transform=proj_lat_lon,
    #             size=10)

    #         divisor = eael_set[c]['divisor']
    #         significance_ndigits = eael_set[c]['significance']
    #         max_sig = []
    #         for (i, ((nmin, nmax), line_style)) in enumerate(width_by_range.items()):
    #             if round(nmin/divisor, significance_ndigits) < round(nmax/divisor, significance_ndigits):
    #                 max_sig.append(significance_ndigits)
    #             elif round(nmin/divisor, significance_ndigits+1) < round(nmax/divisor, significance_ndigits+1):
    #                 max_sig.append(significance_ndigits+1)
    #             elif round(nmin/divisor, significance_ndigits+2) < round(nmax/divisor, significance_ndigits+2):
    #                 max_sig.append(significance_ndigits+2)
    #             else:
    #                 max_sig.append(significance_ndigits+3)

    #         significance_ndigits = max(max_sig)
    #         for (i, ((nmin, nmax), width)) in enumerate(width_by_range.items()):
    #             y = base_y - (i*y_step)
    #             line = LineString([(x_l, y), (x_r, y)]).buffer(width)
    #             ax.add_geometries(
    #                 [line],
    #                 crs=proj_lat_lon,
    #                 linewidth=0,
    #                 edgecolor='#000000',
    #                 facecolor='#000000',
    #                 zorder=2)
    #             if nmin == max_weight:
    #                 value_template = '>{:.' + str(significance_ndigits) + 'f}'
    #                 label = value_template.format(
    #                     round(max_weight/divisor, significance_ndigits))
    #             else:
    #                 value_template = '{:.' + str(significance_ndigits) + \
    #                     'f}-{:.' + str(significance_ndigits) + 'f}'
    #                 label = value_template.format(
    #                     round(nmin/divisor, significance_ndigits), round(nmax/divisor, significance_ndigits))

    #             ax.text(
    #                 x_r + x_text_nudge,
    #                 y - y_text_nudge,
    #                 label,
    #                 horizontalalignment='left',
    #                 transform=proj_lat_lon,
    #                 size=10)

    #         title = 'National rail ({}) {} {} {}'.format(eael_set[c]['title'],name,climate_scenario,year)
    #         print ('* Plotting ',title)

    #         plt.title(title, fontsize=14)
    #         legend_from_style_spec(ax, styles,loc='center left')

    #         # output
    #         output_file = os.path.join(
    #             config['paths']['figures'], 'national-rails-{}-{}-{}-{}.png'.format(name,climate_scenario.replace('.',''),year,eael_set[c]['column']))
    #         save_fig(output_file)
    #         plt.close()

    # all_edge_fail_scenarios = fail_scenarios[['edge_id','min_eael','max_eael']]
    # all_edge_fail_scenarios = all_edge_fail_scenarios.groupby(['edge_id'])[['min_eael','max_eael']].max().reset_index()
    # edges_vals = pd.merge(region_file,all_edge_fail_scenarios,how='left',on=['edge_id']).fillna(0)

    # for c in range(len(adapt_set)):
    #     ax = get_axes()
    #     plot_basemap(ax, config['paths']['data'], highlight_region=[])
    #     scale_bar(ax, location=(0.8, 0.05))
    #     plot_basemap_labels(ax, config['paths']['data'],plot_international_left=False)
    #     proj_lat_lon = ccrs.PlateCarree()

    #     # generate weight bins
    #     column = adapt_set[c]['column']
    #     weights = [record[column] for iter_, record in edges_vals.iterrows()]


    #     max_weight = max(weights)
    #     width_by_range = generate_weight_bins(weights)

    #     rail_geoms_by_category = {
    #         '1': [],
    #         '2': []
    #     }

    #     for iter_,record in edges_vals.iterrows():
    #         geom = record.geometry
    #         val = record[column]
    #         if val == 0:
    #             cat = '2'
    #         else:
    #             cat = '1'

    #         buffered_geom = None
    #         for (nmin, nmax), width in width_by_range.items():
    #             if nmin <= val and val < nmax:
    #                 buffered_geom = geom.buffer(width)

    #         if buffered_geom is not None:
    #             rail_geoms_by_category[cat].append(buffered_geom)
    #         else:
    #             print("Feature was outside range to plot", iter_)

    #     styles = OrderedDict([
    #         ('1',  Style(color='#006d2c', zindex=9, label='Hazard failure effect')),  # green
    #         ('2', Style(color='#969696', zindex=7, label='No hazard exposure/effect'))
    #     ])

    #     for cat, geoms in rail_geoms_by_category.items():
    #         cat_style = styles[cat]
    #         ax.add_geometries(
    #             geoms,
    #             crs=proj_lat_lon,
    #             linewidth=0,
    #             facecolor=cat_style.color,
    #             edgecolor='none',
    #             zorder=cat_style.zindex
    #         )


    #     x_l = 102.3
    #     x_r = x_l + 0.4
    #     base_y = 14
    #     y_step = 0.4
    #     y_text_nudge = 0.1
    #     x_text_nudge = 0.1

    #     ax.text(
    #         x_l,
    #         base_y + y_step - y_text_nudge,
    #         adapt_set[c]['legend_label'],
    #         horizontalalignment='left',
    #         transform=proj_lat_lon,
    #         size=10)

    #     divisor = adapt_set[c]['divisor']
    #     significance_ndigits = adapt_set[c]['significance']
    #     max_sig = []
    #     for (i, ((nmin, nmax), line_style)) in enumerate(width_by_range.items()):
    #         if round(nmin/divisor, significance_ndigits) < round(nmax/divisor, significance_ndigits):
    #             max_sig.append(significance_ndigits)
    #         elif round(nmin/divisor, significance_ndigits+1) < round(nmax/divisor, significance_ndigits+1):
    #             max_sig.append(significance_ndigits+1)
    #         elif round(nmin/divisor, significance_ndigits+2) < round(nmax/divisor, significance_ndigits+2):
    #             max_sig.append(significance_ndigits+2)
    #         else:
    #             max_sig.append(significance_ndigits+3)

    #     significance_ndigits = max(max_sig)
    #     for (i, ((nmin, nmax), width)) in enumerate(width_by_range.items()):
    #         y = base_y - (i*y_step)
    #         line = LineString([(x_l, y), (x_r, y)]).buffer(width)
    #         ax.add_geometries(
    #             [line],
    #             crs=proj_lat_lon,
    #             linewidth=0,
    #             edgecolor='#000000',
    #             facecolor='#000000',
    #             zorder=2)
    #         if nmin == max_weight:
    #             value_template = '>{:.' + str(significance_ndigits) + 'f}'
    #             label = value_template.format(
    #                 round(max_weight/divisor, significance_ndigits))
    #         else:
    #             value_template = '{:.' + str(significance_ndigits) + \
    #                 'f}-{:.' + str(significance_ndigits) + 'f}'
    #             label = value_template.format(
    #                 round(nmin/divisor, significance_ndigits), round(nmax/divisor, significance_ndigits))

    #         ax.text(
    #             x_r + x_text_nudge,
    #             y - y_text_nudge,
    #             label,
    #             horizontalalignment='left',
    #             transform=proj_lat_lon,
    #             size=10)


    #     # plot
    #     title = 'National rail ({})'.format(adapt_set[c]['title'])
    #     print(" * Plotting", title)
    #     plt.title(title, fontsize=14)
    #     legend_from_style_spec(ax, styles,loc='center left')

    #     # output
    #     output_file = os.path.join(
    #         config['paths']['figures'], 'national_rail-{}-values.png'.format(column))
    #     save_fig(output_file)
    #     plt.close()

if __name__ == '__main__':
    main()
