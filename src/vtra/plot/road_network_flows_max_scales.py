"""Road network flow maps
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
    mode_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                                   'road_edges.shp')
    flow_file_path = os.path.join(config['paths']['output'], 'flow_mapping_combined',
                                   'weighted_flows_national_road_100_percent.csv')


    mode_file = gpd.read_file(mode_file_path,encoding='utf-8')
    flow_file = pd.read_csv(flow_file_path)
    mode_file = pd.merge(mode_file,flow_file,how='left', on=['edge_id']).fillna(0)

    plot_sets = [
        {
            'file_tag': 'vehicle_count',
            'legend_label': "AADT ('000 vehicles/day)",
            'divisor': 1000,
            'columns': ['vehicle_co'],
            'title_cols': ['Vehicle Count']
        },
        {
            'file_tag': 'commodities',
            'legend_label': "AADF ('000 tons/day)",
            'divisor': 1000,
            'columns': ['max_rice', 'max_cash', 'max_cass', 'max_teas', 'max_maiz',
                        'max_rubb', 'max_swpo', 'max_acof', 'max_rcof', 'max_pepp',
                        'max_sugar', 'max_wood', 'max_steel', 'max_constructi', 'max_cement',
                        'max_fertilizer', 'max_coal', 'max_petroluem', 'max_manufactur', 'max_fishery',
                        'max_meat', 'max_tons'],
            'title_cols': ['Rice', 'Cashew', 'Cassava', 'Teas', 'Maize', 'Rubber',
                           'Sweet Potatoes', 'Coffee Arabica', 'Coffee Robusta',
                           'Pepper', 'Sugar', 'Wood', 'Steel', 'Construction materials',
                           'Cement', 'Fertilizer', 'Coal', 'Petroleum', 'Manufacturing',
                           'Fishery', 'Meat', 'Total tonnage']
        },
    ]

    for plot_set in plot_sets:
        for c in range(len(plot_set['columns'])):
            ax = get_axes()
            plot_basemap(ax, config['paths']['data'],highlight_region=[])
            scale_bar(ax, location=(0.8, 0.05))
            plot_basemap_labels(ax, config['paths']['data'])
            proj_lat_lon = ccrs.PlateCarree()

            # generate weight bins
            if plot_set['columns'][c] == 'vehicle_co':
                column = plot_set['columns'][c]
            else:
                column = 'max_tons'

            weights = [
                record[column]
                for iter_, record in mode_file.iterrows()
            ]
            max_weight = max(weights)
            width_by_range = generate_weight_bins(weights)

            road_geoms_by_category = {
                '1': [],
                '2': [],
                '3': [],
                '4': [],
                '5': [],
                '6': []
            }

            column = plot_set['columns'][c]
            for iter_, record in mode_file.iterrows():
                cat = str(record['road_class'])
                if cat not in road_geoms_by_category:
                    raise Exception
                geom = record.geometry

                val = record[column]

                buffered_geom = None
                for (nmin, nmax), width in width_by_range.items():
                    if nmin <= val and val < nmax:
                        buffered_geom = geom.buffer(width)

                if buffered_geom is not None:
                    road_geoms_by_category[cat].append(buffered_geom)
                else:
                    print("Feature was outside range to plot", record.attributes)

            styles = OrderedDict([
                ('1',  Style(color='#000004', zindex=9, label='Class 1')),  # red
                ('2', Style(color='#2c115f', zindex=8, label='Class 2')),  # orange
                ('3', Style(color='#721f81', zindex=7, label='Class 3')),  # blue
                ('4',  Style(color='#b73779', zindex=6, label='Class 4')),  # green
                ('5', Style(color='#f1605d', zindex=5, label='Class 5')),  # black
                ('6', Style(color='#feb078', zindex=4, label='Class 6')),  # grey
            ])

            for cat, geoms in road_geoms_by_category.items():
                cat_style = styles[cat]
                ax.add_geometries(
                    geoms,
                    crs=proj_lat_lon,
                    linewidth=0,
                    facecolor=cat_style.color,
                    edgecolor='none',
                    zorder=cat_style.zindex
                )

            x_l = 102.3
            x_r = x_l + 0.4
            base_y = 14
            y_step = 0.4
            y_text_nudge = 0.1
            x_text_nudge = 0.1

            ax.text(
                x_l,
                base_y + y_step - y_text_nudge,
                plot_set['legend_label'],
                horizontalalignment='left',
                transform=proj_lat_lon,
                size=10)

            divisor = plot_set['divisor']
            for (i, ((nmin, nmax), width)) in enumerate(width_by_range.items()):
                y = base_y - (i*y_step)
                line = LineString([(x_l, y), (x_r, y)]).buffer(width)
                ax.add_geometries(
                    [line],
                    crs=proj_lat_lon,
                    linewidth=0,
                    edgecolor='#000000',
                    facecolor='#000000',
                    zorder=2)
                if nmin == max_weight:
                    label = '>{:.2f}'.format(max_weight/divisor)
                else:
                    label = '{:.2f}-{:.2f}'.format(nmin/divisor, nmax/divisor)
                ax.text(
                    x_r + x_text_nudge,
                    y - y_text_nudge,
                    label,
                    horizontalalignment='left',
                    transform=proj_lat_lon,
                    size=10)

            plt.title(plot_set['title_cols'][c], fontsize=14)
            legend_from_style_spec(ax, styles)
            output_file = os.path.join(
                config['paths']['figures'],
                'road_flow-map-{}-{}-max-scale.png'.format(plot_set['file_tag'], column))
            save_fig(output_file)
            plt.close()


if __name__ == '__main__':
    main()
