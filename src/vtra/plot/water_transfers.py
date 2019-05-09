"""Water network transfers maps
"""
import os
import sys
from collections import OrderedDict

import numpy as np
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString
from vtra.utils import *


def main(mode):
    config = load_config()
    if mode == 'road':
        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                               'single_edge_failures_transfers_national_road_10_percent_shift.csv')
    elif mode == 'rail':
        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                                   'single_edge_failures_transfers_national_rail_100_percent_shift.csv')
    else:
        raise ValueError("Mode must be road or rail")

    region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                               'coastal_edges.shp')

    region_file = gpd.read_file(region_file_path,encoding='utf-8')
    flow_file = pd.read_csv(flow_file_path)
    region_file = pd.merge(region_file,flow_file,how='left', on=['edge_id']).fillna(0)
    del flow_file



    region_file = region_file[region_file['edge_id'] != 0]

    color = '#045a8d'
    color_by_type = {'Coastal Line': color}

    columns = ['min_tons', 'max_tons']
    column_label_divisors = {c: 1 for c in columns}

    legend_label = "AADF (tons/day)"
    title_cols = ['Total tonnage (min)', 'Total tonnage (max)']

    for c in range(len(columns)):
        ax = get_axes()
        plot_basemap(ax, config['paths']['data'], highlight_region=[])
        scale_bar(ax, location=(0.8, 0.05))
        plot_basemap_labels(ax, config['paths']['data'],plot_international_left=False)
        proj_lat_lon = ccrs.PlateCarree()

        column = columns[c]
        weights = [
            record[column]
            for iter_,record in region_file.iterrows()
        ]
        max_weight = max(weights)
        width_by_range = generate_weight_bins(weights)

        water_geoms_by_category = {
            '1': [],
            '2': []
        }

        geoms_by_range = {}
        for value_range in width_by_range:
            geoms_by_range[value_range] = []

        for iter_,record in region_file.iterrows():
            geom = record.geometry
            val = record[column]
            if val == 0:
                cat = '2'
            else:
                cat = '1'

            buffered_geom = None
            for (nmin, nmax), width in width_by_range.items():
                if nmin <= val and val < nmax:
                    buffered_geom = geom.buffer(width)

            if buffered_geom is not None:
                water_geoms_by_category[cat].append(buffered_geom)
            else:
                print("Feature was outside range to plot", iter_)


        styles = OrderedDict([
            ('1',  Style(color='#045a8d', zindex=9, label='Transfers')),  # blue
            ('2', Style(color='#969696', zindex=7, label='No transfer'))
        ])

        for cat, geoms in water_geoms_by_category.items():
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
            legend_label,
            horizontalalignment='left',
            transform=proj_lat_lon,
            size=10)

        divisor = column_label_divisors[column]
        for (i, ((nmin, nmax), width)) in enumerate(width_by_range.items()):
            y = base_y - (i*y_step)
            line = LineString([(x_l, y), (x_r, y)])
            ax.add_geometries(
                [line.buffer(width)],
                crs=proj_lat_lon,
                linewidth=0,
                edgecolor=color,
                facecolor=color,
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

        plt.title('{} transfers - '.format(mode.title()) + title_cols[c], fontsize=14)
        legend_from_style_spec(ax, styles)
        print ('* Plotting ',title_cols[c])
        if mode == 'road':
            output_file = os.path.join(
                config['paths']['figures'], 'water_flow-map-transfer-{}-10-shift-{}.png'.format(mode,column))
        elif mode == 'rail':
            output_file = os.path.join(
                config['paths']['figures'], 'water_flow-map-transfer-{}-100-shift-{}.png'.format(mode,column))
        else:
            raise ValueError("Mode must be road or rail")
        save_fig(output_file)
        plt.close()


if __name__ == '__main__':
    ok_values = ('road', 'rail')
    for ok in ok_values:
        main(ok)
