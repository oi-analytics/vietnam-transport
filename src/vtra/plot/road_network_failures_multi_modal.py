"""Road network failure maps (multi-modal)
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

    region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                               'road_edges.shp')

    flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                               'single_edge_failures_minmax_national_road_10_percent_modal_shift.csv')

    region_file = gpd.read_file(region_file_path,encoding='utf-8')
    flow_file = pd.read_csv(flow_file_path)
    region_file = pd.merge(region_file,flow_file,how='left', on=['edge_id']).fillna(0)
    del flow_file

    plot_set = [
        {
            'column': 'min_econ_impact',
            'title': 'Min Total Economic loss',
            'legend_label': "Economic Loss (million USD/day)",
            'divisor': 1000000,
            'significance': 0
        },
        {
            'column': 'max_econ_impact',
            'title': 'Max Total Economic loss',
            'legend_label': "Economic Loss (million USD/day)",
            'divisor': 1000000,
            'significance': 0
        }
    ]

    for c in range(len(plot_set)):
        ax = get_axes()
        plot_basemap(ax, config['paths']['data'], highlight_region=[])
        scale_bar(ax, location=(0.8, 0.05))
        plot_basemap_labels(ax, config['paths']['data'],plot_international_left=False)
        proj_lat_lon = ccrs.PlateCarree()

        # generate weight bins
        column = plot_set[c]['column']
        weights = [
            record[column]
            for iter_,record in region_file.iterrows()
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

        for iter_,record in region_file.iterrows():
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
            plot_set[c]['legend_label'],
            horizontalalignment='left',
            transform=proj_lat_lon,
            size=10)

        divisor = plot_set[c]['divisor']
        significance_ndigits = plot_set[c]['significance']
        max_sig = []
        for (i, ((nmin, nmax), line_style)) in enumerate(width_by_range.items()):
            if round(nmin/divisor, significance_ndigits) < round(nmax/divisor, significance_ndigits):
                max_sig.append(significance_ndigits)
            elif round(nmin/divisor, significance_ndigits+1) < round(nmax/divisor, significance_ndigits+1):
                max_sig.append(significance_ndigits+1)
            elif round(nmin/divisor, significance_ndigits+2) < round(nmax/divisor, significance_ndigits+2):
                max_sig.append(significance_ndigits+2)
            else:
                max_sig.append(significance_ndigits+3)

        significance_ndigits = max(max_sig)
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
                value_template = '>{:.' + str(significance_ndigits) + 'f}'
                label = value_template.format(
                    round(max_weight/divisor, significance_ndigits))
            else:
                value_template = '{:.' + str(significance_ndigits) + \
                    'f}-{:.' + str(significance_ndigits) + 'f}'
                label = value_template.format(
                    round(nmin/divisor, significance_ndigits), round(nmax/divisor, significance_ndigits))
            ax.text(
                x_r + x_text_nudge,
                y - y_text_nudge,
                label,
                horizontalalignment='left',
                transform=proj_lat_lon,
                size=10)

        plt.title(plot_set[c]['title'], fontsize=14)
        legend_from_style_spec(ax, styles)
        output_file = os.path.join(
            config['paths']['figures'],
            'road_failure-map-{}-multi-modal-options.png'.format(
                column))
        save_fig(output_file)
        plt.close()


if __name__ == '__main__':
    main()
