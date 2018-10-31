"""Road network flows
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
from shapely.geometry import LineString
from vtra.utils import *


def main(mode):
    config = load_config()
    if mode == 'road':
        region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                                   'road_edges.shp')
        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                                   'single_edge_failures_minmax_national_road_10_percent_modal_shift.csv')
    elif mode == 'rail':
        region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                                   'rail_edges.shp')
        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                                   'single_edge_failures_minmax_national_rail_100_percent_disrupt_multi_modal.csv')
    else:
        raise ValueError("Mode must be road or rail")

    region_file = gpd.read_file(region_file_path,encoding='utf-8')
    flow_file = pd.read_csv(flow_file_path)
    region_file = pd.merge(region_file,flow_file,how='left', on=['edge_id']).fillna(0)
    del flow_file

    plot_sets = [
        {
            'file_tag': 'loss',
            'no_access': [0, 1],
            'legend_label': "Economic loss (million USD/day)",
            'divisor': 1000000,
            'columns': ['min_econ_impact', 'max_econ_impact'],
            'title_cols': ['Economic impact (min)', 'Economic impact (max)']
        },
    ]

    for plot_set in plot_sets:
        for c in range(len(plot_set['columns'])):
            column = plot_set['columns'][c]

            ax = get_axes()
            plot_basemap(ax, config['paths']['data'], highlight_region=[])
            scale_bar(ax, location=(0.8, 0.05))
            plot_basemap_labels(ax, config['paths']['data'],plot_international_left=False)
            proj_lat_lon = ccrs.PlateCarree()


            weights = [
                record[column]
                for iter_,record in region_file.iterrows()
            ]

            min_weight = min(weights)
            max_weight = max(weights)
            abs_max_weight = max([abs(w) for w in weights])

            width_by_range = OrderedDict()
            colors_by_range = {}
            n_steps = 8

            positive_colors = [
                '#f4a582',
                '#d6604d',
                '#b2182b',
                '#67001f',
            ]
            negative_colors = [
                '#92c5de',
                '#4393c3',
                '#2166ac',
                '#053061',
            ]
            width_step = 0.01

            mins = np.linspace(0, abs_max_weight, n_steps/2)

            maxs = list(mins)
            maxs.append(abs_max_weight*10)
            maxs = maxs[1:]

            assert len(maxs) == len(mins)

            # positive
            for i, (min_, max_) in reversed(list(enumerate(zip(mins, maxs)))):
                width_by_range[(min_, max_)] = (i + 2) * width_step
                colors_by_range[(min_, max_)] = positive_colors[i]

            # negative
            for i, (min_, max_) in enumerate(zip(mins, maxs)):
                width_by_range[(-max_, -min_)] = (i + 2) * width_step
                colors_by_range[(-max_, -min_)] = negative_colors[i]


            geoms_by_range = {}
            for value_range in width_by_range:
                geoms_by_range[value_range] = []

            zero_value_geoms = []
            for iter_,record in region_file.iterrows():
                val = record[column]
                geom = record.geometry
                if val != 0:
                    for nmin, nmax in geoms_by_range:
                        if nmin <= val and val < nmax:
                            geoms_by_range[(nmin, nmax)].append(geom)
                else:
                    zero_value_geoms.append(geom)

            # plot
            for range_, width in width_by_range.items():
                ax.add_geometries(
                    [geom.buffer(width) for geom in geoms_by_range[range_]],
                    crs=proj_lat_lon,
                    edgecolor='none',
                    facecolor=colors_by_range[range_],
                    zorder=2)

            width_min = min([width for range_, width in width_by_range.items()])
            ax.add_geometries(
                [geom.buffer(width_min) for geom in zero_value_geoms],
                crs=proj_lat_lon,
                edgecolor='none',
                facecolor='#969696',
                zorder=1)

            x_l = 102.3
            x_r = x_l + 0.4
            base_y = 14
            y_step = 0.4
            y_text_nudge = 0.1
            x_text_nudge = 0.1

            ax.text(
                x_l - x_text_nudge,
                base_y + y_step - y_text_nudge,
                plot_set['legend_label'],
                horizontalalignment='left',
                transform=proj_lat_lon,
                size=8)

            divisor = plot_set['divisor']

            i = 0
            for (nmin, nmax), width in width_by_range.items():
                if not geoms_by_range[(nmin, nmax)]:
                    continue
                y = base_y - (i*y_step)
                i = i + 1
                line = LineString([(x_l, y), (x_r, y)])
                ax.add_geometries(
                    [line.buffer(width)],
                    crs=proj_lat_lon,
                    linewidth=0,
                    edgecolor=colors_by_range[(nmin, nmax)],
                    facecolor=colors_by_range[(nmin, nmax)],
                    zorder=2)
                if nmin == max_weight:
                    label = '>{:.2f}'.format(max_weight/divisor)
                elif nmax == -abs_max_weight:
                    label = '<{:.2f}'.format(-abs_max_weight/divisor)
                else:
                    label = '{:.2f} to {:.2f}'.format(nmin/divisor, nmax/divisor)
                ax.text(
                    x_r + x_text_nudge,
                    y - y_text_nudge,
                    label,
                    horizontalalignment='left',
                    transform=proj_lat_lon,
                    size=8)

            styles = OrderedDict([
                ('1',  Style(color='#b2182b', zindex=9, label='Economic loss effect')),  # green
                ('2',  Style(color='#2166ac', zindex=9, label='Economic gain effect')),
                ('3', Style(color='#969696', zindex=9, label='No hazard exposure/effect'))
            ])
            plt.title(plot_set['title_cols'][c], fontsize=14)
            legend_from_style_spec(ax, styles,loc='center left')

            print ('* Plotting {} {}'.format(mode,plot_set['title_cols'][c]))
            if mode == 'road':
                output_file = os.path.join(
                    config['paths']['figures'], 'road_failure-map-{}-{}-multi-modal-options-10-shift.png'.format(plot_set['file_tag'], column))
            elif mode == 'rail':
                output_file = os.path.join(
                    config['paths']['figures'], 'rail_failure-map-{}-{}-multi-modal-options.png'.format(plot_set['file_tag'], column))
            else:
                raise ValueError("Mode must be road or rail")
            save_fig(output_file)
            plt.close()
            print(" >", output_file)


if __name__ == '__main__':
    ok_values = ('road', 'rail')
    for ok in ok_values:
        main(ok)
