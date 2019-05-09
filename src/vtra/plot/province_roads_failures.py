"""Provincial road network loss maps
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
    plot_set = [
        {
            'column': 'min_tr_loss',
            'title': 'Min Rerouting loss',
            'legend_label': "Rerouting Loss (USD/day)",
            'divisor': 1,
            'significance': 0
        },
        {
            'column': 'max_tr_loss',
            'title': 'Max Rerouting loss',
            'legend_label': "Rerouting Loss (USD/day)",
            'divisor': 1,
            'significance': 0
        },
        {
            'column': 'min_netrev',
            'title': 'Min Net revenue disrupted',
            'legend_label': "Net revenue disrupted ('000 USD/day)",
            'divisor': 1000,
            'significance': 0
        },
        {
            'column': 'max_netrev',
            'title': 'Max Net revenue disrupted',
            'legend_label': "Net revenue disrupted ('000 USD/day)",
            'divisor': 1000,
            'significance': 0
        },
        {
            'column': 'min_econ_impact',
            'title': 'Min Economic loss',
            'legend_label': "Economic Loss ('000 USD/day)",
            'divisor': 1000,
            'significance': 0
        },
        {
            'column': 'max_econ_impact',
            'title': 'Max Economic loss',
            'legend_label': "Economic Loss ('000 USD/day)",
            'divisor': 1000,
            'significance': 0
        },
        {
            'column': 'min_croptons',
            'title': 'Min Daily Crop Tons disrupted',
            'legend_label': "Crops transport disrupted (tons/day)",
            'divisor': 1,
            'significance': 0
        },
        {
            'column': 'max_croptons',
            'title': 'Max Daily Crop Tons disrupted',
            'legend_label': "Crops transport disrupted (tons/day)",
            'divisor': 1,
            'significance': 0
        }
    ]

    for region in regions:
        region_file_path = os.path.join(config['paths']['data'], 'post_processed_networks',
                                   '{}_roads_edges.shp'.format(region.lower().replace(' ', '')))
        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                                   'single_edge_failures_minmax_{}_5_tons_100_percent_disrupt.csv'.format(region.lower().replace(' ', '')))

        region_file = gpd.read_file(region_file_path,encoding='utf-8')
        flow_file = pd.read_csv(flow_file_path)
        region_file = pd.merge(region_file,flow_file,how='left', on=['edge_id']).fillna(0)

        plot_settings = get_region_plot_settings(region)

        for c in range(len(plot_set)):
            ax = get_axes(plot_settings['bbox'], figsize=plot_settings['figure_size'])

            if region == 'Binh Dinh':
                plot_basemap(ax, config['paths']['data'], country_border='none',
                             plot_states=False, plot_districts=True, highlight_region=region)
            else:
                plot_basemap(ax, config['paths']['data'], country_border='none',
                             plot_states=True, plot_districts=True, highlight_region=region)

            scale_bar(ax, location=(0.8, 0.05), length=plot_settings['scale_legend'])
            proj_lat_lon = ccrs.PlateCarree()

            # generate weight bins
            column = plot_set[c]['column']
            weights = [record[column] for iter_, record in region_file.iterrows()]

            max_weight = max(weights)
            if column in ('min_bc_rat', 'max_bc_rat') and region == 'Lao Cai':
                width_by_range = generate_weight_bins_with_colour_gradient(
                    weights, n_steps=6, width_step=0.001)
            else:
                width_by_range = generate_weight_bins_with_colour_gradient(
                    weights, width_step=0.001)

            road_geoms_by_category = {
                region: []
            }

            styles = OrderedDict([
                ('1',  Style(color='#252525', zindex=6, label='Hazard failure: No access')),  # black
                ('2', Style(color='#f1605d', zindex=8, label='Hazard failure: Reroutiing possible')),  # orange
                ('3', Style(color='#969696', zindex=7, label='No hazard exposure/effect')) # grey
            ])

            for iter_, record in region_file.iterrows():
                cat = region
                geom = record.geometry
                noaccess = record.no_access

                val = record[column]

                buffered_geom = None
                for (nmin, nmax), line_style in width_by_range.items():
                    if nmin <= val and val < nmax:
                        buffered_geom = geom.buffer(line_style[1])

                        if buffered_geom is not None:
                            if val == 0 and noaccess == 0:
                                ax.add_geometries(
                                    [buffered_geom],
                                    crs=proj_lat_lon,
                                    linewidth=0,
                                    facecolor='#969696',
                                    edgecolor='none',
                                    zorder=2 + line_style[0]
                                )

                            elif noaccess == 0:
                                ax.add_geometries(
                                    [buffered_geom],
                                    crs=proj_lat_lon,
                                    linewidth=0,
                                    facecolor=str(line_style[2]),
                                    edgecolor='none',
                                    zorder=3 + line_style[0]
                                )
                            else:
                                ax.add_geometries(
                                    [buffered_geom],
                                    crs=proj_lat_lon,
                                    linewidth=0,
                                    facecolor='#252525',
                                    edgecolor='none',
                                    zorder=4 + line_style[0]
                                )

                        else:
                            print("Feature was outside range to plot", iter_)

            x_l = plot_settings['weight_legend']['x_l']
            x_r = plot_settings['weight_legend']['x_r']
            base_y = plot_settings['weight_legend']['base_y']
            y_step = plot_settings['weight_legend']['y_step']
            y_text_nudge = plot_settings['weight_legend']['y_text_nudge']
            x_text_nudge = plot_settings['weight_legend']['x_text_nudge']

            # text above weight legend
            ax.text(
                x_l,
                base_y + y_step - y_text_nudge,
                plot_set[c]['legend_label'],
                horizontalalignment='left',
                transform=proj_lat_lon,
                size=10)

            # weight legend
            divisor = plot_set[c]['divisor']

            for (i, ((nmin, nmax), line_style)) in enumerate(width_by_range.items()):
                y = base_y - (i*y_step)
                line = LineString([(x_l, y), (x_r, y)]).buffer(line_style[1])
                ax.add_geometries(
                    [line],
                    crs=proj_lat_lon,
                    linewidth=0,
                    edgecolor='#969696',
                    facecolor='#969696',
                    zorder=2)
                significance_ndigits = plot_set[c]['significance']
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

            # district labels
            plot_district_labels(ax, config['paths']['data'], highlight_region=region)

            # plot
            title = '{} ({})'.format(region, plot_set[c]['title'])
            print(" * Plotting", title)
            plt.title(title, fontsize=14)
            legend_from_style_spec(ax, styles,loc='upper right')

            # output
            output_file = os.path.join(
                config['paths']['figures'], 'commune_center-{}-{}-failures.png'.format(region.lower().replace(' ', ''), column))
            save_fig(output_file)
            plt.close()


if __name__ == '__main__':
    main()
