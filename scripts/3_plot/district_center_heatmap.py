"""Road network flows
"""
import os
import sys

from collections import OrderedDict

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
            'column': 'min_netrev',
            'title': 'Min Net Rev',
            'legend_label': "MNR ('000 USD/day)",
            'divisor': 1000,
            'significance': 0
        },
        {
            'column': 'max_netrev',
            'title': 'Max Net Rev',
            'legend_label': "MNR ('000 USD/day)",
            'divisor': 1000,
            'significance': 0
        },
        {
            'column': 'min_tons',
            'title': 'Min Tons',
            'legend_label': "MT (tons/day)",
            'divisor': 1,
            'significance': 0
        },
        {
            'column': 'max_tons',
            'title': 'Max Tons',
            'legend_label': "MT (tons/day)",
            'divisor': 1,
            'significance': 0
        }
    ]

    for region in regions:

        region_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_commune_center_flows_' + region.lower().replace(' ', '') + '_5_tons.shp')
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
            weights = [
                record.attributes[column]
                for record in shpreader.Reader(region_file).records()
            ]
            max_weight = max(weights)
            width_by_range = generate_weight_bins_with_colour_gradient(weights, width_step=0.001)

            road_geoms_by_category = {
                region: []
            }

            styles = OrderedDict([
                (region,  Style(color='#ba0f03', zindex=6, label=region))
            ])

            for record in shpreader.Reader(region_file).records():
                cat = region
                geom = record.geometry

                val = record.attributes[column]

                buffered_geom = None
                for (nmin, nmax), line_style in width_by_range.items():
                    if nmin <= val and val < nmax:
                        buffered_geom = geom.buffer(line_style[1])

                        if buffered_geom is not None:
                            ax.add_geometries(
                                [buffered_geom],
                                crs=proj_lat_lon,
                                linewidth=0,
                                facecolor=str(line_style[2]),
                                edgecolor='none',
                                zorder=3 + line_style[0]
                            )
                        else:
                            print("Feature was outside range to plot", record.attributes)

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
                    edgecolor=str(line_style[2]),
                    facecolor=str(line_style[2]),
                    zorder=2)
                significance_ndigits = plot_set[c]['significance']
                if nmin == max_weight:
                    value_template = '>{:.' + str(significance_ndigits) + 'f}'
                    label = value_template.format(round(max_weight/divisor, significance_ndigits))
                else:
                    value_template = '{:.' + str(significance_ndigits) + 'f}-{:.' + str(significance_ndigits) + 'f}'
                    label = value_template.format(round(nmin/divisor, significance_ndigits), round(nmax/divisor, significance_ndigits))
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
            plt.title(title, fontsize = 14)

            # output
            output_file = os.path.join(config['paths']['figures'], 'commune_center-{}-{}.png'.format(region.lower().replace(' ', ''), column))
            save_fig(output_file)
            plt.close()

if __name__ == '__main__':
    main()
