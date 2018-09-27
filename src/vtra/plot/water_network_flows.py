"""Water network flows
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString


from vtra.utils import *

def main():
    config = load_config()
    water_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'wateredges_flows.shp')

    color = '#045a8d'
    color_by_type = {'Waterway route': color}

    crop_cols = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']
    ind_cols = ['sugar','wood','steel','constructi','cement','fertilizer','coal','petroluem','manufactur','fishery','meat']

    columns = crop_cols + ind_cols + ['total_flow']
    column_label_divisors = {c: 1000 for c in columns}

    legend_label = "AADF ('000 tons/day)"
    title_cols = ['Rice','Cashew','Cassava','Teas','Maize','Rubber','Sweet Potatoes','Coffee Arabica','Coffee Robusta','Pepper',
                'Sugar','Wood','Steel','Construction materials','Cement','Fertilizer','Coal','Petroleum',
                'Manufacturing','Fishery','Meat','Total tonnage']

    for c in range(len(columns)):
        ax = get_axes()
        plot_basemap(ax, config['paths']['data'])
        scale_bar(ax, location=(0.8, 0.05))
        plot_basemap_labels(ax, config['paths']['data'])
        proj_lat_lon = ccrs.PlateCarree()

        column = columns[c]
        # generate weight bins
        weights = [
            record.attributes[column]
            for record
            in shpreader.Reader(water_edge_file).records()
        ]
        max_weight = max(weights)
        width_by_range = generate_weight_bins(weights)

        geoms_by_range = {}
        for value_range in width_by_range:
            geoms_by_range[value_range] = []

        for record in shpreader.Reader(water_edge_file).records():
            val = record.attributes[column]
            geom = record.geometry
            for nmin, nmax in geoms_by_range:
                if nmin <= val and val < nmax:
                    geoms_by_range[(nmin, nmax)].append(geom)

                    # plot
        for range_, width in width_by_range.items():
            ax.add_geometries(
                [geom.buffer(width) for geom in geoms_by_range[range_]],
                crs=proj_lat_lon,
                edgecolor='none',
                facecolor=color,
                zorder=2)

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
        plt.title(title_cols[c], fontsize = 14)
        output_file = os.path.join(config['paths']['figures'], 'water_flow-map-{}.png'.format(column))
        save_fig(output_file)
        plt.close()

if __name__ == '__main__':
    main()
