"""Rail network flows map
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.utils import *

def main():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'rail_flow-map.png')
    water_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'railnetworkedges_flows.shp')

    color_by_type = {'Rail Line': '#006d2c'}
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    column = 'total_flow'
    column_label_divisors = {
        "total_flow": 1000,
        "cost": 1000000,
        "centr": 1,
    }
    legend_label = "AADF ('000 tons/day)"
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
            facecolor='#006d2c',
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

    save_fig(output_file)

if __name__ == '__main__':
    main()
