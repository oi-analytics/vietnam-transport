"""Rail network map
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt


from vtra.utils import *

def main():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'rail-map.png')
    rail_edge_file = os.path.join(config['paths']['data'], 'Railways', 'national_rail', 'railway_2009_edges.shp')
    rail_node_file = os.path.join(config['paths']['data'], 'Railways', 'national_rail', 'railwaystations_2009_nodes.shp')

    color_by_type = {'Rail line': '#006d2c','Rail stop': '#000000'}
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    for record in shpreader.Reader(rail_edge_file).records():
        geom = record.geometry
        ax.add_geometries(
            geom,
            crs=proj_lat_lon,
            linewidth=1.5,
            edgecolor='#006d2c',
            facecolor='none',
            zorder=3,
            label = 'Rail line'
        )

    # Stations
    xs = []
    ys = []
    for record in shpreader.Reader(rail_node_file).records():
        node_type = record.attributes['NAME']
        if node_type != '0':
            geom = record.geometry
            x = geom.x
            y = geom.y
            xs.append(x)
            ys.append(y)
            name = record.attributes['NAME']

    ax.scatter(xs, ys,transform=proj_lat_lon, facecolor='#000000', s=4, zorder=5, label = 'Rail station')
    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=line)
        for line, color in color_by_type.items()
]
    plt.legend(handles=legend_handles,loc = 'lower left')
    save_fig(output_file)

if __name__ == '__main__':
    main()
