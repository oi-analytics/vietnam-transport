"""Coastal network map
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
    output_file = os.path.join(config['paths']['figures'], 'coastal-map.png')
    coastal_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_flows_national_coastal.shp')
    coastal_node_file = os.path.join(config['paths']['data'], 'Waterways', 'waterways', 'ports_nodes.shp')

    color_by_type = {'Coastal route': '#045a8d','Domestic hub': '#d95f0e','International hub': '#d90e23'}
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    for record in shpreader.Reader(coastal_edge_file).records():
        flow = record.attributes['max_tons']
        if flow > 0:
            geom = record.geometry
            ax.add_geometries(
                geom,
                crs=proj_lat_lon,
                linewidth=1.5,
                edgecolor='#045a8d',
                facecolor='none',
                zorder=3
            )

    # Stations
    xs_domestic = []
    ys_domestic = []

    xs_international = []
    ys_international = []

    for record in shpreader.Reader(coastal_node_file).records():
        port_type = record.attributes['port_class']
        if port_type == 'class_1':
            geom = record.geometry
            x = geom.x
            y = geom.y
            xs_domestic.append(x)
            ys_domestic.append(y)
        elif port_type == 'class_1A':
            geom = record.geometry
            x = geom.x
            y = geom.y
            xs_international.append(x)
            ys_international.append(y)

    ax.scatter(xs_domestic, ys_domestic,transform=proj_lat_lon, facecolor='#d95f0e', s=12, zorder=5)
    ax.scatter(xs_international, ys_international,transform=proj_lat_lon, facecolor='#d90e23', s=12, zorder=5)

    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=line)
        for line, color in color_by_type.items()
]
    plt.legend(handles=legend_handles,loc = 'lower left')
    save_fig(output_file)

if __name__ == '__main__':
    main()
