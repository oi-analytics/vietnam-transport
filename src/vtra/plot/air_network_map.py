"""Air network map
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
    output_file = os.path.join(config['paths']['figures'], 'air-map.png')
    air_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'airport_edges_flows.shp')
    air_node_file = os.path.join(config['paths']['data'], 'Airports', 'airnetwork', 'airportnodes.shp')

    color_by_type = {'Air route': '#252525','Airport': '#d95f0e'}
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    for record in shpreader.Reader(air_edge_file).records():
        flow = record.attributes['total_flow']
        if flow > 0:
            geom = record.geometry
            ax.add_geometries(
                geom,
                crs=proj_lat_lon,
                linewidth=1.5,
                edgecolor='#252525',
                facecolor='none',
                zorder=3
            )

    # Stations
    xs = []
    ys = []
    for record in shpreader.Reader(air_node_file).records():
        geom = record.geometry
        x = geom.x
        y = geom.y
        xs.append(x)
        ys.append(y)

    ax.scatter(xs, ys,transform=proj_lat_lon, facecolor='#d95f0e', s=12, zorder=5)
    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=line)
        for line, color in color_by_type.items()
]
    plt.legend(handles=legend_handles,loc = 'lower left')
    save_fig(output_file)

if __name__ == '__main__':
    main()
