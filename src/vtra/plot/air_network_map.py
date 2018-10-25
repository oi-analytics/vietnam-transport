"""Air network map
"""
import os
import sys
from collections import OrderedDict

import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from vtra.utils import *


def main():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'air-map.png')
    air_edge_file_path = os.path.join(
        config['paths']['data'], 'post_processed_networks', 'air_edges.shp')
    air_flow_file_path = os.path.join(config['paths']['output'], 'flow_mapping_combined',
                                   'weighted_flows_national_air_100_percent.csv')
    air_node_file = os.path.join(config['paths']['data'],
                                 'post_processed_networks', 'air_nodes.shp')


    air_edge_file = gpd.read_file(air_edge_file_path,encoding='utf-8')
    air_flow_file = pd.read_csv(air_flow_file_path)
    air_edge_file = pd.merge(air_edge_file,air_flow_file,how='left', on=['edge_id']).fillna(0)

    color_by_type = {'Air route': '#252525', 'Airport': '#d95f0e'}
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'],highlight_region=[])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    geoms = []
    for iter_,record in air_edge_file.iterrows():
        flow = record['max_tons']
        if flow > 0:
            geom = record.geometry
            geoms.append(geom)

    ax.add_geometries(
        geoms,
        crs=proj_lat_lon,
        linewidth=1.5,
        edgecolor='#252525',
        facecolor='none',
        zorder=3)

    # Stations
    xs = []
    ys = []
    for record in shpreader.Reader(air_node_file).records():
        geom = record.geometry
        x = geom.x
        y = geom.y
        xs.append(x)
        ys.append(y)

    ax.scatter(xs, ys, transform=proj_lat_lon, facecolor='#d95f0e', s=12, zorder=5)
    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=line)
        for line, color in color_by_type.items()
    ]
    plt.legend(handles=legend_handles, loc='lower left')
    save_fig(output_file)


if __name__ == '__main__':
    main()
