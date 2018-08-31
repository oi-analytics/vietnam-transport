"""Multimodal network map
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.utils import *

def main():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'multimodal-zoom-map.png')
    rail_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_flows_national_rail.shp')
    roads_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_flows_national_road.shp')
    coastal_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_flows_national_coastal.shp')
    multimodal_edge_file = os.path.join(config['paths']['data'], 'Multi', 'multi_edges', 'multimodal_edges.shp')
    multimodal_node_file = os.path.join(config['paths']['data'], 'Multi', 'multi_edges', 'multimodal_nodes.shp')

    color_by_type = {'Rail link and station': '#006d2c','Road link and junction': '#000000', 'Waterway link and port': '#045a8d', 'Multi-Modal link': '#8d3704'}
    ax = get_axes([108.05, 108.35, 16.35, 15.95])
    plot_basemap(ax, config['paths']['data'], country_border='none', plot_states=False)
    scale_bar(ax, location=(0.8, 0.05), length=5)
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

    for record in shpreader.Reader(roads_edge_file).records():
        geom = record.geometry
        ax.add_geometries(
            geom,
            crs=proj_lat_lon,
            linewidth=1.5,
            edgecolor='#000000',
            facecolor='none',
            zorder=3,
            label = 'Roads'
        )

    for record in shpreader.Reader(coastal_edge_file).records():
        geom = record.geometry
        ax.add_geometries(
            geom,
            crs=proj_lat_lon,
            linewidth=1.5,
            edgecolor='#045a8d',
            facecolor='none',
            zorder=3,
            label = 'Waterway'
        )

    for record in shpreader.Reader(multimodal_edge_file).records():
        geom = record.geometry
        ax.add_geometries(
            geom,
            crs=proj_lat_lon,
            linewidth=1.5,
            edgecolor='#8d3704',
            facecolor='none',
            zorder=3,
            label = 'MultiModal'
        )

    geoms = {
        'roads': {
            'xs': [],
            'ys': []
        },
        'railways': {
            'xs': [],
            'ys': []
        },
        'waterways': {
            'xs': [],
            'ys': []
        },
    }

    for record in shpreader.Reader(multimodal_node_file).records():
        mode = record.attributes['mode']
        geom = record.geometry
        x = geom.x
        y = geom.y
        geoms[mode]['xs'].append(x)
        geoms[mode]['ys'].append(y)

    ax.scatter(geoms['roads']['xs'], geoms['roads']['ys'], transform=proj_lat_lon, facecolor='#000000', s=25, zorder=5)
    ax.scatter(geoms['railways']['xs'], geoms['railways']['ys'], transform=proj_lat_lon, facecolor='#006d2c', s=25, zorder=5)
    ax.scatter(geoms['waterways']['xs'], geoms['waterways']['ys'], transform=proj_lat_lon, facecolor='#045a8d', s=25, zorder=5)

    # Legend
    legend_handles = [
        mpatches.Patch(color=color, label=line)
        for line, color in color_by_type.items()
]
    plt.legend(handles=legend_handles,loc = 'lower left')
    save_fig(output_file)

if __name__ == '__main__':
    main()
