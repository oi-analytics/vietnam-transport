"""Inland network map
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
    output_file = os.path.join(config['paths']['figures'], 'inland-map.png')
    inland_edge_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_flows_national_inland.shp')
    inland_node_file = os.path.join(config['paths']['data'], 'Waterways', 'waterways', 'ports_nodes.shp')

    remove_routes_ids = [
        ('watern_149', 'watern_429'), 
        ('watern_429', 'watern_520'), 
        ('watern_700', 'watern_520'),
        ('watern_210', 'watern_700'),
        ('watern_209', 'watern_210'),
        ('watern_1057', 'watern_1050'),
        ('watern_1050', 'watern_1051'),
        ('watern_1051', 'watern_183'),
        ('watern_183', 'watern_354'),
        ('watern_176', 'watern_354'),
    ]

    color_by_type = {'Inland route': '#0689d7','Inland port': '#d95f0e'}
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    for record in shpreader.Reader(inland_edge_file).records():
        flow = record.attributes['max_tons']
        edge_id = (record.attributes['from_node'], record.attributes['to_node'])
        if flow > 0 and (edge_id not in remove_routes_ids):
            geom = record.geometry
            ax.add_geometries(
                geom,
                crs=proj_lat_lon,
                linewidth=1.5,
                edgecolor='#0689d7',
                facecolor='none',
                zorder=3
            )

    # Stations
    xs = []
    ys = []
    for record in shpreader.Reader(inland_node_file).records():
        port_type = record.attributes['PORT_TYPE']
        if port_type == 'inland':
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
