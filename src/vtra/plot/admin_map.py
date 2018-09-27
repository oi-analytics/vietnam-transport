"""Admin map
"""
import os
import sys

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt


from vtra.utils import *

def plot_admin_map():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'admin-map.png')
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    save_fig(output_file)

def plot_admin_map_with_regions_highlighted():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'admin-map-regions-highlighted.png')
    ax = get_axes()
    plot_basemap(ax, config['paths']['data'], highlight_region = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    save_fig(output_file)

if __name__ == '__main__':
    plot_admin_map()
    plot_admin_map_with_regions_highlighted()
