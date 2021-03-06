"""Rail network routes map
"""
import csv
import os
import sys
from collections import OrderedDict, defaultdict
from pprint import pprint

import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from vtra.utils import *


def main():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'rail-map-routes.png')
    rails_file = os.path.join(
        config['paths']['data'], 'post_processed_networks', 'rail_edges.shp')
    rail_descriptions = os.path.join(
        config['paths']['output'], 'network_stats', 'national_rail_routes.xlsx')

    ax = get_axes()
    plot_basemap(ax, config['paths']['data'],highlight_region=[])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    rail_geoms_by_category = defaultdict(list)

    for record in shpreader.Reader(rails_file).records():
        cat = record.attributes['name']
        geom = record.geometry
        rail_geoms_by_category[cat].append(geom)

    colours = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
               '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
               '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
               '#17becf', '#9edae5']

    styles = OrderedDict([])

    reader = pd.read_excel(rail_descriptions,encoding='utf-8')
    for idx, row in reader.iterrows():
        styles.update({row[0]: Style(color=colours[idx], zindex=4, label=row[1])})

    for cat, geoms in rail_geoms_by_category.items():
        cat_style = styles[cat]
        ax.add_geometries(
            geoms,
            crs=proj_lat_lon,
            linewidth=3,
            edgecolor=cat_style.color,
            facecolor='none',
            zorder=cat_style.zindex
        )

    legend_from_style_spec(ax, styles, loc='center left')
    save_fig(output_file)


if __name__ == '__main__':
    main()
