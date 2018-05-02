"""Road network map
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
    output_file = os.path.join(config['paths']['figures'], 'road-map.png')
    roads_file = os.path.join(config['paths']['data'], 'Roads', 'roads2009', 'roads2009edges.shp')

    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    road_geoms_by_category = {
        'National_Road': [],
        'Provincial_Road': [],
        'District_Road': [],
        'Other': []
    }

    for record in shpreader.Reader(roads_file).records():
        cat = record.attributes['type']
        if cat not in road_geoms_by_category:
            cat = 'Other'
        geom = record.geometry
        road_geoms_by_category[cat].append(geom)

    styles = OrderedDict([
        ('National_Road',  Style(color='#ba0f03', zindex=6, label='National')),
        ('Provincial_Road', Style(color='#e0881f', zindex=5, label='Provincial')),
        ('District_Road', Style(color='#1f99e0', zindex=4, label='District')),
        ('Other', Style(color='#777777', zindex=3, label='Other')),
    ])

    for cat, geoms in road_geoms_by_category.items():
        cat_style =styles[cat]
        ax.add_geometries(
            geoms,
            crs=proj_lat_lon,
            linewidth=1,
            edgecolor=cat_style.color,
            facecolor='none',
            zorder=cat_style.zindex
        )

    legend_from_style_spec(ax, styles)
    save_fig(output_file)

if __name__ == '__main__':
    main()
