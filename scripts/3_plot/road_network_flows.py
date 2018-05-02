"""Road network flows
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from shapely.geometry import LineString

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.utils import *

def main():
    config = load_config()
    output_file = os.path.join(config['paths']['figures'], 'road_flow-map.png')
    edges_file = os.path.join(config['paths']['data'], 'Roads', 'roads2009', 'roads2009edges.shp')
    flows_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'road2009edges_flows.shp')

    # edge_id => type
    edge_type_by_id = {}
    for record in shpreader.Reader(edges_file).records():
        id_ = record.attributes['edge_id']
        type_ = record.attributes['type']
        edge_type_by_id[id_] = type_

    ax = get_axes()
    plot_basemap(ax, config['paths']['data'])
    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, config['paths']['data'])
    proj_lat_lon = ccrs.PlateCarree()

    # generate weight bins
    column = 'total_flow'
    weights = [
        record.attributes[column]
        for record in shpreader.Reader(flows_file).records()
    ]
    max_weight = max(weights)
    width_by_range = generate_weight_bins(weights)

    road_geoms_by_category = {
        'National_Road': [],
        'Provincial_Road': [],
        'District_Road': [],
        'Other': []
    }

    for record in shpreader.Reader(flows_file).records():
        cat = edge_type_by_id[record.attributes['edge_id']]
        if cat not in road_geoms_by_category:
            cat = 'Other'
        geom = record.geometry

        val = record.attributes[column]

        buffered_geom = None
        for (nmin, nmax), width in width_by_range.items():
            if nmin <= val and val < nmax:
                buffered_geom = geom.buffer(width)

        if buffered_geom is not None:
            road_geoms_by_category[cat].append(buffered_geom)
        else:
            print("Feature was outside range to plot", record.attributes)

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
            linewidth=0,
            facecolor=cat_style.color,
            edgecolor='none',
            zorder=cat_style.zindex
        )

    # Flow weight legend
    column_label_divisors = {
        "total_flow": 1000,
        "cost": 1000000,
        "centr": 1,
    }
    legend_label = "AADF ('000 tons/day)"

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
        line = LineString([(x_l, y), (x_r, y)]).buffer(width)
        ax.add_geometries(
            [line],
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

    legend_from_style_spec(ax, styles)
    save_fig(output_file)



if __name__ == '__main__':
    main()
