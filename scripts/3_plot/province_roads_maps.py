"""Road network flows
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from shapely.geometry import LineString


from vtra.utils import *

def main():
    config = load_config()

    regions = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']

    styles = OrderedDict([
        ('National_Road',  Style(color='#ba0f03', zindex=6, label='National')),
        ('Provincial_Road', Style(color='#e0881f', zindex=5, label='Provincial')),
        ('Local_Road', Style(color='#1f99e0', zindex=4, label='Local')),
        ('Other', Style(color='#777777', zindex=3, label='Other')),
    ])

    for region in regions:

        region_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'weighted_edges_commune_center_flows_' + region.lower().replace(' ', '') + '_5_tons.shp')
        plot_settings = get_region_plot_settings(region)

        ax = get_axes(plot_settings['bbox'], figsize=plot_settings['figure_size'])

        if region == 'Binh Dinh':
            plot_basemap(ax, config['paths']['data'], country_border='none', plot_states=False)
        else:
            plot_basemap(ax, config['paths']['data'], country_border='none', plot_states=True)

        scale_bar(ax, location=(0.8, 0.05), length=plot_settings['scale_legend'])
        proj_lat_lon = ccrs.PlateCarree()

        road_geoms_by_category = {
            'National_Road': [],
            'Provincial_Road': [],
            'Local_Road': [],
            'Other': []
        }
        road_geoms_by_category_keys = list(road_geoms_by_category.keys())

        for record in shpreader.Reader(region_file).records():
            geom = record.geometry
            val = record.attributes['level']
            if val > 3 or val < 0:
                val = 3

            road_geoms_by_category[road_geoms_by_category_keys[val]].append(geom)


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

        # output
        output_file = os.path.join(config['paths']['figures'], 'province_roads-{}.png'.format(region.lower().replace(' ', '')))
        save_fig(output_file)
        plt.close()

if __name__ == '__main__':
    main()
