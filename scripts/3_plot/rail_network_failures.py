"""Rail network flows map
"""
import os
import sys

from collections import OrderedDict

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.utils import *

def main():
	config = load_config()
	flows_file = os.path.join(config['paths']['data'], 'Results', 'Failure_shapefiles', 'weighted_edges_failures_national_rail_2.shp')
	
	plot_sets = [
		{
			'file_tag': 'reroute',
			'no_access': [-1,0],
			'legend_label': "(million USD/day)",
			'divisor': 1000000,
			'columns': ['min_tr_los','max_tr_los'],
			'title_cols': ['Rerouting costs (min)','Rerouting costs (max)']
		},
		{
			'file_tag': 'commodities',
			'no_access':[-1,1],
			'legend_label': "AADF (tons/day)",
			'divisor': 1,
			'columns': ['min_tons','max_tons'],
			'title_cols': ['Total tonnage (min)','Total tonnage (max)']
		},
		{
			'file_tag': 'economic',
			'no_access':[-1,1],
			'legend_label': "(million USD/day)",
			'divisor': 1000000,
			'columns': ['min_econ_l','max_econ_l'],
			'title_cols': ['Economic losses (min)','Economic losses (max)']
		},
		{
			'file_tag': 'total',
			'no_access':[0,1],
			'legend_label': "(million USD/day)",
			'divisor': 1000000,
			'columns': ['min_loss','max_loss'],
			'title_cols': ['Total economic impact (min)','Total economic impact (max)']
		}

	]

	color = '#006d2c'
	color_by_type = {'Rail Line': color}
	
	for plot_set in plot_sets:
		for c in range(len(plot_set['columns'])):
			ax = get_axes()
			plot_basemap(ax, config['paths']['data'], highlight_region = [])
			scale_bar(ax, location=(0.8, 0.05))
			plot_basemap_labels(ax, config['paths']['data'])
			proj_lat_lon = ccrs.PlateCarree()
			
			column = plot_set['columns'][c]
			weights = [
				record.attributes[column]
				for record in shpreader.Reader(flows_file).records()
				if int(record.attributes['no_access']) in plot_set['no_access']
			]
			
			max_weight = max(weights)
			width_by_range = generate_weight_bins(weights)
		
			geoms_by_range = {}
			for value_range in width_by_range:
				geoms_by_range[value_range] = []

			for record in [rec for rec in shpreader.Reader(flows_file).records() if int(rec.attributes['no_access']) in plot_set['no_access']]: 
				val = record.attributes[column]
				geom = record.geometry
				for nmin, nmax in geoms_by_range:
					if nmin <= val and val < nmax:
						geoms_by_range[(nmin, nmax)].append(geom)

						# plot
			for range_, width in width_by_range.items():
				ax.add_geometries(
					[geom.buffer(width) for geom in geoms_by_range[range_]],
					crs=proj_lat_lon,
					edgecolor='none',
					facecolor=color,
					zorder=2)
	
			x_l = 102.3
			x_r = x_l + 0.4
			base_y = 14
			y_step = 0.4
			y_text_nudge = 0.1
			x_text_nudge = 0.1

			ax.text(
				x_l,
				base_y + y_step - y_text_nudge,
				plot_set['legend_label'],
				horizontalalignment='left',
				transform=proj_lat_lon,
				size=10)

			divisor = plot_set['divisor']
			for (i, ((nmin, nmax), width)) in enumerate(width_by_range.items()):
				y = base_y - (i*y_step)
				line = LineString([(x_l, y), (x_r, y)])
				ax.add_geometries(
					[line.buffer(width)],
					crs=proj_lat_lon,
					linewidth=0,
					edgecolor=color,
					facecolor=color,
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
			
			plt.title(plot_set['title_cols'][c], fontsize = 14)
			output_file = os.path.join(config['paths']['figures'], 'rail_failure-map-{}-{}.png'.format(plot_set['file_tag'],column))
			save_fig(output_file)
			plt.close()

if __name__ == '__main__':
	main()
