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
from vtra.utils import *

def main():
	config = load_config()
	flows_file = os.path.join(config['paths']['data'], 'Results', 'Failure_shapefiles', 'weighted_edges_failures_national_rail_transfer_from_road_10_shift.shp')

	plot_sets = [
		{
			'file_tag': 'commodities',
			'no_access':[0,1],
			'legend_label': "AADF (tons/day)",
			'divisor': 1,
			'columns': ['min_tons','max_tons'],
			'title_cols': ['Total tonnage (min)','Total tonnage (max)']
		},
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
			]

			max_weight = max(weights)
			width_by_range = generate_weight_bins(weights)

			geoms_by_range = {}
			for value_range in width_by_range:
				geoms_by_range[value_range] = []

			for record in shpreader.Reader(flows_file).records():
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
			output_file = os.path.join(config['paths']['figures'], 'rail_failure-map-transfer-road-10-shift-{}-{}.png'.format(plot_set['file_tag'],column))
			save_fig(output_file)
			plt.close()

if __name__ == '__main__':
	main()
