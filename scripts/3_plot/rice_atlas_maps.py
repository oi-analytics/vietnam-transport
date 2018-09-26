"""Generate crop maps for Vietnam
"""
# pylint: disable=C0103
import os
import sys

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl

mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 11.

#mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 14.
mpl.rcParams['xtick.labelsize'] = 11.
mpl.rcParams['ytick.labelsize'] = 11.
mpl.rcParams['savefig.pad_inches'] = 0.05


from vtra.utils import *

def main():
	config = load_config()
	crop_month_fields = ['P_Jan','P_Feb','P_Mar','P_Apr','P_May','P_Jun','P_Jul','P_Aug','P_Sep','P_Oct','P_Nov','P_Dec']
	title_cols = ['January','February','March','April','May','June','July','August','September','October','November','December']
	for cr in range(len(crop_month_fields)):
		rice_month_file = os.path.join(config['paths']['data'],'rice_atlas_vietnam','rice_production.shp')
		ax = get_axes()
		plot_basemap(ax, config['paths']['data'])
		scale_bar(ax, location=(0.8, 0.05))
		plot_basemap_labels(ax, config['paths']['data'])
		proj = ccrs.PlateCarree()
		for record in shpreader.Reader(rice_month_file).records():
			geom = record.geometry
			region_val = 100.0*record.attributes[crop_month_fields[cr]]/record.attributes['P_total']
			if region_val:
				if region_val > 0 and region_val <= 20:
					color = '#ffffcc' # TODO
					ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor=color,label = '0 to 20')
				elif region_val > 20 and region_val <= 40:
					color = '#c2e699' # TODO
					ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor=color,label = '20 to 40')
				if region_val > 40 and region_val <= 60:
					color = '#78c679' # TODO
					ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor=color,label = '40 to 60')
				elif region_val > 60 and region_val <= 80:
					color = '#31a354' # TODO
					ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor=color,label = '60 to 80')
				elif region_val > 80 and region_val <= 100:
					color = '#31a354' # TODO
					ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor=color,label = '80 to 100')

			else:
				ax.add_geometries([geom], crs=proj, edgecolor='#ffffff', facecolor='#ffffcc',label = '0 to 20')


		colors = ['#ffffcc','#c2e699','#78c679','#31a354','#31a354']
		labels = ['0 to 20','20 to 40','40 to 60','60 to 80','80 to 100']
		# Legend
		legend_handles = []
		for c in range(len(colors)):
			legend_handles.append(mpatches.Patch(color=colors[c], label=labels[c]))

		ax.legend(
			handles=legend_handles,
			title = 'Percentage production',
			loc='center left'
			)


		plt.title(title_cols[cr], fontsize = 14)
		output_file = os.path.join(config['paths']['figures'], 'rice_production_{}.png'.format(title_cols[cr]))
		save_fig(output_file)
		plt.close()

if __name__ == '__main__':
	main()
