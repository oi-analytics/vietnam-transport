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

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.utils import *

def main():
	config = load_config()
	crop_cols = ['rice','cash','cass','teas','maiz','rubb','swpo','acof','rcof','pepp']
	title_cols = ['Rice','Cashew','Cassava','Teas','Maize','Rubber','Sweet Potatoes','Coffee Arabica','Coffee Robusta','Pepper']
	crop_file_path = os.path.join(config['paths']['data'], 'Agriculture_crops', 'crop_data')
	for file in os.listdir(crop_file_path):
		if file.endswith('.tif'):
			crop_file = os.path.join(crop_file_path,file)
			crop_name = [cr for cr in crop_cols if cr in file.lower().strip()]
			if crop_name:
				crop_name = crop_name[0]
				crop_title = title_cols[crop_cols.index(crop_name)]
				ax = get_axes()
				plot_basemap(ax, config['paths']['data'])
				scale_bar(ax, location=(0.8, 0.05))
				plot_basemap_labels(ax, config['paths']['data'])
				proj_lat_lon = ccrs.PlateCarree()


				# Create color map
				colors = plt.get_cmap('YlGn')
				#colors.colors[0] = (1, 1, 1, 0)

				# Read in raster data
				data, lat_lon_extent = get_data(crop_file)
				data[data <= 0] = np.nan
				max_val = np.nanmax(data)
				norm=mpl.colors.Normalize(vmin=0, vmax=max_val)
				
				# Plot population data
				im = ax.imshow(data, extent=lat_lon_extent,transform=proj_lat_lon, cmap=colors,norm =norm, zorder=5)

				# Add colorbar
				cbar = plt.colorbar(im, ax=ax,fraction=0.1, shrink=0.87,pad=0.01, drawedges=False, orientation='horizontal',
									norm=mpl.colors.Normalize(vmin=0, vmax=max_val), ticks=list(np.linspace(0,max_val,3)))
				cbar.set_clim(vmin=0,vmax=max_val)


				cbar.outline.set_color("none")
				cbar.ax.yaxis.set_tick_params(color='black')
				cbar.ax.set_xlabel('Crop Annual Production(tons)',fontsize=12,color='black')
				
				plt.title(crop_title, fontsize = 14)
				output_file = os.path.join(config['paths']['figures'], '{}_production.png'.format(crop_name))
				save_fig(output_file)
				plt.close()

if __name__ == '__main__':
	main()
