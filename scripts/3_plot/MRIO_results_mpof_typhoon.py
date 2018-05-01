# -*- coding: utf-8 -*-
"""
Created on Tue May  1 11:03:17 2018

@author: cenv0574
"""
import sys
import os
import json

import pandas as pd
import geopandas as gpd
import numpy as np

import cartopy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches


sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from scripts.utils import plot_basemap,plot_basemap_labels_large_region,scale_bar,get_axes,save_fig,set_ax_bg

def get_bounds(ctry_list,world):

    subset = world.loc[world['ADM0_A3'].isin(ctry_list)]

    tot_bounds = list(subset.total_bounds)
    return [tot_bounds[0]-0.2,tot_bounds[2]+0.2,tot_bounds[1]-0.2,tot_bounds[3]+0.2]

def main():
    None

if __name__ == "__main__":
    
    # Define current directory and data directory
    config_path = os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', '..', 'config.json')
    )
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    data_path = config['paths']['data']

    events = ['level13','level14','level15','level16']
    for event in events:
#        if event != 'level13':
#            break
        # load typhoon results
        outcome_event = pd.concat(pd.read_pickle(os.path.join(data_path,'Results','IO_Analysis','losses_{}.pkl'.format(event)))[event],axis=1)
    
        # get list of countries
        ctry_list =  list(np.unique([x[0] for x in list(outcome_event.index)]))
        
        # add China, Myanmar and Indonesia to have a map that makes sense
        ctry_list_plot = [e for e in ctry_list if e not in ('CHN','IND','IDN','BGD','MMR','MYS')]
        
        # load dataframe with countries
        states_filename = os.path.join(    data_path,    'Global_boundaries',    'ne_10m_admin_0_countries_lakes.shp'    )
        world = gpd.read_file(states_filename)
    
        # set figure 
        plt.figure(figsize=(6, 10), dpi=300)
    
        # set projection
        proj_lat_lon = cartopy.crs.PlateCarree()
        ax = plt.axes([0.025, 0.025, 0.95, 0.90],facecolor='#D0E3F4', projection=proj_lat_lon)
    
        # set bounds and extent
        ax.set_extent(get_bounds(ctry_list_plot,world), crs=proj_lat_lon)
        
        # load background 
        set_ax_bg(ax)
    
    #    ax = get_axes(get_bounds(ctry_list_plot,world))
        plot_basemap(ax, config['paths']['data'], focus='VNM', neighbours=ctry_list,plot_regions=False)
    
        # create bin set
        bins = [-1, 0.5, 1, 2.5, 5, 10,20,40]
        
        outcome_countries = outcome_event.groupby(outcome_event.index.get_level_values(0)).sum()
    
        outcome_countries.reset_index(inplace=True)
        subset_viz = world.loc[world['ADM0_A3'].isin(ctry_list),['ADM0_A3','geometry']].merge(outcome_countries,left_on='ADM0_A3',right_on='from_region')
    
        # bin the data for easier plotting with colorscales
        subset_viz['binned_losses'] = pd.cut(subset_viz['total_loss'], bins=bins, labels=[0,1,2,3,4,5,6])
    
        # create cmap
        cmap = cm.get_cmap('Reds', len(bins)) # Colour map (there are many others)
        cmaplist = ['#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d']
        cmap = cmap.from_list('Custom cmap', cmaplist, len(cmaplist))
    
        # plot figure
        subset_viz.plot(ax=ax,column='binned_losses',cmap=cmap,zorder=2,edgecolor='white',lw=1)
        
        # plot labels
        plot_basemap_labels_large_region(ax, config['paths']['data'])
    
        # create legend
        handles = []
        lnames = ['0-0.5 m$','0.5-1 m$','1-2.5 m$','2.5-5 m$','5-10 m$','10-20 m$','20-40 m$'] #,'50-75 km',
        l = 0
        for color in cmaplist:
            handles.append(mpatches.Patch(color=color, label=lnames[l]))
            l += 1
    
        ax.legend(handles=handles,loc=4, prop={'size': 15}, framealpha=0.5) 
    
        # set figure title
        plt.title('Total daily losses in Southeast Asia \n for Typhoon %s' % event,fontweight='bold',fontsize=16)    
    
        output_file = os.path.join(config['paths']['figures'], 'mpof-map-typhoon_{}.png'.format(event))
        save_fig(output_file)

