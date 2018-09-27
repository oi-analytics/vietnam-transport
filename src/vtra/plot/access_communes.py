# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 15:42:31 2018

@author: cenv0574
"""

import pandas as pd
import geopandas as gpd
import os
import json
import numpy as np

import cartopy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

from vtra.utils import plot_basemap,plot_basemap_labels,scale_bar,set_ax_bg


mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 13.
mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 14.
mpl.rcParams['xtick.labelsize'] = 13.
mpl.rcParams['ytick.labelsize'] = 13.

def main():

    # =============================================================================
    #     # Define current directory and data directory
    # =============================================================================
    config_path = os.path.realpath(
        os.path.join(os.path.dirname(__file__), '..', '..', 'config.json')
    )
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    data_path = config['paths']['data']
    output_path = config['paths']['output']
    figure_path = config['paths']['figures']

    # =============================================================================
    #     # Load commune data
    # =============================================================================
    commune_affected_path = os.path.join(output_path,'Thanh_Hoa_affected_communes_flood.shp')
    communes_affected = gpd.read_file(commune_affected_path)

    # =============================================================================
    #     # prepare some lists for the final output
    # =============================================================================

    sectors = ['nongnghiep', 'khaikhoang', 'chebien', 'detmay','gogiay', 'sanxuat', 'xaydung', 'thuongmai', 'dichvu']
    sectors_eng = ['agriculture', 'mining', 'processing', 'textile','wood & paper', 'manufacture', 'construction', 'trade', 'service']
    sectors_eng = [x.capitalize() for x in sectors_eng]

    rps = ['no_flood','1_freq','5_freq','10_freq','20_freq']
    rp_names = {'no_flood': 'no flooding','1_freq': 'a 1/100 flood','5_freq' : 'a 1/20 flood',
                '10_freq' : 'a 1/10 flood','20_freq': 'a 1/5 flood'}

    # =============================================================================
    #     # create all figures per return period
    # =============================================================================

    for rp in rps:
        # create map of affected communes in province
        map_distance_per_commune(rp,communes_affected,figure_path,rp_names,data_path)

        if rp != 'no_flood':

            # create pie chart of affected firms and employees in affected communes in province
            create_pie_plots(rp,communes_affected,sectors,sectors_eng,rp_names,figure_path)

            # create bar plot of affected firms and employees in province
            create_bar_plots(rp,communes_affected,sectors,sectors_eng,rp_names,figure_path)


def map_distance_per_commune(rp,communes_affected,figure_path,rp_names,data_path):

    # create bin set
    bins = [-1, 1, 5, 10, 25, 50,75,np.float('inf')]

    # bin the data for easier plotting with colorscales
    communes_affected['binned_{}'.format(rp)] = pd.cut(communes_affected[rp], bins=bins, labels=[0,1,2,3,4,5,6])
    communes_affected.loc[communes_affected.geometry.area.idxmin(),'binned_no_flood'] = 6

    # create figure
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(8,6))

    # set projection
    proj_lat_lon = cartopy.crs.PlateCarree()
    ax = plt.axes([0.0,0.0,1.0, 1.0] ,facecolor='#D0E3F4', projection=proj_lat_lon)

    # set bounds and extent
    tot_bounds = list(communes_affected.total_bounds)
    ax.set_extent([tot_bounds[0]-0.1,tot_bounds[2]+0.1,tot_bounds[1]-0.1,tot_bounds[3]+0.1] , crs=proj_lat_lon)

    # load background
#    world = gpd.read_file(os.path.join(data_path,'Vietnam_boundaries','who_boundaries','who_provinces.shp'))
#    world.plot(ax=ax,color='#FEF9E0',lw=0.3,edgecolor='k')
    set_ax_bg(ax)

    plot_basemap(ax, data_path,country_border=None)
#    scale_bar(ax, location=(0.8, 0.05))
    plot_basemap_labels(ax, data_path,province_zoom=True)

    # create cmap
    cmap = cm.get_cmap('Reds', len(bins)) # Colour map (there are many others)
    cmaplist = ['#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15','#442c2d']
    cmap = cmap.from_list('Custom cmap', cmaplist, len(cmaplist))

    # plot figure
    communes_affected.plot(ax=ax,column='binned_{}'.format(rp),cmap=cmap,zorder=2)


    # create legend
    handles = []
    lnames = ['0-1 km','1-5 km','5-10 km','10-25 km','25-50 km','50-75 km', 'No Access'] #,'50-75 km',
    l = 0
    for color in cmaplist:
        handles.append(mpatches.Patch(color=color, label=lnames[l]))
        l += 1

    ax.legend(handles=handles,loc=3, prop={'size': 13})

    # set figure title
    plt.title('Distance to major towns in Thanh Hoa for %s' % rp_names[rp],fontweight='bold',fontsize=16)

    # save figure
    figure_out= os.path.join(figure_path,'Dist_Major_Towns_Thanh_Hoa_%s.png' % rp)
    plt.savefig(figure_out,dpi=600,bbox_inches='tight')


def create_pie_plots(rp,communes_affected,sectors,sectors_eng,rp_names,figure_path):

    # estimate some summary statistics for final plotting. Amount of firms per sector in this case
    only_communes_affected = communes_affected.loc[communes_affected['no_flood'] != communes_affected[rp]]
    only_communes_affected = only_communes_affected.copy()
    only_communes_affected[sectors] = only_communes_affected[sectors].multiply(only_communes_affected['nfirm'],axis='index')

    firms_affected = only_communes_affected[sectors].sum()

    # and get english sector names for plotting
    firms_affected.index = sectors_eng

    # =============================================================================
    #         # create figure for total number of firms affected
    # =============================================================================
    fig2,ax2 = plt.subplots(nrows=1, ncols=1,figsize=(6,6))

    # create cmap
    cmap = cm.get_cmap('Reds') # Colour map (there are many others)
    cmaplist = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']
    cmap = cmap.from_list('Custom cmap', cmaplist, len(cmaplist))

    # plot figure
    firms_affected.plot(kind='pie',ax=ax2, fontsize=12,cmap=cmap,startangle=90, pctdistance=0.85)

    ax2.axis('equal')
    plt.axis('off')

    # set figure title
    plt.title('Relative distribution of industries affected in \n affected communes Thanh Hoa for %s' % rp_names[rp],fontweight='bold',fontsize=17)

    # save figure
    figure_out= os.path.join(figure_path,'Share_industries_affected_communes_Thanh_Hoa_%s.png' % rp)
    plt.savefig(figure_out,dpi=600,bbox_inches='tight')


    # =============================================================================
    #         # create figure for total number of employees affected
    # =============================================================================

    # estimate employees per sector, furthe rmultiplcation of number of firms per sector
    only_communes_affected[sectors] = only_communes_affected[sectors].multiply(only_communes_affected['labor'],axis='index')
    firms_affected = only_communes_affected[sectors].sum()
    firms_affected.index = sectors_eng

    # set figure and axis
    fig2,ax2 = plt.subplots(nrows=1, ncols=1,figsize=(6,6))

    # plot figure
    firms_affected.plot(kind='pie',ax=ax2, fontsize=12,cmap=cmap,startangle=90, pctdistance=0.85)

    ax2.axis('equal')
    plt.axis('off')

    # set figure title
    plt.title('Relative distribution of employees affected in \n affected communes Thanh Hoa for %s' % rp_names[rp],fontweight='bold',fontsize=17)

    # save figure
    figure_out= os.path.join(figure_path,'Share_employees_affected_communes_Thanh_Hoa_%s.png' % rp)
    plt.savefig(figure_out,dpi=600,bbox_inches='tight')


def create_bar_plots(rp,communes_affected,sectors,sectors_eng,rp_names,figure_path):

    # =============================================================================
    #         # create figure relative impact industries for Thanh Hao
    # =============================================================================

    # create summary statistics for relative stats for total province
    only_communes_affected = communes_affected.loc[communes_affected['no_flood'] != communes_affected[rp]]
    firms_affected = only_communes_affected[sectors].multiply(only_communes_affected['nfirm'],axis='index').sum()


    firms_affected_TH = communes_affected[sectors].multiply(communes_affected['nfirm'],axis='index').sum()
    firms_affected_TH = (firms_affected/firms_affected_TH)*100

    firms_affected_TH.index = sectors_eng

    # set figure and axis
    fig3,ax3 = plt.subplots(nrows=1, ncols=1,figsize=(6,6))

    # set colormap
    cmaplist = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']

    # plot figure
    firms_affected_TH.plot(kind='bar',ax=ax3, fontsize=12,color=cmaplist)

    # set titles
    ax3.set_ylabel("Percentage affected", fontweight='bold',fontsize=15)
    ax3.set_xlabel("Industrial sector", fontweight='bold',fontsize=15)

    plt.title('Relative share of firms affected \n in Thanh Hoa for %s' % rp_names[rp],fontweight='bold',fontsize=17)

    # save figure
    figure_out= os.path.join(figure_path,'Share_firms_affected_Thanh_Hoa_%s.png' % rp)
    plt.savefig(figure_out,dpi=600,bbox_inches='tight')



if __name__ == "__main__":
    main()
