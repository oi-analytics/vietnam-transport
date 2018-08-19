# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 13:52:00 2017

@author: cenv0574
"""

# Import
import pandas as pd 
import geopandas as gpd
import numpy as np
import sys
import os

import cartopy
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config,plot_basemap,plot_basemap_labels,set_ax_bg
from table import io_basic
from model import MRIA_IO as MRIA


if __name__ == '__main__':

    data_path = load_config()['paths']['data']
   
    ''' Specify file path '''
    filepath =  os.path.join(data_path,'input_data','IO_VIETNAM.xlsx')

    '''Create data input'''
    DATA = io_basic('Vietnam',filepath,2010)
    DATA.prep_data()
            
    '''Create model '''    
    MRIA_model = MRIA(DATA.name, DATA.countries,DATA.sectors,DATA.FD_cat)
    MRIA_model.create_sets(FD_SET=['FinDem'])
    MRIA_model.create_alias()

    '''Run model and create some output'''
    output = pd.DataFrame()
 
    '''Specify disruption'''
    disruption = 0.90
    disrupted_ctry =  ['Ha_Noi','Ho_Chi_Minh']
    disrupted_sctr = ['secA', 'secB', 'secC', 'secD', 'secE', 'secF', 'secI', 'secG', 'secH']

    disr = pd.DataFrame(columns = disrupted_ctry,index =disrupted_sctr)
    disr.loc[disrupted_sctr] = disruption
    disr_dict_sup = {} #{(k,r): v for r, kv in disr.iterrows() for k,v in kv.to_dict().items()}

    disrupted_org = ['Ha_Noi']
    disrupted_des =  ['Ho_Chi_Minh','Bac_Ninh','Ha_Nam']
    disrupted_sctr = ['secA','secB','secC', 'secD', 'secE', 'secF', 'secI', 'secG', 'secH']

    disr = pd.DataFrame(columns = disrupted_des,index =disrupted_sctr)
    disr.loc[disrupted_sctr] = np.random.randint(0,10,size=disr.shape)/100
    disr_dict_fd = {(disrupted_org[0],k,r): v for r, kv in disr.iterrows() for k,v in kv.to_dict().items()}

    '''Create model'''
    MRIA_RUN = MRIA(DATA.name,DATA.countries,DATA.sectors,EORA=False,list_fd_cats=['FinDem'])
    
    '''Define sets and alias'''
    # CREATE SETS
    MRIA_RUN.create_sets()
    
    # CREATE ALIAS
    MRIA_RUN.create_alias()
    
    ''' Define tables and parameters'''
    Regmaxcap = 0.98
    
    MRIA_RUN.baseline_data(DATA,disr_dict_sup,disr_dict_fd)
    MRIA_RUN.impact_data(DATA,disr_dict_sup,disr_dict_fd)

    output['x_in'] = pd.Series(MRIA_RUN.X.get_values())*43
#   
    MRIA_RUN.run_impactmodel()
#  
    output['x_out'] = pd.Series(MRIA_RUN.X.get_values())*43
    output['loss'] = (output['x_out'] - output['x_in'])

    prov_path = os.path.join(data_path,'Vietnam_boundaries','boundaries_stats','province_level_stats.shp')
    provinces = gpd.read_file(prov_path)[['name_eng','geometry']]
    provinces.name_eng = provinces.name_eng.apply(lambda x: x.replace(' ','_').replace('-','_'))
    
    prov_impact = output.groupby(level=0,axis=0).sum()
    
    prov_impact = provinces.merge(prov_impact,left_on='name_eng',right_index=True)
    # create bin set
    bins = [-1, 1, 5, 10, 25, 50,75,np.float('inf')]

    # bin the data for easier plotting with colorscales
    prov_impact['binned_loss'] = pd.cut((prov_impact['loss']*-1), bins=bins, labels=[0,1,2,3,4,5,6])
    prov_impact.loc[prov_impact.geometry.area.idxmin(),'binned_no_flood'] = 6

    # create figure        
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(9,7))

    # set projection
    proj_lat_lon = cartopy.crs.PlateCarree()
    ax = plt.axes(facecolor='#D0E3F4', projection=proj_lat_lon)

    # set bounds and extent
    ax.set_extent([102.2, 109.5, 8.5, 23.3], crs=proj_lat_lon)
    
    # load background 
    set_ax_bg(ax)

    plot_basemap(ax, data_path,country_border=None)
    plot_basemap_labels(ax, data_path,province_zoom=False)
    
    # create cmap
    cmap = cm.get_cmap('Reds', len(bins)) # Colour map (there are many others)
    cmaplist = ['#fee5d9','#fcbba1','#fc9272','#fb6a4a','#de2d26','#a50f15','#442c2d']
    cmap = cmap.from_list('Custom cmap', cmaplist, len(cmaplist))

    # plot figure
    prov_impact.plot(ax=ax,column='binned_loss',cmap=cmap,zorder=2)

    # create legend
    handles = []
    lnames = ['0-1 M-USD','1-5 M-USD','5-10 M-USD','10-25 M-USD','25-50 M-USD','50-75 M-USD', '>75 M-USD'] #,'50-75 km',
    l = 0
    for color in cmaplist:
        handles.append(mpatches.Patch(color=color, label=lnames[l]))
        l += 1

    ax.legend(handles=handles,loc=6, prop={'size': 10}) 

#    # set figure title
#    plt.title('Distance to major towns in Thanh Hoa for %s' % rp_names[rp],fontweight='bold',fontsize=16)
#
    # save figure
    figure_out= os.path.join(data_path,'losses_mria.png')
    plt.savefig(figure_out,dpi=600,bbox_inches='tight')