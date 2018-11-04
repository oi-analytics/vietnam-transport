"""Road network flows
"""
import sys
import os
import ast
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from vtra.utils import *
from vtra.transport_flow_and_failure_functions import *

mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 10.
mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 10.
mpl.rcParams['xtick.labelsize'] = 9.
mpl.rcParams['ytick.labelsize'] = 9.

def plot_cumsum_ranges(input_data_list, division_factor,x_label, y_label,plot_title,plot_color,plot_file_path):
    fig, ax = plt.subplots(figsize=(8, 4))
    for i in range(len(input_data_list)):
        input_data = input_data_list[i]
        vals_min_max = list(zip(*list(h for h in input_data.itertuples(index=False))))
        ax.plot(1.0*np.array(vals_min_max[0])/division_factor,
            1.0*np.array(vals_min_max[1])/division_factor,
            linewidth=0.5,
            color=plot_color[i]
        )
    # ax.fill_between(percentlies,
    #     1.0*np.array(vals_min_max[0])/division_factor,
    #     1.0*np.array(vals_min_max[1])/division_factor,
    #     alpha=0.5,
    #     edgecolor=None,
    #     facecolor=plot_color
    # )

    # ax.tick_params(axis='x', rotation=45)
    plt.xlabel(x_label, fontweight='bold')
    plt.ylabel(y_label, fontweight='bold')
    plt.title(plot_title)
    
    plt.tight_layout()
    plt.savefig(plot_file_path, dpi=500)
    plt.close()

def plot_many_ranges(input_dfs, division_factor,x_label, y_label,plot_title,plot_color,plot_labels,plot_file_path):
    fig, ax = plt.subplots(figsize=(8, 4))
    
    length = []
    for i in range(len(input_dfs)):
        input_data = input_dfs[i]

        vals_min_max = []
        for a, b in input_data.itertuples(index=False):
            if a < b:
                min_, max_ = a, b
            else:
                min_, max_ = b, a
            vals_min_max.append((min_, max_))

        vals_min_max.sort(key=lambda el: el[1])

        vals_min_max = list(zip(*vals_min_max))

        percentlies = 100.0*np.arange(0,len(vals_min_max[0]))/len(vals_min_max[0])
        length.append(len(vals_min_max[0]))
        ax.plot(percentlies,
            1.0*np.array(vals_min_max[0])/division_factor,
            linewidth=0.5,
            color=plot_color[i]
        )
        ax.plot(percentlies,
            1.0*np.array(vals_min_max[1])/division_factor,
            linewidth=0.5,
            color=plot_color[i]
        )
        ax.fill_between(percentlies,
            1.0*np.array(vals_min_max[0])/division_factor,
            1.0*np.array(vals_min_max[1])/division_factor,
            alpha=0.5,
            edgecolor=None,
            facecolor=plot_color[i],
            label = plot_labels[i]
        )

    length = max(length)
    if 'BCR' in y_label:
        ax.plot(np.arange(0,100),
            np.array([1]*100),
            linewidth=0.5,
            color='red',
            label = 'BCR = 1'
        )
        ax.set_yscale('log')
    
    # ax.tick_params(axis='x', rotation=45)
    ax.legend(loc='upper left')
    plt.xlabel(x_label, fontweight='bold')
    plt.ylabel(y_label, fontweight='bold')
    plt.title(plot_title)
    
    plt.tight_layout()
    plt.savefig(plot_file_path, dpi=500)
    plt.close()

def main():
    config = load_config()
    regions = ['Lao Cai','Binh Dinh','Thanh Hoa','National roads']
    modes_colors = ['#000004','#006d2c','#0689d7','#045a8d']
    flood_colors = ['#252525','#54278f','#08519c']
    flood_labels = ['Current','RCP 4.5','RCP 8.5']
    adapt_cols = ['min_benefit','min_ini_adap_cost','min_tot_adap_cost','min_bc_ratio','max_benefit','max_ini_adap_cost','max_tot_adap_cost','max_bc_ratio']
    adapt_groups = [['min_benefit','max_benefit'],['min_ini_adap_cost','max_ini_adap_cost'],['min_tot_adap_cost','max_tot_adap_cost'],['min_bc_ratio','max_bc_ratio']]
    adapt_names = ['benefit','ini_adap_cost','tot_adap_cost','bc_ratio']
    adapt_labels = ['Benefits','Initial investments','Total investments', 'BCR']
    adapt_units = ['million USD','million USD','million USD','ratio']
    adapt_divisor = [1000000,1000000,1000000,1]
    duration = 10
    for region in regions:
        if region in ['Lao Cai','Binh Dinh','Thanh Hoa']:
            adaptation_df =  pd.read_csv(os.path.join(config['paths']['output'],
                'network_stats',
                '{}_adapt_summary_fixed_parameter.csv'.format(region.replace(' ','').lower())
                )
            )
        else:
            adaptation_df = pd.read_csv(os.path.join(config['paths']['output'],
                'network_stats',
                'national_roads_adapt_summary_fixed_parameters.csv')
            )


        adaptation_df = adaptation_df[adapt_cols]
        for cols in adapt_groups:
            adaptation_df['swap'] = adaptation_df.apply(lambda x: swap_min_max(x,cols[0],cols[1]),axis=1)
            adaptation_df[cols] = adaptation_df['swap'].apply(pd.Series)
            adaptation_df.drop('swap', axis=1, inplace=True)

        adaptation_df = adaptation_df.sort_values(['max_bc_ratio'], ascending=False)
        print (adaptation_df)
        adaptation_df = adaptation_df.cumsum()

        # print (adaptation_df)
        plot_values = [adaptation_df[['min_tot_adap_cost','min_benefit']],adaptation_df[['max_tot_adap_cost','max_benefit']]]
        plt_file_path = os.path.join(config['paths']['figures'],'{}-cost-benefit-returns-ranges.png'.format(region))
        plot_cumsum_ranges(plot_values,1000000, "Cumulative investments (million USD)", 
                "Cumulative benefits (million USD)","Investments vs Benefits returns for {}".format(region),['#000004','#006d2c'],plt_file_path)
if __name__ == '__main__':
    main()