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

mpl.style.use('ggplot')
mpl.rcParams['font.size'] = 10.
mpl.rcParams['font.family'] = 'tahoma'
mpl.rcParams['axes.labelsize'] = 10.
mpl.rcParams['xtick.labelsize'] = 9.
mpl.rcParams['ytick.labelsize'] = 9.

def plot_ranges(input_data, division_factor,x_label, y_label,plot_title,plot_color,plot_file_path):
    fig, ax = plt.subplots(figsize=(8, 4))
    vals_min_max = list(zip(*list(h for h in input_data.itertuples(index=False))))

    percentlies = 100.0*np.arange(0,len(vals_min_max[0]))/len(vals_min_max[0])
    ax.plot(percentlies,
        1.0*np.array(vals_min_max[0])/division_factor,
        linewidth=0.5,
        color=plot_color
    )
    ax.plot(percentlies,
        1.0*np.array(vals_min_max[1])/division_factor,
        linewidth=0.5,
        color=plot_color
    )
    ax.fill_between(percentlies,
        1.0*np.array(vals_min_max[0])/division_factor,
        1.0*np.array(vals_min_max[1])/division_factor,
        alpha=0.5,
        edgecolor=None,
        facecolor=plot_color
    )

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
        # ax.set_yscale('log')
    
    ax.set_yscale('log')
    # ax.set_xscale('log')
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
    modes = ['Lao Cai','Binh Dinh','Thanh Hoa']
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
    for m in range(len(modes)):
        flow_file_path = os.path.join(config['paths']['output'], 'flow_mapping_combined',
                                       'weighted_edges_commune_center_access_flows_{}_5_tons_100_percent.csv'.format(modes[m].replace(' ','').lower()))
        flow_file = pd.read_csv(flow_file_path).fillna(0)
        flow_file = flow_file.sort_values(['max_netrev'], ascending=True)
        plt_file_path = os.path.join(config['paths']['figures'],'{}-aadf-ranges.png'.format(modes[m].replace(' ','').lower()))
        plot_ranges(flow_file[['min_netrev','max_netrev']],1000, "Percentile rank", 
                "Net Revenue ('000 USD/day)","Range of Net revenue flows on links",modes_colors[m],plt_file_path)

        flow_file_path = os.path.join(config['paths']['output'], 'failure_results','minmax_combined_scenarios',
                                   'single_edge_failures_minmax_{}_5_tons_100_percent_disrupt.csv'.format(modes[m].replace(' ','').lower()))
        flow_file = pd.read_csv(flow_file_path).fillna(0)
        flow_file = flow_file.sort_values(['max_econ_impact'], ascending=True)
        plt_file_path = os.path.join(config['paths']['figures'],'{}-economic_impact-ranges.png'.format(modes[m].replace(' ','').lower()))
        plot_ranges(flow_file[['min_econ_impact','max_econ_impact']],1000, 
                "Percentile rank", 
                "Economic losses ('000 USD/day)","{} - Range of Economic losses due to link failures".format(modes[m]),
                modes_colors[m],plt_file_path)

                    
        flow_file_path = os.path.join(config['paths']['output'], 'hazard_scenarios',
                           'roads_hazard_intersections_{}_risks.csv'.format(modes[m].replace(' ','').lower()))

        fail_scenarios = pd.read_csv(flow_file_path)
        fail_scenarios = pd.merge(fail_scenarios,flow_file,how='left',on=['edge_id']).fillna(0)
        fail_scenarios['min_eael'] = duration*fail_scenarios['max_duration_wt']*fail_scenarios['risk_wt']*fail_scenarios['min_econ_impact']
        fail_scenarios['max_eael'] = duration*fail_scenarios['max_duration_wt']*fail_scenarios['risk_wt']*fail_scenarios['max_econ_impact']

        fail_rcp45 = fail_scenarios[(fail_scenarios['hazard_type'] == 'flooding') & (fail_scenarios['year'] > 2016) & (fail_scenarios['climate_scenario'] == 'rcp 4.5')]
        fail_rcp45.rename(columns={'min_eael':'min_eael_rcp45','max_eael':'max_eael_rcp45'},inplace=True)
        fail_rcp45_min = fail_rcp45.groupby(['edge_id'])['min_eael_rcp45'].min().reset_index()
        fail_rcp45_max = fail_rcp45.groupby(['edge_id'])['max_eael_rcp45'].max().reset_index()
        fail_rcp45 = pd.merge(fail_rcp45_min,fail_rcp45_max,how='left',on=['edge_id']).fillna(0)
        fail_rcp45 = fail_rcp45.sort_values(['max_eael_rcp45'], ascending=True)

        fail_rcp85 = fail_scenarios[(fail_scenarios['hazard_type'] == 'flooding') & (fail_scenarios['year'] > 2016) & (fail_scenarios['climate_scenario'] == 'rcp 8.5')]
        fail_rcp85.rename(columns={'min_eael':'min_eael_rcp85','max_eael':'max_eael_rcp85'},inplace=True)
        fail_rcp85_min = fail_rcp85.groupby(['edge_id'])['min_eael_rcp85'].min().reset_index()
        fail_rcp85_max = fail_rcp85.groupby(['edge_id'])['max_eael_rcp85'].max().reset_index()
        fail_rcp85 = pd.merge(fail_rcp85_min,fail_rcp85_max,how='left',on=['edge_id']).fillna(0)
        fail_rcp85 = fail_rcp85.sort_values(['max_eael_rcp85'], ascending=True)

        fail_cur = fail_scenarios[(fail_scenarios['hazard_type'] == 'flooding') & (fail_scenarios['year'] == 2016)]
        fail_min = fail_cur.groupby(['edge_id'])['min_eael'].min().reset_index()
        fail_max = fail_cur.groupby(['edge_id'])['max_eael'].max().reset_index()
        fail_cur = pd.merge(fail_min,fail_max,how='left',on=['edge_id']).fillna(0)
        fail_cur = fail_cur.sort_values(['max_eael'], ascending=True)
        
        fail_dfs = [fail_cur[['min_eael','max_eael']],fail_rcp45[['min_eael_rcp45','max_eael_rcp45']],fail_rcp85[['min_eael_rcp85','max_eael_rcp85']]]

        plt_file_path = os.path.join(config['paths']['figures'],'{}-flood-eael-ranges.png'.format(modes[m].replace(' ','').lower()))
        plot_many_ranges(fail_dfs,1000, "Percentile rank", 
                "EAEL ('000 USD)",
                "{} - Range of Expected Annual economic losses to link failures".format(modes[m]),
                flood_colors,flood_labels,plt_file_path)

        flow_file_path = os.path.join(config['paths']['output'], 'adaptation_results',
                               'output_adaptation_{}_10_days_max_disruption_fixed_parameters.csv'.format(modes[m].replace(' ','').lower()))
        fail_scenarios = pd.read_csv(flow_file_path)
        fail_scenarios = fail_scenarios[fail_scenarios['max_econ_impact'] > 0]
        
        for cols in ['min_ini_adap_cost','max_ini_adap_cost']:
            fail_scenarios[cols] = fail_scenarios[cols].apply(lambda x: np.max(np.array(ast.literal_eval(x))))
    
        # fail_scenarios = fail_scenarios.groupby(['edge_id'])[adapt_cols].max().reset_index()
        
        for c in range(len(adapt_groups)):
            cols = adapt_groups[c]
            new_cols = ['{}_rcp45'.format(cols[0]),'{}_rcp45'.format(cols[1])]
            fail_rcp45 = fail_scenarios[(fail_scenarios['hazard_type'] == 'flooding') & (fail_scenarios['year'] > 2016) & (fail_scenarios['climate_scenario'] == 'rcp 4.5')]
            fail_rcp45 = fail_rcp45.groupby(['edge_id'])[cols].max().reset_index()
            fail_rcp45.rename(columns={cols[0]:new_cols[0],cols[1]:new_cols[1]},inplace=True)
            fail_rcp45_min = fail_rcp45.groupby(['edge_id'])[new_cols[0]].min().reset_index()
            fail_rcp45_max = fail_rcp45.groupby(['edge_id'])[new_cols[1]].max().reset_index()
            fail_rcp45 = pd.merge(fail_rcp45_min,fail_rcp45_max,how='left',on=['edge_id']).fillna(0)
            fail_rcp45 = fail_rcp45.sort_values([new_cols[1]], ascending=True)
            fail_rcp45 = fail_rcp45[new_cols]

            new_cols = ['{}_rcp85'.format(cols[0]),'{}_rcp85'.format(cols[1])]
            fail_rcp85 = fail_scenarios[(fail_scenarios['hazard_type'] == 'flooding') & (fail_scenarios['year'] > 2016) & (fail_scenarios['climate_scenario'] == 'rcp 4.5')]
            fail_rcp85 = fail_rcp85.groupby(['edge_id'])[cols].max().reset_index()
            fail_rcp85.rename(columns={cols[0]:new_cols[0],cols[1]:new_cols[1]},inplace=True)
            fail_rcp85_min = fail_rcp85.groupby(['edge_id'])[new_cols[0]].min().reset_index()
            fail_rcp85_max = fail_rcp85.groupby(['edge_id'])[new_cols[1]].max().reset_index()
            fail_rcp85 = pd.merge(fail_rcp85_min,fail_rcp85_max,how='left',on=['edge_id']).fillna(0)
            fail_rcp85 = fail_rcp85.sort_values([new_cols[1]], ascending=True)
            fail_rcp85 = fail_rcp85[new_cols]

            fail_cur = fail_scenarios[(fail_scenarios['hazard_type'] == 'flooding') & (fail_scenarios['year'] == 2016)]
            fail_cur = fail_cur.groupby(['edge_id'])[cols].max().reset_index()
            fail_min = fail_cur.groupby(['edge_id'])[cols[0]].min().reset_index()
            fail_max = fail_cur.groupby(['edge_id'])[cols[1]].max().reset_index()
            fail_cur = pd.merge(fail_min,fail_max,how='left',on=['edge_id']).fillna(0)
            fail_cur = fail_cur.sort_values([cols[1]], ascending=True)
            fail_cur = fail_cur[cols]
        
            fail_dfs = [fail_cur,fail_rcp45,fail_rcp85]
            plt_file_path = os.path.join(config['paths']['figures'],'{}-{}-flooding-ranges-fixed-parameters.png'.format(modes[m].replace(' ','').lower(),adapt_names[c]))
            plot_many_ranges(fail_dfs,adapt_divisor[c], 
                    "Percentile rank", 
                    "{} ({})".format(adapt_labels[c],adapt_units[c]),"{} - Range of {} of adaptation".format(modes[m],adapt_labels[c]),
                    flood_colors,flood_labels,plt_file_path)

        if modes[m] == 'Lao Cai':
            for c in range(len(adapt_groups)):
                cols = adapt_groups[c]
                new_cols = ['{}_rcp45'.format(cols[0]),'{}_rcp45'.format(cols[1])]
                fail_rcp45 = fail_scenarios[(fail_scenarios['hazard_type'] == 'landslide') & (fail_scenarios['year'] == 2050) & (fail_scenarios['climate_scenario'] == 'rcp 4.5')]
                fail_rcp45 = fail_rcp45.groupby(['edge_id'])[cols].max().reset_index()
                fail_rcp45.rename(columns={cols[0]:new_cols[0],cols[1]:new_cols[1]},inplace=True)
                fail_rcp45_min = fail_rcp45.groupby(['edge_id'])[new_cols[0]].min().reset_index()
                fail_rcp45_max = fail_rcp45.groupby(['edge_id'])[new_cols[1]].max().reset_index()
                fail_rcp45 = pd.merge(fail_rcp45_min,fail_rcp45_max,how='left',on=['edge_id']).fillna(0)
                fail_rcp45 = fail_rcp45.sort_values([new_cols[1]], ascending=True)
                fail_rcp45 = fail_rcp45[new_cols]

                new_cols = ['{}_rcp85'.format(cols[0]),'{}_rcp85'.format(cols[1])]
                fail_rcp85 = fail_scenarios[(fail_scenarios['hazard_type'] == 'landslide') & (fail_scenarios['year'] == 2050) & (fail_scenarios['climate_scenario'] == 'rcp 4.5')]
                fail_rcp85 = fail_rcp85.groupby(['edge_id'])[cols].max().reset_index()
                fail_rcp85.rename(columns={cols[0]:new_cols[0],cols[1]:new_cols[1]},inplace=True)
                fail_rcp85_min = fail_rcp85.groupby(['edge_id'])[new_cols[0]].min().reset_index()
                fail_rcp85_max = fail_rcp85.groupby(['edge_id'])[new_cols[1]].max().reset_index()
                fail_rcp85 = pd.merge(fail_rcp85_min,fail_rcp85_max,how='left',on=['edge_id']).fillna(0)
                fail_rcp85 = fail_rcp85.sort_values([new_cols[1]], ascending=True)
                fail_rcp85 = fail_rcp85[new_cols]

                fail_cur = fail_scenarios[(fail_scenarios['hazard_type'] == 'landslide') & (fail_scenarios['year'] == 2016)]
                fail_cur = fail_cur.groupby(['edge_id'])[cols].max().reset_index()
                fail_min = fail_cur.groupby(['edge_id'])[cols[0]].min().reset_index()
                fail_max = fail_cur.groupby(['edge_id'])[cols[1]].max().reset_index()
                fail_cur = pd.merge(fail_min,fail_max,how='left',on=['edge_id']).fillna(0)
                fail_cur = fail_cur.sort_values([cols[1]], ascending=True)
                fail_cur = fail_cur[cols]
            
                fail_dfs = [fail_cur,fail_rcp45,fail_rcp85]
                plt_file_path = os.path.join(config['paths']['figures'],'{}-{}-landslide-ranges-fixed-parameters.png'.format(modes[m].replace(' ','').lower(),adapt_names[c]))
                plot_many_ranges(fail_dfs,adapt_divisor[c], 
                        "Percentile rank", 
                        "{} ({})".format(adapt_labels[c],adapt_units[c]),"{} - Range of {} of adaptation".format(modes[m],adapt_labels[c]),
                        flood_colors,flood_labels,plt_file_path)
if __name__ == '__main__':
    main()