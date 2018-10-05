# -*- coding: utf-8 -*-
"""Assign commodity flows on the road network
"""
import ast
import copy
import csv
import itertools
import math
import operator
import os
import sys
from collections import Counter

import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
import psycopg2
from sqlalchemy import create_engine
from vtra.transport_network_creation import *
from vtra.utils import *


def main():
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    types = ['min', 'max']
    path_types = ['min_edge_path', 'max_edge_path']
    tons_types = ['min_tons', 'max_tons']
    dist_types = ['min_distance', 'max_distance']
    time_types = ['min_time', 'max_time']
    cost_types = ['min_gcost', 'max_gcost']
    vehicle_types = ['min_vehicle_nums', 'max_vehicle_nums']
    rice_type = ['min_rice', 'max_rice']
    ind_crop_cols = ['sugar', 'wood', 'steel', 'constructi', 'cement', 'fertilizer', 'coal', 'petroluem', 'manufactur',
                     'fishery', 'meat', 'cash', 'cass', 'teas', 'maiz', 'rubb', 'swpo', 'acof', 'rcof', 'pepp']

    shp_output_path = os.path.join(output_path, 'failure_shapefiles')
    flow_paths_data = os.path.join(output_path, 'flow_mapping_paths',
                                   'national_scale_flow_paths_2.xlsx')
    fail_scenarios_data = os.path.join(
        output_path, 'hazard_scenarios', 'national_scale_hazard_intersections.xlsx')

    """
    Get the modal shares
    """
    # modes_file_paths = [('Roads','national_roads'), ('Railways','national_rail'), ('Waterways','waterways'), ('Waterways','waterways')]
    # modes_file_paths = [('Roads','national_roads'), ('Railways','national_rail')]
    modes_file_paths = [('Roads', 'national_roads'), ('Railways', 'national_rail'),
                        ('Waterways', 'waterways'), ('Waterways', 'waterways'), ('Multi', 'multi_edges')]
    # modes_file_paths = [('Roads','national_roads')]
    modes = ['road', 'rail', 'inland', 'coastal', 'multi']
    veh_wt = [20, 800, 800, 1200]
    usage_factors = [(0, 0), (0, 0), (0.2, 0.25), (0.2, 0.25), (0, 0)]
    speeds = [(0, 0), (40, 60), (9, 20), (9, 20), (0, 0)]
    multi_md_len = 3.0

    md_prop_file = os.path.join(data_path, 'mode_properties', 'mode_costs.xlsx')
    rd_prop_file = os.path.join(data_path, 'mode_properties', 'road_properties.xlsx')

    G_multi_df = []
    for m in range(len(modes_file_paths)):
        mode_data_path = os.path.join(
            data_path, modes_file_paths[m][0], modes_file_paths[m][1])
        for file in os.listdir(mode_data_path):
            try:
                if file.endswith(".shp") and 'edges' in file.lower().strip():
                    edges_in = os.path.join(mode_data_path, file)
            except:
                return ('Network nodes and edge files necessary')

        if modes[m] == 'road':
            G_df = national_road_shapefile_to_dataframe(edges_in, rd_prop_file)
        elif modes[m] == 'multi':
            G_df = multi_modal_shapefile_to_dataframe(
                edges_in, md_prop_file, modes[m], multi_md_len)
        else:
            G_df = network_shapefile_to_dataframe(
                edges_in, md_prop_file, modes[m], speeds[m][0], speeds[m][1])
            if modes[m] in ('inland', 'coastal'):
                G_df['min_time_cost'] = (1 + usage_factors[m][0])*G_df['min_time_cost']
                G_df['max_time_cost'] = (1 + usage_factors[m][1])*G_df['max_time_cost']
                G_df['min_tariff_cost'] = (1 + usage_factors[m][0])*G_df['min_tariff_cost']
                G_df['max_tariff_cost'] = (1 + usage_factors[m][1])*G_df['max_tariff_cost']

        G_multi_df.append(G_df[['edge_id', 'g_id', 'from_node', 'to_node', 'length', 'min_time',
                                'max_time', 'min_time_cost', 'max_time_cost', 'min_tariff_cost', 'max_tariff_cost']])

    G_multi_df = pd.concat(G_multi_df, axis=0, sort='False', ignore_index=True)
    # G_multi_df = G_multi_df.reindex(list(G_multi_df.columns)[2:]+list(G_multi_df.columns)[:2], axis=1)
    G_multi_df = G_multi_df[['from_node', 'to_node', 'edge_id', 'g_id', 'length', 'min_time',
                             'max_time', 'min_time_cost', 'max_time_cost', 'min_tariff_cost', 'max_tariff_cost']]
    G_multi = ig.Graph.TupleList(G_multi_df.itertuples(
        index=False), edge_attrs=list(G_multi_df.columns)[2:])
    print(G_multi)
    G_multi_df.to_csv(os.path.join(output_path, 'failure_results', 'multi_modal_network.csv'))

    # modes_file_paths = [('Roads','national_roads'), ('Railways','national_rail'), ('Waterways','waterways'), ('Waterways','waterways'), ('Multi','multi_edges')]
    # modes_file_paths = [('Railways','national_rail')]
    # # modes = ['road','rail','inland','coastal','multi']
    # modes = ['rail']
    # veh_wt = [800]
    # usage_factors = [(0, 0)]
    # speeds = [(40, 60)]
    # fail_mode = ['road','rail']
    # fail_mode = ['rail']
    fail_mode = ['road']
    for m in range(len(fail_mode)):
        # mode_data_path = os.path.join(data_path, modes_file_paths[m][0], modes_file_paths[m][1])
        # for file in os.listdir(mode_data_path):
        #     try:
        #         if file.endswith(".shp") and 'edges' in file.lower().strip():
        #             edges_in = os.path.join(mode_data_path, file)
        #     except:
        #         return ('Network nodes and edge files necessary')

        # if fail_mode[m] == 'road':
        #     G_df =  national_road_shapefile_to_dataframe(edges_in, rd_prop_file)
        # elif fail_mode[m] == 'multi':
        #     G_df = multi_modal_shapefile_to_dataframe(edges_in, md_prop_file, multi_md_len)
        # else:
        #     G_df = network_shapefile_to_dataframe(edges_in, md_prop_file, fail_mode[m], speeds[m][0], speeds[m][1])

        fail_df = pd.read_excel(fail_scenarios_data, sheet_name=modes[m])
        single_ef_list = list(set(fail_df['edge_id'].values.tolist()))
        print('scenarios in national {0} are {1}'.format(modes[m], len(single_ef_list)))
        # """
        # Select individual edge first
        # columns of failure excel are
        # band_name
        # band_num
        # climate_scenario
        # commune_id
        # commune_name
        # district_id
        # district_name
        # edge_id
        # hazard_type
        # max_val
        # min_val
        # probability
        # province_id
        # province_name
        # sector
        # year
        # length
        # """

        # """
        # First do single edge failures
        # """
        flow_df = pd.read_excel(flow_paths_data, sheet_name=modes[m])
        for perct in range(10, 20, 10):
            edge_fail_ranges = []
            for t in range(len(types)):
                flow_df[tons_types[t]] = 0.01*perct*flow_df[tons_types[t]]
                ef_list = []
                for edge in single_ef_list:
                    # ef_dict = igraph_scenario_edge_failures_multi(G_multi_df, [edge], flow_df, veh_wt[m], usage_factors[m][0], usage_factors[m][1], path_types[t], cost_types[t], time_types[t])
                    ef_dict = igraph_scenario_edge_failures_changing_tonnages(G_multi_df, [
                                                                              edge], flow_df, veh_wt[m], usage_factors[m], path_types[t], tons_types[t], cost_types[t], time_types[t])
                    if ef_dict:
                        ef_list += ef_dict

                    print('Done with mode {0} edge {1} type {2}'.format(
                        modes[m], edge, types[t]))

                df = pd.DataFrame(ef_list)
                # df.to_csv(os.path.join(output_path,'failure_results','single_edge_failures_all_paths_national_{0}_{1}_multi_modal_options.csv'.format(modes[m], types[t])), index = False)

                select_cols = ['origin', 'destination', 'o_region', 'd_region', dist_types[t], time_types[t],
                               cost_types[t], vehicle_types[t]] + ind_crop_cols + [rice_type[t], tons_types[t]]
                flow_df_select = flow_df[select_cols]
                flow_df_select = pd.merge(flow_df_select, df, on=[
                                          'origin', 'destination'], how='left').fillna(0)
                flow_df_select = flow_df_select[(flow_df_select[tons_types[t]] > 0) & (
                    flow_df_select['edge_id'] != 0)]

                flow_df_select['dist_diff'] = flow_df_select['new_distance'] - \
                    flow_df_select[dist_types[t]]
                flow_df_select['time_diff'] = flow_df_select['new_time'] - \
                    flow_df_select[time_types[t]]
                # flow_df_select['tr_loss'] = (1 - flow_df_select['no_access'])*flow_df_select[vehicle_types[t]]*(flow_df_select['new_cost'] - flow_df_select[cost_types[t]])
                flow_df_select['tr_loss'] = (
                    1 - flow_df_select['no_access'])*(flow_df_select['new_cost'] - flow_df_select[cost_types[t]])
                df_path = os.path.join(output_path, 'failure_results', 'single_edge_failures_all_path_impacts_national_{0}_{1}_multi_modal_options_{2}_shift.csv'.format(
                    modes[m], types[t], perct))
                flow_df_select.to_csv(df_path, index=False)
                # flow_df_select = pd.read_csv(df_path).fillna(0)
                flow_df_select.rename(columns={'transport_loss': 'tr_loss'}, inplace=True)

                # select_cols = ['edge_id','o_region','d_region','no_access'] + ind_crop_cols + [rice_type[t], tons_types[t]]
                # edge_impact = flow_df_select[select_cols]
                # edge_impact = edge_impact[edge_impact['no_access'] == 1]
                # edge_impact = edge_impact.groupby(['edge_id', 'o_region','d_region'])[ind_crop_cols + [rice_type[t], tons_types[t]]].sum().reset_index()
                # df_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_national_{0}_{1}_multi_modal_options.csv'.format(modes[m], types[t]))
                # edge_impact.to_csv(df_path, index = False)
                # edge_fail_ranges.append(edge_impact)
                # edge_impact = flow_df_select[select_cols+['tr_loss']]
                # edge_impact = edge_impact[edge_impact['no_access'] == 0]
                # edge_impact = edge_impact.groupby(['edge_id', 'o_region','d_region'])['tr_loss'].sum().reset_index()
                # df_path = os.path.join(output_path,'failure_results','single_edge_failures_tr_loss_national_{0}_{1}_multi_modal_options_{2}_shift.csv'.format(modes[m], types[t], perct))
                # edge_impact.to_csv(df_path, index = False)

                # select_cols = ['edge_id','no_access','tr_loss'] + ind_crop_cols + [rice_type[t], tons_types[t]]
                select_cols = ['edge_id', 'no_access', 'tr_loss'] + [tons_types[t]]
                edge_impact = flow_df_select[select_cols]
                edge_impact = edge_impact.groupby(['edge_id', 'no_access'])[
                    select_cols[2:]].sum().reset_index()
                # for col_name in [c for c in select_cols if c not in ['edge_id','no_access', rice_type[t], tons_types[t]]]:
                #     edge_impact.rename(columns={col_name: '{}_'.format(types[t])+col_name}, inplace=True)

                edge_fail_ranges.append(edge_impact)
                del edge_impact

            edge_impact = edge_fail_ranges[0]
            edge_impact = pd.merge(edge_impact, edge_fail_ranges[1], how='left', on=[
                                   'edge_id', 'no_access']).fillna(0)
            df_path = os.path.join(output_path, 'failure_results',
                                   'single_edge_failures_all_losses_national_{0}_multi_modal_options_{1}_shift.csv'.format(modes[m], perct))
            edge_impact.to_csv(df_path, index=False)

            # network_failure_assembly(edge_impact, veh_wt[m], modes[m], G_df, save_edges = True, output_path =shp_output_path)


if __name__ == "__main__":
    main()
