"""Failure analysis of national-scale networks
For transport modes at national scale:
    ['road', 'rail']

Input data requirements
-----------------------
1. Correct paths to all files and correct input parameters
2. MINIMUM MANDATORY DATA SPECIFICATION
    A. Excel file with mode sheets containing network graph structure and attributes
        Should contain following column names and attributes:
            edge_id - String/Integer/Float Edge ID
            from_node - String/Integer/Float node ID that should be present in node_id column
            to_node - String/Integer/Float node ID that should be present in node_id column
            length - Float length of edge in km
            min_time - Float minimum time of travel in hours on edge 
            max_time - Float maximum time of travel in hours on edge  
            min_time_cost - Float minimum cost of time in USD on edge   
            max_time_cost - Float maximum cost of time in USD on edge 
            min_tariff_cost - Float minimum tariff cost in USD on edge   
            max_tariff_cost - Float maximum tariff cost in USD on edge
             
    B. Edge shapefiles for all national-scale networks
        Should contain following column names and attributes:
            edge_id - String/Integer/Float Edge ID
            geometry - Shapely LineString geometry of edges

    C. Excel file with mode sheets containing node-level OD values
        Should contain following column names and attributes:
            origin - String node ID of Origin
            destination - String node ID of Destination
            min_tons -  Float values of minimum daily OD in tons
            max_tons - Float values of maximum daily OD in tons
        Should also contain names of the industry columns specified in the inputs
            
Results
-------
1. Excel sheets with results of flow mapping based on MIN-MAX generalised costs estimates:
        origin - String node ID of Origin
        destination - String node ID of Destination
        min_edge_path - List of string of edge ID's for paths with minimum generalised cost flows
        max_edge_path - List of string of edge ID's for paths with maximum generalised cost flows
        min_distance - Float values of estimated distance for paths with minimum generalised cost flows
        max_distance - Float values of estimated distance for paths with maximum generalised cost flows
        min_time - Float values of estimated time for paths with minimum generalised cost flows
        max_time - Float values of estimated time for paths with maximum generalised cost flows
        min_gcost - Float values of estimated generalised cost for paths with minimum generalised cost flows
        max_gcost - Float values of estimated generalised cost for paths with maximum generalised cost flows
        min_vehicle_nums - Float values of estimated vehicle numbers for paths with minimum generalised cost flows
        max_vehicle_nums - Float values of estimated vehicle numbers for paths with maximum generalised cost flows
        industry_columns - All daily tonnages of industry columns given in the OD matrix data
2. Shapefiles
        edge_id - String/Integer/Float Edge ID
        geometry - Shapely LineString geomtry of edges
        min_{industry} - Float values of estimated minimum daily industries/commodities/total volumes in tons on edges 
        max_{industry} - Float values of estimated maximum daily industries/commodities/total volumes in tons on edges
            For all industry/commodities in the list: 
            ['cement', 'coal', 'constructi', 'fertilizer', 'fishery', 
            'manufactur', 'acof', 'cash', 'cass', 'maiz', 'pepp', 'rcof',
            'rubb', 'swpo', 'teas', 'meat', 'rice', 'petroluem', 'steel', 
            'sugar', 'wood', 'tons'] 

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


def igraph_scenario_edge_failures_changing_tonnages(network_df_in, edge_failure_set, flow_dataframe, vehicle_wt, utilization_factors, path_criteria, tons_criteria, cost_criteria, time_criteria):
    network_graph_df = copy.deepcopy(network_df_in)
    edge_fail_dictionary = []
    edge_path_index = []
    for edge in edge_failure_set:
        network_graph_df = network_graph_df[network_graph_df.edge_id != edge]
        edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains(
            "'{}'".format(edge))].index.tolist()

    edge_path_index = list(set(edge_path_index))
    # print (edge_path_index)
    if edge_path_index:
        network_graph = ig.Graph.TupleList(network_graph_df.itertuples(
            index=False), edge_attrs=list(network_graph_df.columns)[2:])
        # only keep connected network
        # network_graph = network_graph.clusters().giant()

        for e in edge_path_index:
            origin = flow_dataframe.iloc[e]['origin']
            destination = flow_dataframe.iloc[e]['destination']
            origin_node = [x for x in network_graph.vs if x['name'] == origin]
            destination_node = [x for x in network_graph.vs if x['name'] == destination]

            if not origin_node or not destination_node:
                """
                no alternative path exists
                """
                edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destination,
                                             'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

            else:
                tons = flow_dataframe.iloc[e][tons_criteria]
                vh_nums = math.ceil(1.0*tons/vehicle_wt)
                # print (tons, vh_nums, utilization_factors[0], utilization_factors[1])
                network_graph = add_igraph_generalised_costs_network(
                    network_graph, vh_nums, tons, utilization_factors[0], utilization_factors[1])
                new_route = network_graph.get_shortest_paths(
                    origin_node[0], destination_node[0], weights=cost_criteria, output='epath')[0]
                if not new_route:
                    """
                    no alternative path exists
                    """
                    # path_index = path_index_list[e]
                    edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destination,
                                                 'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

                else:
                    new_dist = 0
                    new_time = 0
                    new_gcost = 0
                    for n in new_route:
                        new_dist += network_graph.es[n]['length']
                        new_time += network_graph.es[n][time_criteria]
                        new_gcost += network_graph.es[n][cost_criteria]

                    edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destination,
                                                 'new_distance': new_dist, 'new_time': new_time, 'new_cost': new_gcost, 'no_access': 0})

    return edge_fail_dictionary


def network_od_path_estimations(graph, source, target, cost_criteria, time_criteria):
    # compute min cost paths and values
    paths = graph.get_shortest_paths(source, target, weights=cost_criteria, output="epath")

    edge_path_list = []
    path_dist_list = []
    path_time_list = []
    path_gcost_list = []

    for path in paths:
        edge_path = []
        path_dist = 0
        path_time = 0
        path_gcost = 0
        if path:
            for n in path:
                edge_path.append(graph.es[n]['edge_id'])
                path_dist += graph.es[n]['length']
                path_time += graph.es[n][time_criteria]
                path_gcost += graph.es[n][cost_criteria]

        edge_path_list.append(edge_path)
        path_dist_list.append(path_dist)
        path_time_list.append(path_time)
        path_gcost_list.append(path_gcost)

    return edge_path_list, path_dist_list, path_time_list, path_gcost_list


def igraph_scenario_edge_failures(network_df_in, edge_failure_set, flow_dataframe, vehicle_wt, min_factor, max_factor, path_criteria, cost_criteria, time_criteria):
    network_graph_df = copy.deepcopy(network_df_in)
    edge_fail_dictionary = []
    edge_path_index = []
    for edge in edge_failure_set:
        network_graph_df = network_graph_df[network_graph_df.edge_id != edge]
        edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains(
            "'{}'".format(edge))].index.tolist()

    edge_path_index = list(set(edge_path_index))
    # print (edge_path_index)
    if edge_path_index:
        network_graph = ig.Graph.TupleList(network_graph_df.itertuples(
            index=False), edge_attrs=list(network_graph_df.columns)[2:])
        # only keep connected network
        # network_graph = network_graph.clusters().giant()
        network_graph = add_igraph_generalised_costs_network(
            network_graph, 1, vehicle_wt, min_factor, max_factor)
        nodes_name = np.asarray([x['name'] for x in network_graph.vs])

        select_flows = flow_dataframe[flow_dataframe.index.isin(edge_path_index)]
        no_access = select_flows[(~select_flows['origin'].isin(nodes_name)) | (
            ~select_flows['destination'].isin(nodes_name))]
        if len(no_access.index) > 0:
            for iter_, value in no_access.iterrows():
                edge_dict = {'edge_id': edge,
                             'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1}

                edge_fail_dictionary.append({'edge_id': edge, 'origin': value['origin'], 'destination': value['destination'],
                                             'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

        po_access = select_flows[(select_flows['origin'].isin(nodes_name)) & (
            select_flows['destination'].isin(nodes_name))]
        if len(po_access.index) > 0:
            po_access = po_access.set_index('origin')
            origins = list(set(po_access.index.values.tolist()))
            # print (origins)
            # print (po_access)
            for origin in origins:
                # destinations = po_access.loc[origin,'destination']
                # print (destinations, type(destinations))
                # if isinstance(destinations, str):
                if len(po_access.loc[origin].shape) == 1:
                    destinations = po_access.loc[origin, 'destination']

                    # compute min cost paths and values
                    paths = network_graph.get_shortest_paths(
                        origin, destinations, weights=cost_criteria, output="epath")[0]
                    if len(paths) > 0:
                        new_dist = 0
                        new_time = 0
                        new_gcost = 0
                        for n in paths:
                            new_dist += network_graph.es[n]['length']
                            new_time += network_graph.es[n][time_criteria]
                            new_gcost += network_graph.es[n][cost_criteria]

                        edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destinations,
                                                     'new_distance': new_dist, 'new_time': new_time, 'new_cost': new_gcost, 'no_access': 0})
                    else:
                        edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destinations,
                                                     'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

                else:
                    destinations = po_access.loc[origin, 'destination'].values.tolist()
                    # compute min cost paths and values
                    paths = network_graph.get_shortest_paths(
                        origin, destinations, weights=cost_criteria, output="epath")
                    for p in range(len(paths)):
                        if len(paths[p]) > 0:
                            new_dist = 0
                            new_time = 0
                            new_gcost = 0
                            for n in paths[p]:
                                new_dist += network_graph.es[n]['length']
                                new_time += network_graph.es[n][time_criteria]
                                new_gcost += network_graph.es[n][cost_criteria]

                            edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destinations[p],
                                                         'new_distance': new_dist, 'new_time': new_time, 'new_cost': new_gcost, 'no_access': 0})
                        else:
                            edge_fail_dictionary.append({'edge_id': edge, 'origin': origin, 'destination': destinations[p],
                                                         'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

    return edge_fail_dictionary


def network_failure_assembly(edge_failure_dataframe, vehicle_wt, transport_mode, gdf_edges, save_edges=True, output_path=''):
    """
    Assign net revenue to roads assets in Vietnam

    Parameters
    ---------
    start_points - GeoDataFrame of start points for shortest path analysis.
    end_points - GeoDataFrame of potential end points for shorest path analysis.
    G - iGraph network of the province.
    save_edges -

    Outputs
    -------
    Shapefile with all edges and the total net reveneu transferred along each edge
    GeoDataFrame of total net revenue transferred along each edge
    """
    # min_ind_cols = []
    # max_ind_cols = []
    # ch_min_ind_cols = []
    # ch_max_ind_cols = []
    # for ind in industry_columns:
    #     min_ind_cols.append('min_{}'.format(ind))
    #     max_ind_cols.append('max_{}'.format(ind))
    #     if ind in ('rice','tons'):
    #         ch_min_ind_cols.append('min_{}'.format(ind))
    #         ch_max_ind_cols.append('max_{}'.format(ind))
    #     else:
    #         ch_min_ind_cols.append(ind)
    #         ch_max_ind_cols.append(ind)

    # print (len(ch_min_ind_cols))
    # print (len(ch_max_ind_cols))
    # print (len(min_ind_cols))
    # print (len(max_ind_cols))

    failure_columns = edge_failure_dataframe.columns.values.tolist()
    failure_columns = [f for f in failure_columns if f != 'edge_id']

    for fc in failure_columns:
        gdf_edges[fc] = 0

    for iter_, row in edge_failure_dataframe.iterrows():
        # print (row[1:])
        gdf_edges.loc[gdf_edges['edge_id'] == row['edge_id'],
                      failure_columns] = row[failure_columns].values

    # gdf_edges[min_ind_cols] = gdf_edges['min_vals'].apply(pd.Series)
    # gdf_edges[max_ind_cols] = gdf_edges['max_vals'].apply(pd.Series)
    # gdf_edges.drop('min_vals', axis=1, inplace=True)
    # gdf_edges.drop('max_vals', axis=1, inplace=True)

    industry_columns = list(set([f.split('min_')[1] for f in failure_columns if 'min' in f]))

    for ind in industry_columns:
        gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(
            x, 'min_{}'.format(ind), 'max_{}'.format(ind)), axis=1)
        gdf_edges[['min_{}'.format(ind), 'max_{}'.format(ind)]
                  ] = gdf_edges['swap'].apply(pd.Series)
        gdf_edges.drop('swap', axis=1, inplace=True)

    if save_edges == True:
        gdf_edges.to_file(os.path.join(
            output_path, 'weighted_edges_failures_national_{0}_2.shp'.format(transport_mode)))

    del gdf_edges, edge_failure_dataframe


def main():
    """
    Specify the paths from where you to read and write:
    1. Input data
    2. Intermediate calcuations data
    3. Output results

    Supply input data and parameters
    1. Names of the three Provinces
        List of string types 
    2. Names of modes
        List of strings
    3. Names of output modes
        List of strings
    4. Names of hazard bands
        List of integers
    5. Names of hazard thresholds
        List of integers
    6. Condition 'Yes' or 'No' is the users wants to process
        Province scale results
        National scale results 

    Give the paths to the input data files:
    1. Commune boundary and stats data shapefile
    2. Hazard datasets description Excel file
    3. String name of sheet in hazard datasets description Excel file
    
    Specify the output files and paths to be created 
    """
    data_path, calc_path, output_path = load_config()['paths']['data'], load_config()[
        'paths']['calc'], load_config()['paths']['output']

    """Supply input data and parameters
    """
    modes = ['road', 'rail']
    veh_wt = [20, 800]
    usage_factors = [(0, 0), (0, 0)]

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
    percentage = 100.0
    
    """Give the paths to the input data files
    """
    network_data_path = os.path.join(data_path,'post_processed_networks')
    network_data_excel = os.path.join(data_path,'post_processed_networks','national_edges.xlsx')
    flow_paths_data = os.path.join(output_path, 'flow_mapping_paths',
                                   'national_scale_flow_paths.xlsx')
    fail_scenarios_data = os.path.join(
        output_path, 'hazard_scenarios', 'national_scale_hazard_intersections.xlsx')

    md_prop_file = os.path.join(data_path, 'mode_properties', 'mode_costs.xlsx')
    rd_prop_file = os.path.join(data_path, 'mode_properties', 'road_properties.xlsx')

    """Specify the output files and paths to be created 
    """
    shp_output_path = os.path.join(output_path, 'failure_shapefiles')
    if os.path.exists(shp_output_path) == False:
        os.mkdir(shp_output_path)

    for m in range(len(modes)):
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
        else:
            G_df = network_shapefile_to_dataframe(
                edges_in, md_prop_file, modes[m], speeds[m][0], speeds[m][1])

        """Load mode igraph network and GeoDataFrame
        """
        print ('* Loading {} igraph network and GeoDataFrame'.format(modes[m]))
        G_df = pd.read_excel(network_data_excel,sheet_name = modes[m],encoding='utf-8')
        gdf_edges = gpd.read_file(os.path.join(network_data_path,'{}_edges.shp'.format(modes[m])),encoding='utf-8')
        gdf_edges = gdf_edges[['edge_id','geometry']]

        """Load mode igraph network and GeoDataFrame
        """
        fail_df = pd.read_excel(fail_scenarios_data, sheet_name=modes[m])
        single_ef_list = list(set(fail_df['edge_id'].values.tolist()))
        print('scenarios in national {0} are {1}'.format(modes[m], len(single_ef_list)))

        """
        First do single edge failures
        """
        flow_df = pd.read_excel(flow_paths_data, sheet_name=modes[m])
        edge_fail_ranges = []
        for t in range(len(types)):
            flow_df[tons_types[t]] = 0.01*percentage*flow_df[tons_types[t]]
            ef_list = []
            for edge in single_ef_list:
                ef_dict = igraph_scenario_edge_failures_changing_tonnages(
                    G_df, [edge], flow_df, veh_wt[m], usage_factors[m], path_types[t], tons_types[t], cost_types[t], time_types[t])
                if ef_dict:
                    ef_list += ef_dict

                print('Done with mode {0} edge {1} type {2}'.format(modes[m], edge, types[t]))

            df = pd.DataFrame(ef_list)

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
            flow_df_select['tr_loss'] = (1 - flow_df_select['no_access']) * \
                (flow_df_select['new_cost'] - flow_df_select[cost_types[t]])

            df_path = os.path.join(output_path, 'failure_results', 'single_edge_failures_all_path_impacts_national_{0}_{1}_{2}_shift.csv'.format(
                modes[m], types[t], 100-percentage))
            flow_df_select.to_csv(df_path, index=False)
            flow_df_select.rename(columns={'transport_loss': 'tr_loss'}, inplace=True)

            select_cols = ['edge_id','o_region','d_region','no_access'] + ind_crop_cols + [rice_type[t], tons_types[t]]
            edge_impact = flow_df_select[select_cols]
            edge_impact = edge_impact[edge_impact['no_access'] == 1]
            edge_impact = edge_impact.groupby(['edge_id', 'o_region','d_region'])[ind_crop_cols + [rice_type[t], tons_types[t]]].sum().reset_index()
            df_path = os.path.join(output_path,'failure_results','single_edge_failures_totals_national_{0}_{1}_2.csv'.format(modes[m], types[t]))
            edge_impact.to_csv(df_path, index = False)
            edge_fail_ranges.append(edge_impact)
            edge_impact = flow_df_select[select_cols+['tr_loss']]
            edge_impact = edge_impact[edge_impact['no_access'] == 0]
            edge_impact = edge_impact.groupby(['edge_id', 'o_region','d_region'])['tr_loss'].sum().reset_index()
            df_path = os.path.join(output_path,'failure_results','single_edge_failures_tr_loss_national_{0}_{1}_2.csv'.format(modes[m], types[t]))
            edge_impact.to_csv(df_path, index = False)

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
                               'single_edge_failures_all_losses_national_{0}_{1}_shift.csv'.format(modes[m], 100 - percentage))
        edge_impact.to_csv(df_path, index=False)

        network_failure_assembly(
            edge_impact, veh_wt[m], modes[m], G_df, save_edges=True, output_path=shp_output_path)


if __name__ == "__main__":
    main()
