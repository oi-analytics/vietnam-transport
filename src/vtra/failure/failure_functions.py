"""Functions used in the provincial and national-scale network failure analysis 
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


def identify_all_failure_paths(network_df_in,edge_failure_set,flow_dataframe,path_criteria):
    """Identify all paths that contain an edge

    Parameters
    ---------
    network_df_in - Pandas DataFrame of network
    edge_failure_set - List of string edge ID's 
    flow_dataframe - Pandas DataFrame of list of edge paths 
    path_criteria - String name of column of edge paths in flow dataframe

    Outputs
    -------
    network_df - Pandas DataFrame of network
        With removed edges
    edge_path_index - List of integer indexes 
        Of locations of paths in flow dataframe
    """

    edge_path_index = []
    network_df = copy.deepcopy(network_df_in)
    for edge in edge_failure_set:
        network_df = network_df[network_df.edge_id != edge]
        edge_path_index += flow_dataframe.loc[flow_dataframe[path_criteria].str.contains(
            "'{}'".format(edge))].index.tolist()
    
    edge_path_index = list(set(edge_path_index))
    return network_df, edge_path_index

def igraph_scenario_edge_failures_changing_tonnages(network_df_in, edge_failure_set, 
    flow_dataframe, vehicle_wt, path_criteria, tons_criteria, cost_criteria, time_criteria):
    """Estimate network impacts of each failures
    When the tariff costs of each path depends on the changing tonnages

    Parameters
    ---------
    network_df_in - Pandas DataFrame of network
    edge_failure_set - List of string edge ID's
    flow_dataframe - Pandas DataFrame of list of edge paths
    vehicle_wt - Float weight of vehcile weight
    path_criteria - String name of column of edge paths in flow dataframe
    tons_criteria - String name of column of path tons in flow dataframe
    cost_criteria - String name of column of path costs in flow dataframe
    time_criteria - String name of column of path travel time in flow dataframe


    Outputs
    -------
    edge_failure_dictionary - List of dictionaries
        With attributes
        edge_id - String name or list of failed edges
        origin - String node ID of Origin of disrupted OD flow
        destination - String node ID of Destination of disrupted OD flow
        no_access - Boolean 1 (no reroutng) or 0 (rerouting)     
        new_cost - Float value of estimated cost of OD journey after disruption 
        new_distance - Float value of estimated distance of OD journey after disruption    
        new_path - List of string edge ID's of estimated new route of OD journey after disruption   
        new_time - Float value of estimated time of OD journey after disruption      
    """
    edge_fail_dictionary = []
    
    network_df,edge_path_index = identify_all_failure_paths(network_df_in,edge_failure_set,flow_dataframe,path_criteria)
    
    if edge_path_index:
        if len(edge_failure_set) == 1:
            edge_failure_set = edge_failure_set[0]

        network_graph = ig.Graph.TupleList(network_df.itertuples(
            index=False), edge_attrs=list(network_df.columns)[2:])

        for e in edge_path_index:
            origin = flow_dataframe.iloc[e]['origin']
            destination = flow_dataframe.iloc[e]['destination']
            origin_node = [x for x in network_graph.vs if x['name'] == origin]
            destination_node = [x for x in network_graph.vs if x['name'] == destination]

            if not origin_node or not destination_node:
                """
                no alternative path exists
                """
                edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': origin, 'destination': destination,
                                             'new_path':[],'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

            else:
                tons = flow_dataframe.iloc[e][tons_criteria]
                vh_nums = math.ceil(1.0*tons/vehicle_wt)
                network_graph = add_igraph_generalised_costs_network(
                    network_graph, vh_nums, tons)
                new_route = network_graph.get_shortest_paths(
                    origin, destination, weights=cost_criteria, output='epath')[0]
                if not new_route:
                    """
                    no alternative path exists
                    """
                    edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': origin, 'destination': destination,
                                                 'new_path':[],'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

                else:
                    new_dist = 0
                    new_time = 0
                    new_gcost = 0
                    new_path = []
                    for n in new_route:
                        new_dist += network_graph.es[n]['length']
                        new_time += network_graph.es[n][time_criteria]
                        new_gcost += network_graph.es[n][cost_criteria]
                        new_path.append(network_graph.es[n]['edge_id'])

                    edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': origin, 'destination': destination,
                                                 'new_path':new_path,'new_distance': new_dist, 'new_time': new_time, 'new_cost': new_gcost, 'no_access': 0})

    return edge_fail_dictionary


def igraph_scenario_edge_failures(network_df_in, edge_failure_set, 
    flow_dataframe, vehicle_wt, path_criteria, 
    tons_criteria, cost_criteria, time_criteria):
    """Estimate network impacts of each failures
    When the tariff costs of each path are fixed by vehicle weight

    Parameters
    ---------
    network_df_in - Pandas DataFrame of network
    edge_failure_set - List of string edge ID's
    flow_dataframe - Pandas DataFrame of list of edge paths
    vehicle_wt - Float weight of vehcile weight
    path_criteria - String name of column of edge paths in flow dataframe
    tons_criteria - String name of column of path tons in flow dataframe
    cost_criteria - String name of column of path costs in flow dataframe
    time_criteria - String name of column of path travel time in flow dataframe


    Outputs
    -------
    edge_failure_dictionary - List of dictionaries
        With attributes
        edge_id - String name or list of failed edges
        origin - String node ID of Origin of disrupted OD flow
        destination - String node ID of Destination of disrupted OD flow
        no_access - Boolean 1 (no reroutng) or 0 (rerouting)     
        new_cost - Float value of estimated cost of OD journey after disruption 
        new_distance - Float value of estimated distance of OD journey after disruption    
        new_path - List of string edge ID's of estimated new route of OD journey after disruption   
        new_time - Float value of estimated time of OD journey after disruption      
    """
    edge_fail_dictionary = []    
    network_df,edge_path_index = identify_all_failure_paths(network_df_in,edge_failure_set,flow_dataframe,path_criteria)
    
    if edge_path_index:
        if len(edge_failure_set) == 1:
            edge_failure_set = edge_failure_set[0]

        network_graph = ig.Graph.TupleList(network_df.itertuples(
            index=False), edge_attrs=list(network_df.columns)[2:])
        network_graph = add_igraph_generalised_costs_network(
            network_graph, 1, vehicle_wt)
        
        nodes_name = np.asarray([x['name'] for x in network_graph.vs])
        select_flows = flow_dataframe[flow_dataframe.index.isin(edge_path_index)]
        
        no_access = select_flows[(~select_flows['origin'].isin(nodes_name)) | (
            ~select_flows['destination'].isin(nodes_name))]
        if len(no_access.index) > 0:
            for iter_, value in no_access.iterrows():
                edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': value['origin'], 'destination': value['destination'],
                                             'new_path':[],'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

        po_access = select_flows[(select_flows['origin'].isin(nodes_name)) & (
            select_flows['destination'].isin(nodes_name))]
        if len(po_access.index) > 0:
            po_access = po_access.set_index('origin')
            origins = list(set(po_access.index.values.tolist()))
            for origin in origins:
                destinations = po_access.loc[[origin], 'destination'].values.tolist()
                paths = network_graph.get_shortest_paths(
                    origin, destinations, weights=cost_criteria, output="epath")
                for p in range(len(paths)):
                    if len(paths[p]) > 0:
                        new_dist = 0
                        new_time = 0
                        new_gcost = 0
                        new_path = []
                        for n in paths[p]:
                            new_dist += network_graph.es[n]['length']
                            new_time += network_graph.es[n][time_criteria]
                            new_gcost += network_graph.es[n][cost_criteria]
                            new_path.append(network_graph.es[n]['edge_id'])
                        edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': origin, 'destination': destinations[p],
                                                     'new_path':new_path,'new_distance': new_dist, 'new_time': new_time, 
                                                     'new_cost': new_gcost, 'no_access': 0})
                    else:
                        edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': origin, 'destination': destinations[p],
                                                     'new_path':[],'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

    return edge_fail_dictionary

def network_failure_assembly(edge_failure_dataframe, gdf_edges, save_edges=True, shape_output_path=''):
    """
    Write results to Shapefiles

    Parameters
    ---------
    edge_failure_dataframe - Pandas DataFrame of edge failure results
    gdf_edges - GeoDataFrame of network edge set with edge ID's and geometry
    save_Edges - Boolean condition to tell code to save created edge shapefile
    shape_output_path - Path where the output shapefile will be stored 

    Outputs
    -------
    gdf_edges - Shapefile 
        With results of edge failure dataframe
    """
    failure_columns = edge_failure_dataframe.columns.values.tolist()
    failure_columns = [f for f in failure_columns if f != 'edge_id']

    for fc in failure_columns:
        gdf_edges[fc] = 0

    for iter_, row in edge_failure_dataframe.iterrows():
        # print (row[1:])
        gdf_edges.loc[gdf_edges['edge_id'] == row['edge_id'],
                      failure_columns] = row[failure_columns].values


    industry_columns = list(set([f.split('min_')[1] for f in failure_columns if 'min' in f]))

    for ind in industry_columns:
        gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(
            x, 'min_{}'.format(ind), 'max_{}'.format(ind)), axis=1)
        gdf_edges[['min_{}'.format(ind), 'max_{}'.format(ind)]
                  ] = gdf_edges['swap'].apply(pd.Series)
        gdf_edges.drop('swap', axis=1, inplace=True)

    if save_edges == True:
        gdf_edges.to_file(shape_output_path)

    del gdf_edges, edge_failure_dataframe

def edge_failure_sampling(failure_scenarios,edge_column):
    """Criteria for selecting failure samples

    Parameters
    ---------
    failure_scenarios - Pandas DataFrame of failure scenarios
    edge_column - String name of column to select failed edge ID's

    Outputs
    -------
    edge_failure_samples - List of lists of failed edge sets    
    """
    edge_failure_samples = list(set(failure_scenarios[edge_column].values.tolist()))

    return edge_failure_samples

def merge_failure_results(flow_df_select,failure_df,tons_col,dist_col,time_col,cost_col,vehicle_col,changing_tonnages=True):
    """Merge failure results with flow results

    Parameters
    ---------
    flow_df_select - Pandas DataFrame of edge flow values
    failure_df - Pandas DataFrame of edge failure values
    tons_col - String name of column of tonnages in flow dataframe
    dist_col - String name of column of distance in flow dataframe
    time_col - String name of column of time in flow dataframe
    cost_col - String name of column of cost in flow dataframe
    vehicle_col - String name of column of vehicle counts in flow dataframe
    changing_tonnages - Boolean True or False 

    Outputs
    -------
    flow_df_select - Pandas DataFrame 
        Of edge flow and failure values merged
    """
    flow_df_select = pd.merge(flow_df_select, failure_df, on=[
                              'origin', 'destination'], how='left').fillna(0)
    flow_df_select = flow_df_select[(flow_df_select[tons_col] > 0) & (flow_df_select['edge_id'] != 0)]

    flow_df_select['dist_diff'] = (1 - flow_df_select['no_access'])*(flow_df_select['new_distance'] - flow_df_select[dist_col])
    flow_df_select['time_diff'] = (1 - flow_df_select['no_access'])*(flow_df_select['new_time'] - flow_df_select[time_col])
    if changing_tonnages == True:
        flow_df_select['tr_loss'] = (1 - flow_df_select['no_access']) * (flow_df_select['new_cost'] - flow_df_select[cost_col])
    else:
        flow_df_select['tr_loss'] = (1 - flow_df_select['no_access'])*flow_df_select[vehicle_col]*(flow_df_select['new_cost'] - flow_df_select[cost_col])

    return flow_df_select
