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

import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from vtra.utils import *

def spatial_scenario_selection(network_shapefile, polygon_shapefile, hazard_dictionary, data_dictionary, network_type ='nodes',name_province =''):
    """Intersect network edges/nodes and boundary Polygons to collect boundary and hazard attributes

    Parameters
        - network_shapefile - Shapefile of edge LineStrings or node Points
        - polygon_shapefile - Shapefile of boundary Polygons
        - hazard_dictionary - Dictionary of hazard attributes
        - data_dictionary - Dictionary of network-hazard-boundary intersection attributes
        - network_type - String value -'edges' or 'nodes' - Default = 'nodes'
        - name_province - String name of province if needed - Default = ''

    Outputs
        data_dictionary - Dictionary of network-hazard-boundary intersection attributes:
            - edge_id/node_id - String name of intersecting edge ID or node ID
            - length - Float length of intersection of edge LineString and hazard Polygon: Only for edges
            - province_id - String/Integer ID of Province
            - province_name - String name of Province in English
            - district_id - String/Integer ID of District
            - district_name - String name of District in English
            - commune_id - String/Integer ID of Commune
            - commune_name - String name of Commune in English
            - hazard_attributes - Dictionary of all attributes from hazard dictionary
    """
    line_gpd = gpd.read_file(network_shapefile)
    poly_gpd = gpd.read_file(polygon_shapefile)


    if len(line_gpd.index) > 0 and len(poly_gpd.index) > 0:
        line_gpd.columns = map(str.lower, line_gpd.columns)
        poly_gpd.columns = map(str.lower, poly_gpd.columns)
        if name_province != '':
            poly_gpd = poly_gpd[poly_gpd['pro_name_e'] == name_province]

        # create spatial index
        poly_sindex = poly_gpd.sindex

        poly_sindex = poly_gpd.sindex
        for l_index, lines in line_gpd.iterrows():
            intersected_polys = poly_gpd.iloc[list(
                poly_sindex.intersection(lines.geometry.bounds))]
            for p_index, poly in intersected_polys.iterrows():
                if (lines['geometry'].intersects(poly['geometry']) is True) and (poly.geometry.is_valid is True) and (lines.geometry.is_valid is True):
                    if network_type == 'edges':
                        value_dictionary = {'edge_id': lines['edge_id'], 'length': 1000.0*line_length(lines['geometry'].intersection(poly['geometry'])),
                                            'province_id': poly['province_i'], 'province_name': poly['pro_name_e'],
                                            'district_id': poly['district_i'], 'district_name': poly['dis_name_e'],
                                            'commune_id': poly['commune_id'], 'commune_name': poly['name_eng']}
                    elif network_type == 'nodes':
                        value_dictionary = {'node_id': lines['node_id'],
                                            'province_id': poly['province_i'], 'province_name': poly['pro_name_e'],
                                            'district_id': poly['district_i'], 'district_name': poly['dis_name_e'],
                                            'commune_id': poly['commune_id'], 'commune_name': poly['name_eng']}

                    data_dictionary.append({**value_dictionary, **hazard_dictionary})

    del line_gpd, poly_gpd
    return data_dictionary

def swap_min_max(x, min_col, max_col):
    """Swap columns if necessary
    """
    if x[min_col] < 0 and x[max_col] < 0:
        if abs(x[min_col]) > abs(x[max_col]):
            return x[max_col], x[min_col]
        else:
            return x[min_col], x[max_col]
    else:
        if x[min_col] > x[max_col]:
            return x[max_col], x[min_col]
        else:
            return x[min_col], x[max_col]

def add_igraph_generalised_costs(G, vehicle_numbers, tonnage):
    # G.es['max_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['max_speed'])))
    # G.es['min_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['min_speed'])))
    # print (G.es['max_time'])
    G.es['max_gcost'] = list(

            vehicle_numbers * np.array(G.es['max_time_cost'])
            + tonnage * np.array(G.es['max_tariff_cost'])
    )
    G.es['min_gcost'] = list(
            vehicle_numbers * np.array(G.es['min_time_cost'])
            + tonnage * np.array(G.es['min_tariff_cost'])
    )

    return G

def network_od_path_estimations(graph,
    source, target, tonnage, vehicle_weight, cost_criteria, time_criteria):
    """Estimate the paths, distances, times, and costs for given OD pair

    Parameters
    ---------
    graph
        igraph network structure
    source
        String/Float/Integer name of Origin node ID
    source
        String/Float/Integer name of Destination node ID
    tonnage : float
        value of tonnage
    vehicle_weight : float
        unit weight of vehicle
    cost_criteria : str
        name of generalised cost criteria to be used: min_gcost or max_gcost
    time_criteria : str
        name of time criteria to be used: min_time or max_time
    fixed_cost : bool

    Returns
    -------
    edge_path_list : list[list]
        nested lists of Strings/Floats/Integers of edge ID's in routes
    path_dist_list : list[float]
        estimated distances of routes
    path_time_list : list[float]
        estimated times of routes
    path_gcost_list : list[float]
        estimated generalised costs of routes

    """
    if vehicle_weight == 0 and tonnage == 0:
        vehicle_weight = 1
        tonnage = 1
    elif vehicle_weight == 0 and tonnage > 0:
        vehicle_weight = tonnage

    graph = add_igraph_generalised_costs(graph, np.ceil(
        tonnage/vehicle_weight), tonnage)

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

def write_flow_paths_to_network_files(save_paths_df,
    industry_columns,min_max_exist,gdf_edges, save_csv=True, save_shapes=True, shape_output_path='',csv_output_path=''):
    """Write results to Shapefiles

    Outputs ``gdf_edges`` - a shapefile with minimum and maximum tonnage flows of all
    commodities/industries for each edge of network.

    Parameters
    ---------
    save_paths_df
        Pandas DataFrame of OD flow paths and their tonnages
    industry_columns
        List of string names of all OD commodities/industries indentified
    min_max_exist
        List of string names of commodity/industry columns for which min-max tonnage column names already exist
    gdf_edges
        GeoDataFrame of network edge set
    save_csv
        Boolean condition to tell code to save created edge csv file
    save_shapes
        Boolean condition to tell code to save created edge shapefile
    shape_output_path
        Path where the output shapefile will be stored
    csv_output_path
        Path where the output csv file will be stored

    """
    if save_shapes == False:
        gdf_edges.drop('geometry', axis=1, inplace=True)

    min_ind_cols = []
    max_ind_cols = []
    ch_min_ind_cols = []
    ch_max_ind_cols = []
    for ind in industry_columns:
        min_ind_cols.append('min_{}'.format(ind))
        max_ind_cols.append('max_{}'.format(ind))
        if ind in min_max_exist:
            ch_min_ind_cols.append('min_{}'.format(ind))
            ch_max_ind_cols.append('max_{}'.format(ind))
        else:
            ch_min_ind_cols.append(ind)
            ch_max_ind_cols.append(ind)

    for i in range(len(min_ind_cols)):
        gdf_edges[min_ind_cols[i]] = 0
        gdf_edges[max_ind_cols[i]] = 0

    for iter_, path in save_paths_df.iterrows():
        min_path = path['min_edge_path']
        max_path = path['max_edge_path']

        gdf_edges.loc[gdf_edges['edge_id'].isin(min_path), min_ind_cols] += path[ch_min_ind_cols].values
        gdf_edges.loc[gdf_edges['edge_id'].isin(max_path), max_ind_cols] += path[ch_max_ind_cols].values


    for ind in industry_columns:
        gdf_edges['swap'] = gdf_edges.apply(lambda x: swap_min_max(x,'min_{}'.format(ind),'max_{}'.format(ind)), axis = 1)
        gdf_edges[['min_{}'.format(ind),'max_{}'.format(ind)]] = gdf_edges['swap'].apply(pd.Series)
        gdf_edges.drop('swap', axis=1, inplace=True)

    if save_shapes == True:
        gdf_edges.to_file(shape_output_path,encoding='utf-8')

    if save_csv == True:
        gdf_edges.to_csv(csv_output_path,index=False,encoding='utf-8')


    del gdf_edges, save_paths_df

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
    flow_dataframe, vehicle_weight, path_criteria, tons_criteria, cost_criteria, time_criteria):
    """Estimate network impacts of each failures
    When the tariff costs of each path depends on the changing tonnages

    Parameters
    ---------
    network_df_in - Pandas DataFrame of network
    edge_failure_set - List of string edge ID's
    flow_dataframe - Pandas DataFrame of list of edge paths
    vehicle_weight - Float weight of vehcile weight
    path_criteria - String name of column of edge paths in flow dataframe
    tons_criteria - String name of column of path tons in flow dataframe
    cost_criteria - String name of column of path costs in flow dataframe
    time_criteria - String name of column of path travel time in flow dataframe


    Returns
    -------
    edge_failure_dictionary : list[dict]
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
                # no alternative path exists
                edge_fail_dictionary.append({'edge_id': edge_failure_set, 'origin': origin, 'destination': destination,
                                             'new_path':[],'new_distance': 0, 'new_time': 0, 'new_cost': 0, 'no_access': 1})

            else:
                tons = flow_dataframe.iloc[e][tons_criteria]
                vh_nums = math.ceil(1.0*tons/vehicle_weight)
                network_graph = add_igraph_generalised_costs(
                    network_graph, vh_nums, tons)
                new_route = network_graph.get_shortest_paths(
                    origin, destination, weights=cost_criteria, output='epath')[0]
                if not new_route:
                    # no alternative path exists
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
    flow_dataframe, vehicle_weight, path_criteria,
    tons_criteria, cost_criteria, time_criteria):
    """Estimate network impacts of each failures
    When the tariff costs of each path are fixed by vehicle weight

    Parameters
    ---------
    network_df_in - Pandas DataFrame of network
    edge_failure_set - List of string edge ID's
    flow_dataframe - Pandas DataFrame of list of edge paths
    vehicle_weight - Float weight of vehcile weight
    path_criteria - String name of column of edge paths in flow dataframe
    tons_criteria - String name of column of path tons in flow dataframe
    cost_criteria - String name of column of path costs in flow dataframe
    time_criteria - String name of column of path travel time in flow dataframe


    Returns
    -------
    edge_failure_dictionary : list[dict]
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
        network_graph = add_igraph_generalised_costs(
            network_graph, 1, vehicle_weight)

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

def rearrange_minmax_values(edge_failure_dataframe):
    """Write results to Shapefiles

    Parameters
    ---------
    edge_failure_dataframe : pandas.DataFrame
        with min-max columns

    Returns
    -------
    edge_failure_dataframe : pandas.DataFrame
        With columns where min < max
    """
    failure_columns = edge_failure_dataframe.columns.values.tolist()
    failure_columns = [f for f in failure_columns if f != ('edge_id','no_access')]

    industry_columns = list(set([f.split('min_')[1] for f in failure_columns if 'min' in f]))

    for ind in industry_columns:
        edge_failure_dataframe['swap'] = edge_failure_dataframe.apply(lambda x: swap_min_max(
            x, 'min_{}'.format(ind), 'max_{}'.format(ind)), axis=1)
        edge_failure_dataframe[['min_{}'.format(ind), 'max_{}'.format(ind)]
                  ] = edge_failure_dataframe['swap'].apply(pd.Series)
        edge_failure_dataframe.drop('swap', axis=1, inplace=True)

    return edge_failure_dataframe

def network_failure_assembly_shapefiles(edge_failure_dataframe, gdf_edges, save_edges=True, shape_output_path=''):
    """Write results to Shapefiles


    Outputs gdf_edges - a Shapefile with results of edge failure dataframe

    Parameters
    ---------
    edge_failure_dataframe
        Pandas DataFrame of edge failure results
    gdf_edges
        GeoDataFrame of network edge set with edge ID's and geometry
    save_edges : bool
        Boolean condition to tell code to save created edge shapefile
    shape_output_path : str
        Path where the output shapefile will be stored

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

    Returns
    -------
    edge_failure_samples - List of lists of failed edge sets
    """
    edge_failure_samples = list(set(failure_scenarios[edge_column].values.tolist()))

    return edge_failure_samples

def merge_failure_results(flow_df_select,failure_df,tons_col,dist_col,time_col,cost_col,vehicle_col,changing_tonnages=True):
    """Merge failure results with flow results

    Parameters
    ---------
    flow_df_select : pandas.DataFrame
        edge flow values
    failure_df : pandas.DataFrame
        edge failure values
    tons_col : str
        name of column of tonnages in flow dataframe
    dist_col : str
        name of column of distance in flow dataframe
    time_col : str
        name of column of time in flow dataframe
    cost_col : str
        name of column of cost in flow dataframe
    vehicle_col : str
        name of column of vehicle counts in flow dataframe
    changing_tonnages : bool

    Returns
    -------
    flow_df_select : pandas.DataFrame
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
