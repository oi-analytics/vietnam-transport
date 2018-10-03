# -*- coding: utf-8 -*-
"""Helper functions to create transport networks
"""
import csv
import os

import geopandas as gpd
import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
import psycopg2
from vtra.utils import line_length


def assign_province_road_conditions(x):
    asset_code = x.code
    asset_level = x.level

    # This is an expressway, national and provincial road
    if asset_code in (17, 303) or asset_level in (0, 1):
        return 'paved'
    else:
        # Anything else not included above
        return 'unpaved'


def assign_assumed_width_to_province_roads_from_file(asset_width, width_range_list):
    """Assign widths to roads assets in Vietnam

    The widths are assigned based on our understanding of:

    1. The reported width in the data which is not reliable
    2. A design specification based understanding of the assumed width based on ranges of
       values

    Parameters
    ----------
    asset_width - Numeric value for width of asset
    width_range_list - List of tuples containing (from_width, to_width, assumed_width)

    Outputs
    -------
    assumed_width - assigned width of the raod asset based on design specifications
    """

    assumed_width = asset_width
    for width_vals in width_range_list:
        if width_vals[0] <= assumed_width <= width_vals[1]:
            assumed_width = width_vals[2]
            break

    return assumed_width


def assign_assumed_width_to_province_roads(x):
    """Assign widths to roads assets in Vietnam

    The widths are assigned based on our understanding of:

    1. The reported width in the data which is not reliable
    2. A design specification based understanding of the assumed width based on ranges of
       values

    Parameters
    ----------

    x : int
        Numeric value for width of asset

    Returns
    -------

    modified_width : int
        Assigned width of the road asset based on design specifications

    """
    if 0 <= x.width < 4.25:
        return 3.5
    elif 4.25 <= x.width < 6.0:
        return 5.0
    elif 6.0 <= x.width < 8.0:
        return 7.0
    elif 8.0 <= x.width < 11.5:
        return 9.0
    elif 11.5 <= x.width < 17.5:
        return 14.0
    elif 17.5 <= x.width < 24.5:
        return 21.0
    elif 24.5 <= x.width < 100:
        return 9.0
    else:
        return x.width


def assign_asset_type_to_province_roads_from_file(asset_code, asset_type_list):
    """
    Assign asset types to roads assets in Vietnam
    The types are assigned based on our understanding of:
    1. The reported asset code in the data

    Parameters
    ---------
    asset code - Numeric value for code of asset

    Outputs
    -------
    asset type - Which is either of (Bridge, Dam, Culvert, Tunnel, Spillway, Road)
    """
    asset_type = 'road'
    for asset in asset_type_list:
        if asset_code == asset[0]:
            asset_type = asset[2]
            break

    return asset_type


def assign_asset_type_to_province_roads(x):
    """
    Assign asset types to roads assets in Vietnam
    The types are assigned based on our understanding of:
    1. The reported asset code in the data

    Parameters
    ---------
    asset code - Numeric value for code of asset

    Outputs
    -------
    asset type - Which is either of (Bridge, Dam, Culvert, Tunnel, Spillway, Road)
    """
    if x.code in (12, 25):
        return 'Bridge'
    elif x.code == (23):
        return 'Dam'
    elif x.code == (24):
        return 'Culvert'
    elif x.code == (26):
        return 'Tunnel'
    elif x.code == (27):
        return 'Spillway'
    else:
        return 'Road'


def assign_minmax_travel_speeds_province_roads_apply(x):
    """
    Assign travel speeds to roads assets in Vietnam
    The speeds are assigned based on our understanding of:
    1. The types of assets
    2. The levels of classification of assets: 0-National, 1-Provinical, 2-Local, 3-Other
    3. The terrain where the assets are located: Flat or Mountain or No information

    Parameters
    ---------
    asset_code - Numeric code for type of asset
    asset_level - Numeric code for level of asset
    asset_terrain - String value of the terrain of asset

    Outputs
    -------
    speed_min - Minimum assigned speed in km/hr
    speed_max - Maximum assigned speed in km/hr
    """
    asset_code = x.code
    asset_level = x.level
    asset_terrain = x.terrain

    if (not asset_terrain) or ('flat' in asset_terrain.lower()):
        if asset_code == 17:
            # This is an expressway
            return 100, 120
        elif asset_code in (15, 4):
            # This is a residential road or a mountain pass
            return 40, 60
        elif asset_level == 0:
            # This is any other national network asset
            return 80, 100
        elif asset_level == 1:
            # This is any other provincial network asset
            return 60, 80
        elif asset_level == 2:
            # This is any other local network asset
            return 40, 60
        else:
            # Anything else not included above
            return 20, 40

    else:
        if asset_level < 3:
            return 40, 60
        else:
            return 20, 40


def assign_minmax_time_costs_province_roads_apply(x, cost_dataframe):
    """
    """
    asset_code = x.code
    asset_level = x.level
    asset_terrain = x.terrain

    min_time_cost = 0
    max_time_cost = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if cost_param.code == asset_code:
            min_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.max_speed)
            max_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.min_speed)
            break
        elif cost_param.level == asset_level and cost_param.terrain == asset_terrain:
            min_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.max_speed)
            max_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.min_speed)
            break

    return min_time_cost, max_time_cost


def assign_minmax_tariff_costs_province_roads_apply(x, cost_dataframe):
    """
    Assign travel speeds to roads assets in Vietnam
    The speeds are assigned based on our understanding of:
    1. The types of assets
    2. The levels of classification of assets: 0-National, 1-Provinical, 2-Local, 3-Other
    3. The terrain where the assets are located: Flat or Mountain or No information

    Parameters
    ---------
    asset_code - Numeric code for type of asset
    asset_level - Numeric code for level of asset
    asset_terrain - String value of the terrain of asset

    Outputs
    -------
    speed_min - Minimum assigned speed in km/hr
    speed_max - Maximum assigned speed in km/hr
    tariff_min_usd    tariff_max_usd
    """
    asset_code = x.code
    asset_level = x.level
    asset_terrain = x.terrain

    min_tariff_cost = 0
    max_tariff_cost = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if cost_param.code == asset_code:
            min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
            max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length
            break
        elif cost_param.level == asset_level and cost_param.terrain == asset_terrain:
            min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
            max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length
            break

    return min_tariff_cost, max_tariff_cost


def province_shapefile_to_dataframe(edges_in, road_terrain, road_properties_file):
    """
    input parameters:
        edges_in : string of path to edges file/network file.

    output:
        SG: connected graph of the shapefile
    """
    add_columns = ['number','name','terrain','level','surface','road_class',
        'road_cond','asset_type','width','length','min_speed','max_speed',
        'min_time','max_time','min_time_cost','max_time_cost','min_tariff_cost',
        'max_tariff_cost','vehicle_co']
    edges = gpd.read_file(edges_in,encoding='utf-8')
    edges.columns = map(str.lower, edges.columns)

    edges['number'] = ''
    edges['name'] = ''
    edges['surface'] = ''
    edges['road_class'] = ''
    edges['vehicle_co'] = 0
    # assgin asset terrain
    edges['terrain'] = road_terrain

    # assign road conditon
    edges['road_cond'] = edges.apply(assign_province_road_conditions, axis=1)

    # assign asset type
    asset_type_list = [
        tuple(x) for x in
        pd.read_excel(road_properties_file, sheet_name='provincial').values
    ]
    edges['asset_type'] = edges.code.apply(
        lambda x: assign_asset_type_to_province_roads_from_file(x, asset_type_list))

    # get the right linelength
    edges['length'] = edges.geometry.apply(line_length)

    # correct the widths of the road assets
    # get the width of edges
    width_range_list = [
        tuple(x) for x in
        pd.read_excel(road_properties_file, sheet_name='widths').values
    ]
    edges['width'] = edges.width.apply(
        lambda x: assign_assumed_width_to_province_roads_from_file(x, width_range_list))

    # assign minimum and maximum speed to network
    edges['speed'] = edges.apply(assign_minmax_travel_speeds_province_roads_apply, axis=1)
    edges[['min_speed', 'max_speed']] = edges['speed'].apply(pd.Series)
    edges.drop('speed', axis=1, inplace=True)

    # assign minimum and maximum travel time to network
    edges['min_time'] = edges['length']/edges['max_speed']
    edges['max_time'] = edges['length']/edges['min_speed']

    cost_values_df = pd.read_excel(road_properties_file, sheet_name='costs')

    # assign minimum and maximum cost of time in USD to the network
    # the costs of time  = (unit cost of time in USD/hr)*(travel time in hr)
    edges['time_cost'] = edges.apply(
        lambda x: assign_minmax_time_costs_province_roads_apply(x, cost_values_df), axis=1)
    edges[['min_time_cost', 'max_time_cost']] = edges['time_cost'].apply(pd.Series)
    edges.drop('time_cost', axis=1, inplace=True)

    # assign minimum and maximum cost of tonnage in USD/ton to the network
    # the costs of time  = (unit cost of tariff in USD/ton-km)*(length in km)
    edges['tariff_cost'] = edges.apply(
        lambda x: assign_minmax_tariff_costs_province_roads_apply(x, cost_values_df), axis=1)
    edges[['min_tariff_cost', 'max_tariff_cost']] = edges['tariff_cost'].apply(pd.Series)
    edges.drop('tariff_cost', axis=1, inplace=True)

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges[['edge_id','g_id','from_node','to_node'] + add_columns + ['geometry']]
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)

    return edges


def province_shapefile_to_network(edges_in, road_terrain, road_properties_file):
    # create network from edge file
    edges = province_shapefile_to_dataframe(edges_in, road_terrain, road_properties_file)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def assign_national_road_terrain(x):
    terrain_type = x.dia_hinh__

    if terrain_type is None:
        return 'flat'
    elif 'flat' in terrain_type.lower().strip():
        # Assume flat for all roads with no terrain
        return 'flat'
    else:
        # Anything else not included above
        return 'mountain'


def assign_national_road_conditions(x):
    road_cond = x.loai_mat__

    if road_cond is None:
        return 'paved'
    elif 'asphalt' in road_cond.lower().strip():
        # Assume asphalt for all roads with no condition
        return 'paved'
    else:
        # Anything else not included above
        return 'unpaved'


def assign_national_road_class(x):
    road_class = x.capkth__ca
    vehicle_numbers = x.vehicle_co

    if road_class is None:
        if vehicle_numbers >= 6000:
            return 1
        elif 3000 <= vehicle_numbers < 6000:
            return 2
        elif 1000 <= vehicle_numbers < 3000:
            return 3
        elif 300 <= vehicle_numbers < 1000:
            return 4
        elif 50 <= vehicle_numbers < 300:
            return 5
        else:
            return 6
    else:
        if ',' in road_class:
            road_class = road_class.split(',')
        else:
            road_class = [road_class]

        class_1 = [rc for rc in road_class if rc == 'i']
        class_2 = [rc for rc in road_class if rc == 'ii']
        class_3 = [rc for rc in road_class if rc == 'iii']
        class_4 = [rc for rc in road_class if rc == 'iv']
        class_5 = [rc for rc in road_class if rc == 'v']
        class_6 = [rc for rc in road_class if rc == 'vi']

        if class_1:
            return 1
        elif class_2:
            return 2
        elif class_3:
            return 3
        elif class_4:
            return 4
        elif class_5:
            return 5
        elif class_6:
            return 6
        elif vehicle_numbers >= 6000:
            return 1
        elif 3000 <= vehicle_numbers < 6000:
            return 2
        elif 1000 <= vehicle_numbers < 3000:
            return 3
        elif 300 <= vehicle_numbers < 1000:
            return 4
        elif 50 <= vehicle_numbers < 300:
            return 5
        else:
            return 6


def assign_assumed_width_to_national_roads_from_file(x, flat_width_range_list, mountain_width_range_list):
    """
    Assign widths to roads assets in Vietnam
    The widths are assigned based on our understanding of:
    1. The class of the road which is not reliable
    2. The number of lanes
    3. The terrain of the road

    Parameters
    ---------
    x - dataframe row
    flat_width_range_list - List of tuples containing (from_width, to_width, assumed_width)

    Outputs
    -------
    assumed_width - assigned width of the raod asset based on design specifications
    """

    road_class = x.road_class
    road_lanes = x.lanenum__s
    if road_lanes is None:
        road_lanes = 0
    else:
        road_lanes = int(road_lanes)

    road_terrain = x.terrain

    assumed_width = 3.5
    if road_terrain == 'flat':
        for vals in flat_width_range_list:
            if road_class == vals.road_class:
                if road_lanes > 0 and road_lanes <= 8:
                    assumed_width = road_lanes * vals.lane_width + \
                        vals.median_strip + 2.0 * vals.shoulder_width
                else:
                    assumed_width = vals.road_width
                break

    else:
        for vals in mountain_width_range_list:
            if road_class == vals.road_class:
                if road_lanes > 0 and road_lanes <= 8:
                    assumed_width = road_lanes * vals.lane_width + \
                        vals.median_strip + 2.0 * vals.shoulder_width
                else:
                    assumed_width = vals.road_width
                break

    return assumed_width


def assign_min_max_speeds_to_national_roads_from_file(x, flat_width_range_list,
                                                      mountain_width_range_list):
    """
    Assign speeds to national roads in Vietnam
    The speeds are assigned based on our understanding of:
    1. The class of the road
    2. The estimated speed from the CVTS data
    3. The terrain of the road

    Parameters
    ---------
    x - dataframe row
    flat_width_range_list - List of tuples containing flat road properties
    mountain_width_range_list - List of tuples containing mountain road properties

    Outputs
    -------
    min and max speeds - assigned speeds of the road asset based on estimated speeds and
    design specifications
    """

    road_class = x.road_class
    road_terrain = x.terrain
    est_speed = x.est_speed

    min_speed = est_speed
    max_speed = est_speed
    if road_terrain == 'flat':
        for vals in flat_width_range_list:
            if road_class == vals.road_class:
                if est_speed == 0:
                    min_speed = vals.design_speed
                    max_speed = vals.design_speed

                elif est_speed >= vals.design_speed:
                    min_speed = vals.design_speed

                else:
                    max_speed = vals.design_speed

                break

    else:
        for vals in mountain_width_range_list:
            if road_class == vals.road_class:
                if est_speed == 0:
                    min_speed = vals.design_speed
                    max_speed = vals.design_speed

                elif est_speed >= vals.design_speed:
                    min_speed = vals.design_speed

                else:
                    max_speed = vals.design_speed

                break

    return min_speed, max_speed


def assign_minmax_time_costs_national_roads_apply(x, cost_dataframe):
    """
    """
    if x.vehicle_co > 2000:
        asset_code = 17
    else:
        asset_code = 1

    asset_level = 1

    asset_terrain = x.terrain

    min_time_cost = 0
    max_time_cost = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if (cost_param.code == asset_code) and (cost_param.road_cond == x.road_cond):
            min_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.max_speed)
            max_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.min_speed)
            break
        elif (cost_param.level == asset_level) and (cost_param.terrain == asset_terrain) and \
                (cost_param.road_cond == x.road_cond):
            min_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.max_speed)
            max_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.min_speed)
            break

    return min_time_cost, max_time_cost


def assign_minmax_tariff_costs_national_roads_apply(x, cost_dataframe):
    """
    Assign travel speeds to roads assets in Vietnam
    The speeds are assigned based on our understanding of:
    1. The types of assets
    2. The levels of classification of assets: 0-National, 1-Provinical, 2-Local, 3-Other
    3. The terrain where the assets are located: Flat or Mountain or No information

    Parameters
    ---------
    asset_code - Numeric code for type of asset
    asset_level - Numeric code for level of asset
    asset_terrain - String value of the terrain of asset

    Outputs
    -------
    speed_min - Minimum assigned speed in km/hr
    speed_max - Maximum assigned speed in km/hr
    tariff_min_usd    tariff_max_usd
    """
    # if x.vehicle_co > 2000:
    #     asset_code = 17
    # else:
    #     asset_code = 1

    # asset_level = 0
    # asset_terrain= x.terrain

    # min_tariff_cost = 0
    # max_tariff_cost = 0
    # cost_list = list(cost_dataframe.itertuples(index=False))
    # for cost_param in cost_list:
    #     if (cost_param.code == asset_code) and (cost_param.road_cond == x.road_cond):
    #         min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
    #         max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length
    #         break
    #     elif (cost_param.level == asset_level) and (cost_param.terrain == asset_terrain) \
    #            and (cost_param.road_cond == x.road_cond):
    #         min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
    #         max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length
    #         break

    min_tariff_cost = 0
    max_tariff_cost = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if cost_param.vehicle_min <= x.vehicle_co < cost_param.vehicle_max:
            min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
            max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length
            break

    return min_tariff_cost, max_tariff_cost


def national_road_shapefile_to_dataframe(edges_in, road_properties_file):
    """
    input parameters:
        edges_in : string of path to edges file/network file.

    output:
        SG: connected graph of the shapefile
    """
    add_columns = ['number','name','terrain','level','surface','road_class',
        'road_cond','asset_type','width','length','min_speed','max_speed',
        'min_time','max_time','min_time_cost','max_time_cost','min_tariff_cost',
        'max_tariff_cost','vehicle_co']

    edges = gpd.read_file(edges_in,encoding='latin1')
    edges.columns = map(str.lower, edges.columns)

    edges['asset_type'] = ''
    edges['level'] = ''
    # assgin asset terrain
    edges['terrain'] = edges.apply(assign_national_road_terrain, axis=1)

    # assign road conditon
    edges['road_cond'] = edges.apply(assign_national_road_conditions, axis=1)

    # assign road class
    edges['road_class'] = edges.apply(assign_national_road_class, axis=1)

    # get the right linelength
    edges['length'] = edges.geometry.apply(line_length)

    # correct the widths of the road assets
    # get the width of edges
    flat_width_range_list = list(pd.read_excel(
        road_properties_file, sheet_name='flat_terrain_designs').itertuples(index=False))
    mountain_width_range_list = list(pd.read_excel(
        road_properties_file, sheet_name='mountain_terrain_designs').itertuples(index=False))

    edges['width'] = edges.apply(lambda x: assign_assumed_width_to_national_roads_from_file(
        x, flat_width_range_list, mountain_width_range_list), axis=1)

    # assign minimum and maximum speed to network
    edges['speed'] = edges.apply(lambda x: assign_min_max_speeds_to_national_roads_from_file(
        x, flat_width_range_list, mountain_width_range_list), axis=1)
    edges[['min_speed', 'max_speed']] = edges['speed'].apply(pd.Series)
    edges.drop('speed', axis=1, inplace=True)

    # assign minimum and maximum travel time to network
    edges['min_time'] = edges['length']/edges['max_speed']
    edges['max_time'] = edges['length']/edges['min_speed']

    cost_values_df = pd.read_excel(road_properties_file, sheet_name='costs')

    # assign minimum and maximum cost of time in USD to the network
    # the costs of time  = (unit cost of time in USD/hr)*(travel time in hr)
    edges['time_cost'] = edges.apply(
        lambda x: assign_minmax_time_costs_national_roads_apply(x, cost_values_df), axis=1)
    edges[['min_time_cost', 'max_time_cost']] = edges['time_cost'].apply(pd.Series)
    edges.drop('time_cost', axis=1, inplace=True)

    # assign minimum and maximum cost of tonnage in USD/ton to the network
    # the costs of time  = (unit cost of tariff in USD/ton-km)*(length in km)
    edges['tariff_cost'] = edges.apply(
        lambda x: assign_minmax_tariff_costs_national_roads_apply(x, cost_values_df), axis=1)
    edges[['min_tariff_cost', 'max_tariff_cost']] = edges['tariff_cost'].apply(pd.Series)
    edges.drop('tariff_cost', axis=1, inplace=True)

    edges.rename(columns={'ten_duong_':'number','ten_doan__':'name','loai_mat__':'surface'},inplace = True)

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges[['edge_id','g_id','from_node','to_node'] + add_columns + ['geometry']]
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)

    return edges


def national_shapefile_to_dataframe():
    # called from analysis.national_flow_mapping.national_industry_flows
    raise NotImplementedError()


def national_road_shapefile_to_network(edges_in, road_properties_file):
    # create network from edge file
    edges = national_road_shapefile_to_dataframe(edges_in, road_properties_file)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def add_igraph_generalised_costs_roads(G, vehicle_numbers, tonnage):
    # G.es['max_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['max_speed'])))
    # G.es['min_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['min_speed'])))
    # print (G.es['max_time'])
    G.es['max_gcost'] = list(vehicle_numbers*(np.array(G.es['max_time_cost'])
                                              ) + tonnage*(np.array(G.es['max_tariff_cost'])))
    G.es['min_gcost'] = list(vehicle_numbers*(np.array(G.es['min_time_cost'])
                                              ) + tonnage*(np.array(G.es['min_tariff_cost'])))

    return G


def add_igraph_generalised_costs_province_roads():
    # called from failure.{adapt_results_process, adaptation_options_multi,
    # failure_estimation_provinces_multi, failure_projections-multi,
    # failure_scenario_generation}
    raise NotImplementedError()


def add_igraph_time_costs_province_roads():
    # called from analysis.province_flow_mapping.{commune_poi_analysis, province_crop_flows}
    raise NotImplementedError()


def add_igraph_generalised_costs_network(G, vehicle_numbers, tonnage, operating_factor_min,
                                         operating_factor_max):
    # G.es['max_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['max_speed'])))
    # G.es['min_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['min_speed'])))
    # print (G.es['max_time'])
    G.es['max_gcost'] = list(
        (1 + operating_factor_max) *
        (
            vehicle_numbers * np.array(G.es['max_time_cost'])
            + tonnage * np.array(G.es['max_tariff_cost'])
        )
    )
    G.es['min_gcost'] = list(
        (1 + operating_factor_min) *
        (
            vehicle_numbers * np.array(G.es['min_time_cost'])
            + tonnage * np.array(G.es['min_tariff_cost'])
        )
    )

    return G


def add_generalised_costs_network_dataframe(network_dataframe, vehicle_numbers, tonnage,
                                            operating_factor_min, operating_factor_max):
    # G.es['max_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['max_speed'])))
    # G.es['min_cost'] = list(cost_param*(np.array(G.es['length'])/np.array(G.es['min_speed'])))
    # print (G.es['max_time'])
    network_dataframe['max_gcost'] = (1 + operating_factor_max) * (
        vehicle_numbers * network_dataframe['max_time_cost']
        + tonnage * network_dataframe['max_tariff_cost']
    )
    network_dataframe['min_gcost'] = (1 + operating_factor_max) * (
        vehicle_numbers * network_dataframe['min_time_cost']
        + tonnage * network_dataframe['min_tariff_cost']
    )

    return network_dataframe


def assign_minmax_time_costs_networks_apply(x, cost_dataframe):
    """
    """
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        min_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.max_speed)
        max_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.min_speed)

    return min_time_cost, max_time_cost


def assign_minmax_tariff_costs_networks_apply(x, cost_dataframe):
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
        max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length

    return min_tariff_cost, max_tariff_cost


def network_shapefile_to_dataframe(edges_in, mode_properties_file, mode_name, speed_min, speed_max):
    """
    input parameters:
        edges_in : string of path to edges file/network file.

    output:
        SG: connected graph of the shapefile
    """
    add_columns = ['number','name','terrain','level',
        'width','length','min_speed','max_speed',
        'min_time','max_time','min_time_cost','max_time_cost','min_tariff_cost',
        'max_tariff_cost','vehicle_co']

    edges = gpd.read_file(edges_in,encoding='utf-8')
    edges.columns = map(str.lower, edges.columns)

    edges['number'] = ''
    edges['terrain'] = ''
    edges['level'] = ''
    edges['width'] = 0
    edges['vehicle_co'] = 0
    if mode_name == 'rail':
        edges.rename(columns={'railwaylin':'name'},inplace = True)
    elif mode_name in ['inland','coastal']:
        edges.rename(columns={'link':'name'},inplace = True)
    else:
        edges['name'] = ''
    # assgin asset terrain

    # get the right linelength
    edges['length'] = edges.geometry.apply(line_length)

    # assign some speeds
    edges['min_speed'] = speed_min
    edges['max_speed'] = speed_max

    # assign minimum and maximum travel time to network
    edges['min_time'] = edges['length']/edges['max_speed']
    edges['max_time'] = edges['length']/edges['min_speed']

    cost_values_df = pd.read_excel(mode_properties_file, sheet_name=mode_name)

    # assign minimum and maximum cost of time in USD to the network
    # the costs of time  = (unit cost of time in USD/hr)*(travel time in hr)
    edges['time_cost'] = edges.apply(
        lambda x: assign_minmax_time_costs_networks_apply(x, cost_values_df), axis=1)
    edges[['min_time_cost', 'max_time_cost']] = edges['time_cost'].apply(pd.Series)
    edges.drop('time_cost', axis=1, inplace=True)

    # assign minimum and maximum cost of tonnage in USD/ton to the network
    # the costs of time  = (unit cost of tariff in USD/ton-km)*(length in km)
    edges['tariff_cost'] = edges.apply(
        lambda x: assign_minmax_tariff_costs_networks_apply(x, cost_values_df), axis=1)
    edges[['min_tariff_cost', 'max_tariff_cost']] = edges['tariff_cost'].apply(pd.Series)
    edges.drop('tariff_cost', axis=1, inplace=True)

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges[['edge_id','g_id','from_node','to_node'] + add_columns + ['geometry']]
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)

    return edges


def network_shapefile_to_network(edges_in, mode_properties_file, mode_name, speed_min, speed_max):
    # create network from edge file
    edges = network_shapefile_to_dataframe(
        edges_in, mode_properties_file, mode_name, speed_min, speed_max)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def assign_minmax_tariff_costs_multi_modal_apply(x, cost_dataframe):
    min_tariff_cost = 0
    max_tariff_cost = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if cost_param.one_mode == x.port_type and cost_param.other_mode == x.to_mode:
            min_tariff_cost = cost_param.tariff_min_usd
            max_tariff_cost = cost_param.tariff_max_usd
            break
        elif cost_param.one_mode == x.to_mode and cost_param.other_mode == x.from_mode:
            min_tariff_cost = cost_param.tariff_min_usd
            max_tariff_cost = cost_param.tariff_max_usd
            break
        elif cost_param.one_mode == x.to_mode and cost_param.other_mode == x.port_type:
            min_tariff_cost = cost_param.tariff_min_usd
            max_tariff_cost = cost_param.tariff_max_usd
            break
        elif cost_param.one_mode == x.from_mode and cost_param.other_mode == x.to_mode:
            min_tariff_cost = cost_param.tariff_min_usd
            max_tariff_cost = cost_param.tariff_max_usd
            break

    return min_tariff_cost, max_tariff_cost


def multi_modal_shapefile_to_dataframe(edges_in, mode_properties_file, mode_name, length_threshold):
    """
    input parameters:
        edges_in : string of path to edges file/network file.

    output:
        SG: connected graph of the shapefile
    """

    edges = gpd.read_file(edges_in,encoding='utf-8')
    edges.columns = map(str.lower, edges.columns)

    # assgin asset terrain

    # get the right linelength
    edges['length'] = edges.geometry.apply(line_length)

    cost_values_df = pd.read_excel(mode_properties_file, sheet_name=mode_name)

    # assign minimum and maximum cost of tonnage in USD/ton to the network
    # the costs of time  = (unit cost of tariff in USD/ton)
    edges['tariff_cost'] = edges.apply(
        lambda x: assign_minmax_tariff_costs_multi_modal_apply(x, cost_values_df), axis=1)
    edges[['min_tariff_cost', 'max_tariff_cost']] = edges['tariff_cost'].apply(pd.Series)
    edges.drop('tariff_cost', axis=1, inplace=True)

    edges['min_time'] = 0
    edges['max_time'] = 0
    edges['min_time_cost'] = 0
    edges['max_time_cost'] = 0

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)
    edges = edges[edges['length'] < length_threshold]

    return edges


def multi_modal_shapefile_to_network(edges_in, mode_properties_file, mode_name, length_threshold):
    # create network from edge file
    edges = multi_modal_shapefile_to_dataframe(
        edges_in, mode_properties_file, mode_name, length_threshold)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def create_port_names(x,port_names_df):
    name = ''
    for iter_,port_names in port_names_df.iterrows():
        if (x.port_type == 'inland') and (port_names.port_type == 'inland') and (x.cangbenid == port_names.cangbenid):
            name = port_names.ten
        elif (x.port_type == 'sea') and (port_names.port_type == 'sea') and (x.objectid == port_names.objectid):
            name = port_names.ten_cang

    return name

def read_waterway_ports(ports_file_with_ids, ports_file_with_names):
    # load data
    ports_with_name = gpd.read_file(ports_file_with_names, encoding='utf-8')
    ports_with_id = gpd.read_file(ports_file_with_ids, encoding='utf-8')

    ports_with_id.columns = map(str.lower, ports_with_id.columns)
    ports_with_name.columns = map(str.lower, ports_with_name.columns)

    ports_with_id['name'] = ports_with_id.apply(lambda x: create_port_names(x,ports_with_name),axis = 1)
    ports_with_id['population'] = 0
    ports_with_id['capacity'] = 1e9
    ports_with_id = ports_with_id[['node_id','name','port_type','port_class','tons','population','capacity','geometry']]

    return ports_with_id

def read_setor_nodes(node_file_with_ids,sector):
    # load data
    add_columns = [('name',''),('tons',0),('population',0),('capacity',1e9)]
    ports_with_id = gpd.read_file(node_file_with_ids, encoding='utf-8')
    ports_with_id.columns = map(str.lower, ports_with_id.columns)

    if sector == 'air':
        ports_with_id.rename(columns={'ten': 'name'}, inplace=True)
    
    for ac in add_columns:
        if ac[0] not in ports_with_id.columns.values.tolist():
            ports_with_id[ac[0]] = ac[1]
    
    ports_with_id = ports_with_id[['node_id','name','tons','population','capacity','geometry']]

    return ports_with_id