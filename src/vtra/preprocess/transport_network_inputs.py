"""
Purpose
-------
Helper functions to create post-processeed networks with attributes from specific types of input datasets

References
----------
1. Pant, R., Koks, E.E., Russell, T., Schoenmakers, R. & Hall, J.W. (2018).
   Analysis and development of model for addressing climate change/disaster risks in multi-modal transport networks in Vietnam.
   Final Report, Oxford Infrastructure Analytics Ltd., Oxford, UK.
2. All input data folders and files referred to in the code below.
"""
import csv
import os

import geopandas as gpd
import igraph as ig
import networkx as nx
import numpy as np
import pandas as pd
from vtra.utils import line_length


def assign_province_road_conditions(x):
    """
    Assign road conditions as paved or unpaved to Province roads

    Parameters
        x - Pandas DataFrame of values
            - code - Numeric code for type of asset
            - level - Numeric code for level of asset

    Returns
        String value as paved or unpaved
    """
    asset_code = x.code
    asset_level = x.level

    # This is an expressway, national and provincial road
    if asset_code in (17, 303) or asset_level in (0, 1):
        return 'paved'
    else:
        # Anything else not included above
        return 'unpaved'


def assign_assumed_width_to_province_roads_from_file(asset_width, width_range_list):
    """
    Assign widths to Province roads assets in Vietnam

    The widths are assigned based on our understanding of:

    1. The reported width in the data which is not reliable
    2. A design specification based understanding of the assumed width based on ranges of
       values

    Parameters
        - asset_width - Numeric value for width of asset
        - width_range_list - List of tuples containing (from_width, to_width, assumed_width)

    Returns
        assumed_width - assigned width of the raod asset based on design specifications
    """

    assumed_width = asset_width
    for width_vals in width_range_list:
        if width_vals[0] <= assumed_width <= width_vals[1]:
            assumed_width = width_vals[2]
            break

    return assumed_width


def assign_assumed_width_to_province_roads(x):
    """
    Assign widths to Province roads assets in Vietnam

    Parameters
        x : int value for width of asset

    Returns
        int assigned width of the road asset based on design specifications
    """
    if float(x.width) == 0:
        return 4.5
    else:
        return float(x.width)


def assign_asset_type_to_province_roads_from_file(asset_code, asset_type_list):
    """
    Assign asset types to roads assets in Vietnam based on values in file

    The types are assigned based on our understanding of:
    1. The reported asset code in the data

    Parameters
        - asset_code - Numeric value for code of asset
        - asset_type_list - List of Strings wiht names of asset types

    Returns
        asset_type - String name of type of asset
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
        x - Pandas DataFrame with numeric asset code

    Returns
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
        x - Pandas dataframe with values
            - code - Numeric code for type of asset
            - level - Numeric code for level of asset
            - terrain - String value of the terrain of asset

    Returns
        - Float minimum assigned speed in km/hr
        - Float maximum assigned speed in km/hr
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
    Assign time costs on Province roads in Vietnam

    The costs are assigned based on our understanding of:

    1. The types of assets
    2. The levels of classification of assets: 0-National, 1-Provinical, 2-Local, 3-Other
    3. The terrain where the assets are located: Flat or Mountain or No information

    Parameters
        - x - Pandas dataframe with values
            - code - Numeric code for type of asset
            - level - Numeric code for level of asset
            - terrain - String value of the terrain of asset
            - length - Float length of edge in km
            - min_speed - Float minimum assigned speed in km/hr
            - max_speed - Float maximum assigned speed in km/hr
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_time_cost - Float minimum assigned cost of time in USD
        - max_time_cost - Float maximum assigned cost of time in USD
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
    Assign tariff costs on Province roads in Vietnam

    The costs are assigned based on our understanding of:

    1. The types of assets
    2. The levels of classification of assets: 0-National, 1-Provinical, 2-Local, 3-Other
    3. The terrain where the assets are located: Flat or Mountain or No information

    Parameters
        - x - Pandas dataframe with values
            - code - Numeric code for type of asset
            - level - Numeric code for level of asset
            - terrain - String value of the terrain of asset
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_tariff_cost - Float minimum assigned tariff cost in USD/ton
        - max_tariff_cost - Float maximum assigned tariff cost in USD/ton
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


def province_shapefile_to_dataframe(edges_in, road_terrain, road_properties_file,usage_factors):
    """
    Create province network dataframe from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - road_terrain - String name of terrain: flat or mountanious
        - road_properties_file - String path to Excel file with road attributes
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        edges - Geopandas DataFrame with network edge topology and attributes
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
    # width_range_list = [
    #     tuple(x) for x in
    #     pd.read_excel(road_properties_file, sheet_name='widths').values
    # ]
    # edges['width'] = edges.width.apply(
    #     lambda x: assign_assumed_width_to_province_roads_from_file(x, width_range_list))

    edges['width'] = edges.apply(assign_assumed_width_to_province_roads,axis=1)

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

    edges['min_time_cost'] = (1 + usage_factors[0])*edges['min_time_cost']
    edges['max_time_cost'] = (1 + usage_factors[1])*edges['max_time_cost']
    edges['min_tariff_cost'] = (1 + usage_factors[0])*edges['min_tariff_cost']
    edges['max_tariff_cost'] = (1 + usage_factors[1])*edges['max_tariff_cost']

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges[['edge_id','g_id','from_node','to_node'] + add_columns + ['geometry']]
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)

    return edges


def province_shapefile_to_network(edges_in, road_terrain, road_properties_file,usage_factors):
    """
    Create province igraph network from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - road_terrain - String name of terrain: flat or mountanious
        - road_properties_file - String path to Excel file with road attributes
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        G - Igraph object with network edge topology and attributes
    """
    edges = province_shapefile_to_dataframe(edges_in, road_terrain, road_properties_file,usage_factors)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    return G


def assign_national_road_terrain(x):
    """
    Assign terrain as flat or mountain to national roads

    Parameters
        x - Pandas DataFrame of values
            - dia_hinh__ - String value of type of terrain

    Returns
        String value of terrain as flat or mountain
    """
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
    """
    Assign road conditions as paved or unpaved to national roads

    Parameters
        x - Pandas DataFrame of values
            - loai_mat__ - String value of road surface

    Returns
        String value of road as paved or unpaved
    """
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
    """
    Assign road speeds to national roads

    Parameters
        x - Pandas DataFrame of values
            - capkth__ca - String value of road class
            - vehicle_co - Float value of number of vehicles on road

    Returns
        - Integer value of road class
    """
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
    Assign widths to national roads assets in Vietnam

    The widths are assigned based on our understanding of:
    1. The class of the road which is not reliable
    2. The number of lanes
    3. The terrain of the road

    Parameters
        - x - Pandas DataFrame row with values
            - road_class - Integer value of road class
            - lanenum__s - Integer value of number of lanes on road
        - flat_width_range_list - List of tuples containing (from_width, to_width, assumed_width)
        - moiuntain_width_range_list - List of tuples containing (from_width, to_width, assumed_width)

    Returns
        assumed_width - Float assigned width of the road asset based on design specifications
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
        x - Pandas DataFrame of values
            - road_class - Integer value of road class
            - terrain - String value of road terrain
            - est_speed - Float value of estimated speed from CVTS data
        - flat_width_range_list - List of tuples containing design speeds
        - moiuntain_width_range_list - List of tuples containing design speeds

    Returns
        - Float minimum assigned speed in km/hr
        - Float maximum assigned speed in km/hr
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
    Assign time costs on national roads in Vietnam

    The costs are assigned based on our understanding of:

    1. The vehicle counts on roads
    2. The levels of classification of assets: 0-National, 1-Provinical, 2-Local, 3-Other
    3. The terrain where the assets are located: Flat or Mountain or No information

    Parameters
        - x - Pandas dataframe with values
            - vehicle_co - Count of number of vehicles on road
            - code - Numeric code for type of asset
            - level - Numeric code for level of asset
            - terrain - String value of the terrain of asset
            - length - Float length of edge in km
            - min_speed - Float minimum assigned speed in km/hr
            - max_speed - Float maximum assigned speed in km/hr
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_time_cost - Float minimum assigned cost of time in USD
        - max_time_cost - Float maximum assigned cost of time in USD
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
    Assign tariff costs on national roads in Vietnam

    The costs are assigned based on our understanding of:

    1. The vehicle counts on roads

    Parameters
        - x - Pandas dataframe with values
            - vehicle_co - Count of number of vehicles on road
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_tariff_cost - Float minimum assigned tariff cost in USD/ton
        - max_tariff_cost - Float maximum assigned tariff cost in USD/ton
    """
    min_tariff_cost = 0
    max_tariff_cost = 0
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        if cost_param.vehicle_min <= x.vehicle_co < cost_param.vehicle_max:
            min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
            max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length
            break

    return min_tariff_cost, max_tariff_cost


def national_road_shapefile_to_dataframe(edges_in, road_properties_file,usage_factors):
    """
    Create national network dataframe from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - road_properties_file - String path to Excel file with road attributes
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        edges: Geopandas DataFrame with network edge topology and attributes
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

    edges['min_time_cost'] = (1 + usage_factors[0])*edges['min_time_cost']
    edges['max_time_cost'] = (1 + usage_factors[1])*edges['max_time_cost']
    edges['min_tariff_cost'] = (1 + usage_factors[0])*edges['min_tariff_cost']
    edges['max_tariff_cost'] = (1 + usage_factors[1])*edges['max_tariff_cost']

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges[['edge_id','g_id','from_node','to_node'] + add_columns + ['geometry']]
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)

    return edges


def national_road_shapefile_to_network(edges_in, road_properties_file,usage_factors):
    """
    Create national igraph network from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - road_properties_file - String path to Excel file with road attributes
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        G - Igraph object with network edge topology and attributes
    """
    edges = national_road_shapefile_to_dataframe(edges_in, road_properties_file,usage_factors)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def assign_minmax_time_costs_networks_apply(x, cost_dataframe):
    """
    Assign time costs on networks in Vietnam

    Parameters
        - x - Pandas dataframe with values
            - length - Float length of edge in km
            - min_speed - Float minimum assigned speed in km/hr
            - max_speed - Float maximum assigned speed in km/hr
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_time_cost - Float minimum assigned cost of time in USD
        - max_time_cost - Float maximum assigned cost of time in USD
    """
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        min_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.max_speed)
        max_time_cost = 1.0*cost_param.time_cost_usd*(x.length/x.min_speed)

    return min_time_cost, max_time_cost


def assign_minmax_tariff_costs_networks_apply(x, cost_dataframe):
    """
    Assign tariff costs on networks in Vietnam

    Parameters
        - x - Pandas dataframe with values
            - length - Float length of edge in km
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_tariff_cost - Float minimum assigned tariff cost in USD/ton
        - max_tariff_cost - Float maximum assigned tariff cost in USD/ton
    """
    cost_list = list(cost_dataframe.itertuples(index=False))
    for cost_param in cost_list:
        min_tariff_cost = 1.0*cost_param.tariff_min_usd*x.length
        max_tariff_cost = 1.0*cost_param.tariff_max_usd*x.length

    return min_tariff_cost, max_tariff_cost


def network_shapefile_to_dataframe(edges_in, mode_properties_file, mode_name, speed_min, speed_max,usage_factors):
    """
    Create network dataframe from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - mode_properties_file - String path to Excel file with mode attributes
        - mode_name - String name of mode
        - speed_min - Float value of minimum assgined speed
        - speed_max - Float value of maximum assgined speed
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        edges - Geopandas DataFrame with network edge topology and attributes
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


    edges['min_time_cost'] = (1 + usage_factors[0])*edges['min_time_cost']
    edges['max_time_cost'] = (1 + usage_factors[1])*edges['max_time_cost']
    edges['min_tariff_cost'] = (1 + usage_factors[0])*edges['min_tariff_cost']
    edges['max_tariff_cost'] = (1 + usage_factors[1])*edges['max_tariff_cost']

    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges[['edge_id','g_id','from_node','to_node'] + add_columns + ['geometry']]
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)

    return edges


def network_shapefile_to_network(edges_in, mode_properties_file, mode_name, speed_min, speed_max,utilization_factors):
    """
    Create igraph network from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - mode_properties_file - String path to Excel file with mode attributes
        - mode_name - String name of mode
        - speed_min - Float value of minimum assgined speed
        - speed_max - Float value of maximum assgined speed
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        G - Igraph object with network edge topology and attributes
    """
    edges = network_shapefile_to_dataframe(
        edges_in, mode_properties_file, mode_name, speed_min, speed_max,utilization_factors)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def assign_minmax_tariff_costs_multi_modal_apply(x, cost_dataframe):
    """
    Assign tariff costs on multi-modal network links in Vietnam

    Parameters
        - x - Pandas dataframe with values
            - port_type - String name of port type
            - from_mode - String name of mode
            - to_mode - String name of mode
            - other_mode - String name of mode
        - cost_dataframe - Pandas Dataframe with costs

    Returns
        - min_tariff_cost - Float minimum assigned tariff cost in USD/ton
        - max_tariff_cost - Float maximum assigned tariff cost in USD/ton
    """
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


def multi_modal_shapefile_to_dataframe(edges_in, mode_properties_file, mode_name, length_threshold,usage_factors):
    """
    Create multi-modal network dataframe from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - mode_properties_file - String path to Excel file with mode attributes
        - mode_name - String name of mode
        - length_threshold - Float value of threshold in km of length of multi-modal links
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        edges - Geopandas DataFrame with network edge topology and attributes
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

    edges['min_time_cost'] = (1 + usage_factors[0])*edges['min_time_cost']
    edges['max_time_cost'] = (1 + usage_factors[1])*edges['max_time_cost']
    edges['min_tariff_cost'] = (1 + usage_factors[0])*edges['min_tariff_cost']
    edges['max_tariff_cost'] = (1 + usage_factors[1])*edges['max_tariff_cost']
    # make sure that From and To node are the first two columns of the dataframe
    # to make sure the conversion from dataframe to igraph network goes smooth
    edges = edges.reindex(list(edges.columns)[2:]+list(edges.columns)[:2], axis=1)
    edges = edges[edges['length'] < length_threshold]

    return edges


def multi_modal_shapefile_to_network(edges_in, mode_properties_file, mode_name, length_threshold,utilization_factors):
    """
    Create multi-modal igraph network dataframe from inputs

    Parameters
        - edges_in - String path to edges file/network Shapefile
        - mode_properties_file - String path to Excel file with mode attributes
        - mode_name - String name of mode
        - length_threshold - Float value of threshold in km of length of multi-modal links
        - usage_factor - Tuple of 2-float values between 0 and 1

    Returns
        G - Igraph object with network edge topology and attributes
    """
    edges = multi_modal_shapefile_to_dataframe(
        edges_in, mode_properties_file, mode_name, length_threshold,utilization_factors)
    G = ig.Graph.TupleList(edges.itertuples(index=False), edge_attrs=list(edges.columns)[2:])

    # only keep connected network
    return G


def create_port_names(x,port_names_df):
    """
    Add port names in Vietnamese to port data

    Parameters
        - x - Pandas DataFrame with values
            - port_type - String type of port
            - cangbenid - Integer ID of inland port
            - objectid - Integer ID of sea port
        - port_names_df - Pandas DataFrame with port names

    Returns
        name - Vietnamese name of port
    """
    name = ''
    for iter_,port_names in port_names_df.iterrows():
        if (x.port_type == 'inland') and (port_names.port_type == 'inland') and (x.cangbenid == port_names.cangbenid):
            name = port_names.ten
        elif (x.port_type == 'sea') and (port_names.port_type == 'sea') and (x.objectid == port_names.objectid):
            name = port_names.ten_cang

    return name

def read_waterway_ports(ports_file_with_ids, ports_file_with_names):
    """
    Create port data with attributes

    Parameters
        - ports_file_with_ids - String path of GeoDataFrame with port IDs
        - ports_file_with_names - String path of GeoDataFrame with port names

    Returns
        ports_with_id - GeoPandas DataFrame with port attributes
    """
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
    """
    Create port data with attributes

    Parameters
        - ports_file_with_ids - String path of GeoDataFrame with port IDs
        - sector - String path of sector

    Returns
        ports_with_id - GeoPandas DataFrame with port attributes
    """
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
