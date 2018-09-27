"""Sum max/min total flow exposed under hazard scenarios at air and water network nodes
"""
import os
import sys

import geopandas as gpd
import pandas as pd


from vtra.utils import *

def main():
    config = load_config()

    # Output to Excel
    # sheets: air, inland, coastal
    flow_output_file = os.path.join(config['paths']['data'], 'Results', 'flow_mapping_paths', 'air-water-node-flows.xlsx')
    impact_output_file = os.path.join(config['paths']['data'], 'Results', 'Failure_results', 'air-water-node-impacts.xlsx')
    impact_summary_file = os.path.join(config['paths']['data'], 'Results', 'Failure_results', 'air-water-node-impacts-summary.xlsx')

    # Output to Shapefiles
    airports_flow_output_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'airport_nodes_flows.shp')
    inland_ports_flow_output_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'inland_ports_nodes_flows.shp')
    coastal_ports_flow_output_file = os.path.join(config['paths']['data'], 'Results', 'Flow_shapefiles', 'coastal_ports_nodes_flows.shp')

    # locate data
    airports_file = os.path.join(
        config['paths']['data'],
        'Airports', 'airnetwork', 'airportnodes.shp')

    ports_file_with_names = os.path.join(
        config['paths']['data'],
        'Waterways', 'waterways', 'ports.shp')
    ports_file_with_ids = os.path.join(
        config['paths']['data'],
        'Waterways', 'waterways', 'ports_nodes.shp')

    flows_file = os.path.join(
        config['paths']['data'],
        'Results', 'flow_mapping_paths', 'national_scale_flow_paths.xlsx')
    hazard_file = os.path.join(
        config['paths']['data'],
        'Results', 'Hazard_network_intersections', 'national_scale_hazard_intersections.xlsx')

    # read data
    print(" * Reading data")
    airports = read_airports(airports_file)
    inland_ports, coastal_ports = read_ports(ports_file_with_ids, ports_file_with_names)
    air_exposure, inland_exposure, coastal_exposure = read_hazards(hazard_file)
    air_flows, water_flows = read_flows(flows_file)

    # aggregate flows to nodes
    print(" * Aggregating flows")
    airports_with_flows = aggregate_flows(airports, air_flows)
    inland_ports_with_flows = aggregate_flows(inland_ports, water_flows)
    coastal_ports_with_flows = aggregate_flows(coastal_ports, water_flows)

    # as gdf
    airports_with_flows_gdf = gpd.GeoDataFrame(airports_with_flows, crs={'init': 'epsg:4326'}, geometry=airports_with_flows.geometry)
    inland_ports_with_flows_gdf = gpd.GeoDataFrame(inland_ports_with_flows, crs={'init': 'epsg:4326'}, geometry=inland_ports_with_flows.geometry)
    coastal_ports_with_flows_gdf = gpd.GeoDataFrame(coastal_ports_with_flows, crs={'init': 'epsg:4326'}, geometry=coastal_ports_with_flows.geometry)

    # save flows
    airports_with_flows_gdf.to_file(airports_flow_output_file)
    inland_ports_with_flows_gdf.to_file(inland_ports_flow_output_file)
    coastal_ports_with_flows_gdf.to_file(coastal_ports_flow_output_file)
    with pd.ExcelWriter(flow_output_file) as writer:
        airports_with_flows.to_excel(writer, sheet_name='air', index=False)
        inland_ports_with_flows.to_excel(writer, sheet_name='inland', index=False)
        coastal_ports_with_flows.to_excel(writer, sheet_name='coastal', index=False)

    # join hazards
    print(" * Joining hazards")
    airports_with_hazards = join_hazards(airports_with_flows, air_exposure)
    inland_ports_with_hazards = join_hazards(inland_ports_with_flows, inland_exposure)
    coastal_ports_with_hazards = join_hazards(coastal_ports_with_flows, coastal_exposure)

    # save hazards
    with pd.ExcelWriter(impact_output_file) as writer:
        airports_with_hazards.to_excel(writer, sheet_name='air', index=False)
        inland_ports_with_hazards.to_excel(writer, sheet_name='inland', index=False)
        coastal_ports_with_hazards.to_excel(writer, sheet_name='coastal', index=False)

    # summarise
    print(" * Summarising")
    airports_summary = summarise(airports_with_hazards)
    inland_ports_summary = summarise(inland_ports_with_hazards)
    coastal_ports_summary = summarise(coastal_ports_with_hazards)

    # save summaries
    with pd.ExcelWriter(impact_summary_file) as writer:
        airports_summary.to_excel(writer, sheet_name='air')
        inland_ports_summary.to_excel(writer, sheet_name='inland')
        coastal_ports_summary.to_excel(writer, sheet_name='coastal')
    print(" * Done")


def read_airports(airports_file):
    return gpd.read_file(airports_file, encoding='utf-8')[[
        'node_id', 'ten', 'ma_iata', 'geometry'
    ]].rename({
        'ten': 'name',
        'ma_iata': 'iata_location_code'
    }, axis=1)


def read_ports(ports_file_with_ids, ports_file_with_names):
    # load data
    ports_with_name = gpd.read_file(ports_file_with_names, encoding='utf-8')
    ports_with_id = gpd.read_file(ports_file_with_ids, encoding='utf-8')

    # merge inland ports
    left = ports_with_name[
        ports_with_name.port_type == 'inland'
    ][
        ['cangbenid', 'ten']
    ]
    right = ports_with_id[
        ports_with_id.PORT_TYPE == 'inland'
    ][
        ['NODE_ID', 'CANGBENID', 'geometry']
    ]
    inland_ports = pd.merge(
        left, right, left_on='cangbenid', right_on='CANGBENID', validate='one_to_one'
    ).drop(
        ['cangbenid', 'CANGBENID'], axis=1
    ).rename({
        'NODE_ID': 'node_id',
        'ten': 'name'
    }, axis=1)

    # merge sea ports (only those with port_class not none)
    left = ports_with_name[
        ports_with_name.port_type == 'sea'
    ][
        ['objectid', 'ten_cang']
    ]
    right = ports_with_id[
        (ports_with_id.PORT_TYPE == 'sea') & (ports_with_id.port_class != 'none')
    ][
        ['NODE_ID', 'OBJECTID', 'port_class', 'geometry']
    ]
    coastal_ports = pd.merge(
        left, right, left_on='objectid', right_on='OBJECTID', validate='one_to_one', how='right'
    ).drop(
        ['objectid', 'OBJECTID'], axis=1
    ).rename({
        'NODE_ID': 'node_id',
        'ten_cang': 'name'
    }, axis=1)

    return inland_ports, coastal_ports


def read_hazards(hazard_file):
    data = pd.read_excel(hazard_file, ['air', 'inland', 'coastal'])
    keep_cols = ['node_id','commune_name','district_name','province_name', 'hazard_type', 'model', 'climate_scenario', 'probability', 'year']
    air_exposure = data['air'][keep_cols]
    inland_exposure = data['inland'][keep_cols]
    coastal_exposure = data['coastal'][keep_cols]
    return air_exposure, inland_exposure, coastal_exposure


def read_flows(flows_file):
    data = pd.read_excel(flows_file, ['air', 'inland', 'coastal'])
    keep_cols = ['origin', 'destination', 'min_tons', 'max_tons']
    air_flows = data['air'][keep_cols]
    inland_flows = data['inland'][keep_cols]
    coastal_flows = data['coastal'][keep_cols]
    # inland and coastal ports may appear in either set, so merge
    water_flows = pd.concat([inland_flows, coastal_flows])
    return air_flows, water_flows


def aggregate_flows(nodes_df, flows_df):
    flow_ids = pd.concat([flows_df.origin, flows_df.destination]).unique()
    flow_nodes = nodes_df[nodes_df.node_id.isin(flow_ids)]

    out_flows = flows_df.groupby('origin').sum()
    with_out = pd.merge(
        flow_nodes, out_flows, left_on='node_id', right_on='origin', how='left'
    ).rename(columns={
        'min_tons': 'min_tons_out',
        'max_tons': 'max_tons_out'
    }).fillna(0)

    in_flows = flows_df.groupby('destination').sum()
    with_in = pd.merge(
        with_out, in_flows, left_on='node_id', right_on='destination', how='left'
    ).rename(columns={
        'min_tons': 'min_tons_in',
        'max_tons': 'max_tons_in'
    }).fillna(0)

    nodes_with_flows = with_in
    nodes_with_flows['min_tons'] = nodes_with_flows.min_tons_in + nodes_with_flows.min_tons_out
    nodes_with_flows['max_tons'] = nodes_with_flows.max_tons_in + nodes_with_flows.max_tons_out
    return nodes_with_flows


def join_hazards(nodes_with_flows_df, hazards_df):
    return pd.merge(
        nodes_with_flows_df, hazards_df, how='inner', validate='one_to_many', on='node_id'
    )


def summarise(nodes_with_hazards_df):
    grouped = nodes_with_hazards_df[
        ['name','commune_name','district_name','province_name', 'min_tons', 'max_tons', 'hazard_type', 'climate_scenario', 'probability']
    ].groupby(
        ['name','commune_name','district_name','province_name', 'min_tons', 'max_tons', 'hazard_type','climate_scenario']
    )

    min_prob = grouped.min(
    ).rename(columns={'probability': 'min_probability'})

    max_prob = grouped.max(
    )[
        ['probability']
    ].rename(columns={'probability': 'max_probability'})

    summary = pd.concat(
        [min_prob, max_prob],
        axis=1
    ).sort_values(
        by=['max_tons', 'hazard_type', 'climate_scenario'], ascending=[False, True, True]
    ).rename(
        columns={
            'min_probability': 'Probability (minimum)',
            'max_probability': 'Probability (maximum)',
        }
    )

    summary.index.names = [
        'Name',  # was 'name'
        'Commune',
        'District',
        'Province',
        'Minimum flow (tons/day)',  # was 'min_tons'
        'Maximum flow (tons/day)',  # was 'max_tons'
        'Hazard type',  # was 'hazard_type'
        'Climate scenario',  # was 'climate_scenario'
    ]

    return summary


if __name__ == '__main__':
    main()
