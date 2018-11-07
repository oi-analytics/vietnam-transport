"""Collect network hazard scenarios
"""
import os
import sys

import pandas as pd
from vtra.utils import load_config
from vtra.failure_scenario_selection.hazard_network_scenarios import (
    combine_hazards_and_network_attributes_and_impacts,
    create_hazard_scenarios_for_adaptation
)


def main():
    """Process results

    1. Specify the paths from where you to read and write:
        - Input data
        - Intermediate calcuations data
        - Output results

    2. Supply input data and parameters
        - Names of the three Provinces - List of string types
        - Names of modes - List of strings
        - Names of output modes - List of strings
        - Names of hazard bands - List of integers
        - Names of hazard thresholds - List of integers
        - Condition 'Yes' or 'No' is the users wants to process results

    3. Give the paths to the input data files:
        - Commune boundary and stats data shapefile
        - Hazard datasets description Excel file
        - String name of sheet in hazard datasets description Excel file

    """
    config = load_config()
    data_path = config['paths']['data']
    calc_path = config['paths']['calc']
    output_path = config['paths']['output']

    # Supply input data and parameters
    modes = ['road', 'rail']
    start_year = 2016
    length_thr = 500.0

    cols = [
        'band_num', 'climate_scenario', 'edge_id', 'hazard_type', 'max_val', 'min_val',
        'model', 'probability', 'year', 'length'
    ]
    index_cols = [
        'edge_id', 'hazard_type', 'model', 'climate_scenario', 'year', 'level', 'terrain',
        'surface', 'road_class', 'road_cond', 'asset_type', 'width', 'road_length'
    ]

    # Give the paths to the input data files
    network_data_excel = os.path.join(data_path, 'post_processed_networks')
    fail_scenarios_data = os.path.join(output_path, 'hazard_scenarios')

    # Specify the output files and paths to be created
    output_dir = os.path.join(output_path, 'hazard_scenarios')

    # Process national scale results
    for m in range(len(modes)):
        if modes[m] != 'road':
            index_cols = ['edge_id', 'hazard_type', 'model',
                          'climate_scenario', 'year', 'road_length']
            sel_cols = ['edge_id', 'length']
        else:
            sel_cols = ['edge_id', 'level', 'terrain', 'surface',
                        'road_class', 'road_cond', 'asset_type', 'width', 'length']

        # Load mode network DataFrame
        print('* Loading {} network DataFrame'.format(modes[m]))
        G_df = pd.read_excel(os.path.join(network_data_excel,
                                          'national_edges.xlsx'), sheet_name=modes[m], encoding='utf-8')
        G_df = G_df[sel_cols]

        # Load failure scenarios
        print('* Loading {} failure scenarios'.format(modes[m]))
        hazard_scenarios = pd.read_excel(os.path.join(
            fail_scenarios_data, 'national_scale_hazard_intersections.xlsx'),
            sheet_name=modes[m])
        hazard_scenarios = hazard_scenarios.drop_duplicates(
            subset=cols, keep='first')

        all_edge_fail_scenarios = combine_hazards_and_network_attributes_and_impacts(
            hazard_scenarios, G_df)

        print('* Creating {} hazard-network scenarios'.format(modes[m]))
        scenarios_df = create_hazard_scenarios_for_adaptation(
            all_edge_fail_scenarios, index_cols, length_thr)
        scenarios_df.rename(
            columns={'road_length': '{}_length'.format(modes[m])}, inplace=True)

        df_path = os.path.join(output_path, 'hazard_scenarios',
                               'national_{}_hazard_intersections_risks.csv'.format(modes[m]))
        scenarios_df.to_csv(df_path, index=False)
        del scenarios_df


if __name__ == "__main__":
    main()
