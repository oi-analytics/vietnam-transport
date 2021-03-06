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
    output_path = config['paths']['output']
    data_path = config['paths']['data']

    # Supply input data and parameters
    provinces = ['Lao Cai', 'Binh Dinh', 'Thanh Hoa']
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

    # Process province scale results
    for province in provinces:
        # Load mode network DataFrame
        province_name = province.replace(' ', '').lower()
        print('* Loading {} network DataFrame'.format(province))
        G_df = pd.read_excel(os.path.join(network_data_excel, 'province_roads_edges.xlsx'),
                             sheet_name=province_name, encoding='utf-8')
        G_df = G_df[['edge_id', 'level', 'terrain', 'surface',
                     'road_class', 'road_cond', 'asset_type', 'width', 'length']]
        # Load failure scenarios
        print('* Loading {} failure scenarios'.format(province))
        hazard_scenarios = pd.read_excel(
            os.path.join(fail_scenarios_data, 'province_roads_hazard_intersections.xlsx'),
            sheet_name=province_name)
        hazard_scenarios = hazard_scenarios.drop_duplicates(
            subset=cols, keep='first')

        all_edge_fail_scenarios = combine_hazards_and_network_attributes_and_impacts(
            hazard_scenarios, G_df)

        print('* Creating {} hazard-network scenarios'.format(province))
        scenarios_df = create_hazard_scenarios_for_adaptation(
            all_edge_fail_scenarios, index_cols, length_thr)

        df_path = os.path.join(output_dir,
                               'roads_hazard_intersections_{}_risks.csv'.format(province_name))
        scenarios_df.to_csv(df_path, index=False)
        del scenarios_df


if __name__ == "__main__":
    main()
