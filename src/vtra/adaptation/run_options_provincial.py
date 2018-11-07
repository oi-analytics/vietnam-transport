# -*- coding: utf-8 -*-
"""Assess provincial adaptation options
"""
import sys

from vtra.adaptation.adaptation_options import run_adaptation_calculation

if __name__ == 'main':
    read_from_file = False
    if len(sys.argv) == 2:
        if sys.argv[1] == '--read_from_file':
            read_from_file = True

    if read_from_file:
        print('Reading param values from file')
    else:
        print('Running with fixed param values,  --read_from_file to use file')

    config = load_config()
    data_path = config['paths']['data']
    calc_path = config['paths']['calc']
    output_path = config['paths']['output']

    duration_list = [10]
    discount_rate = 12
    growth_rate = 6.5
    regions = ['laocai', 'binhdinh', 'thanhhoa']

    for dur in duration_list:
        for file_id in regions:
            run_adaptation_calculation(
                file_id, data_path, output_path, duration_max=dur, discount_rate=discount_rate, growth_rate=growth_rate, read_from_file=read_from_file)
