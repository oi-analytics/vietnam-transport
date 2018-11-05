# -*- coding: utf-8 -*-
"""Create disruption files to be used in the MRIA loss estimations.
"""
import os

import numpy as np
import pandas as pd
from vtra.utils import load_config

pd.options.mode.chained_assignment = None

def create_disruption(input_file, output_dir, min_rice=True, single_point=True):
    """Create disruption file for the economic analysis.

    The input for this disruption file is the outcome of the flow analysis. This function
    translate the failure in transport flows into impacts to the economic sectors.

    Parameters
    ----------
    input_file : str
        name of the path to a flow failure scenario file
    output_dir : str
        name for output directory for mapper file.
    min_rice : bool, optional
        determine whether you want to use the minimal rice value or the maximum rice value from
        the flow analysis. This MUST match the value used when creating the MRIO table. The
        default is **True**.
    single_point : bool, optional
        determine whether you are converting a single-point or multi-point failure analaysis to
        a disruption file. The default is **True**.

    Returns
    -------
    event_dict
        Dictionary of all unique failure events, in terms of percentage disruption per sector
        in each directly affected region due to the flow failure.
    """

     # Define current directory and data directory
    data_path = load_config()['paths']['data']

    # read national IO
    comm_des = df_com_to_ind(pd.read_excel(os.path.join(
        data_path, 'OD_data', 'national_scale_od_matrix.xlsx'), sheet_name='total'), min_rice=min_rice)

   # read vietnam failure results
    vnm_failure_all_failure = pd.read_csv(input_file, index_col=[0])

    if single_point == False:
        mapper = dict(zip(vnm_failure_all_failure.index.unique(), ['multi_{}'.format(
            n) for n in np.arange(1, len(vnm_failure_all_failure.index.unique())+1, 1)]))
        vnm_failure_all_failure.index = vnm_failure_all_failure.index.map(lambda x: mapper[x])
        pd.DataFrame.from_dict(mapper, orient='index').to_csv(
            os.path.join(output_dir, 'mapper.csv'))

    event_dict = {}
    for id_, scenario in vnm_failure_all_failure.groupby(level=0):
        if id_ == 'multi_102':
            break
        disruption = df_com_to_ind(scenario, min_rice=min_rice).div(
            comm_des).dropna(axis=0, how='all').fillna(0)
        disruption = pd.DataFrame(disruption.stack(0))
        disruption.columns = ['value']
        disruption = disruption.loc[disruption.value > 0]
        disruption['value'] = 1 - (disruption['value']*0.2)
        event_dict[id_] = disruption['value'].to_dict()

    return event_dict

def map_comm_ind(x):
    """Map the goods from the flow failure analysis to economic sectors.

    Parameters
    ----------
    x : str
        row in the disruption dataframe.

    Returns
    -------
    x : str
        mapped good to sector for the specific row in the disruption dataframe

    """

    comm_ind_map = {
        'acof': 'Agriculture',
        'cash': 'Agriculture',
        'cass': 'Agriculture',
        'cement': 'Processing',
        'coal': 'Processing',
        'constructi': 'Construction',
        'fertilizer': 'Processing',
        'fishery': 'Agriculture',
        'maiz': 'Agriculture',
        'manufactur': 'Manufacturing',
        'meat': 'Agriculture',
        'min_rice': 'Agriculture',
        'max_rice': 'Agriculture',
        'pepp': 'Agriculture',
        'petroluem': 'Processing',
        'rcof': 'Agriculture',
        'rubb': 'Processing',
        'steel': 'Processing',
        'sugar': 'Agriculture',
        'swpo': 'Processing',
        'teas': 'Agriculture',
        'wood': 'Wood and Paper'
    }
    return comm_ind_map[x]


def map_ind(x):
    """Map the abbreviated names for the industries to their full name.

    Parameters
    ----------
    x : str
        row in the disruption dataframe.

    Returns
    -------
    x : str
        mapped abbrevation to full name for the specific row in the disruption dataframe.

    """
    ind_map = {
        'secA': 'Agriculture',
        'secB': 'Mining',
        'secC': 'Processing',
        'secD': 'Textile and Garment',
        'secE': 'Wood and Paper',
        'secF': 'Manufacturing',
        'secG': 'Construction',
        'secH': 'Trade',
        'secI': 'Services'
    }
    return {v: k for k, v in ind_map.items()}[x]


def df_com_to_ind(comm_des, min_rice=True):
    """Convert the national Origin-Destination matrix from goods to sectors.

    Parameters
    ----------
    comm_dess
        national Origin-Destination matrix, showing origin and destination of goods.
    min_rice : bool
        determine whether to use the minimal rice value or the maximum rice value from the flow
        analysis. This MUST match the value used when creating the MRIO table. The default is
        **True**.

    Returns
    -------
    df_od_ind
        a Pandas dataframe of the national OD matrix on a sector level.

    """
    comm_des.o_region = comm_des.o_region.apply(
        lambda x: x.replace(' ', '_').replace('-', '_'))
    comm_des.d_region = comm_des.d_region.apply(
        lambda x: x.replace(' ', '_').replace('-', '_'))

    comm_des = comm_des.groupby([comm_des.o_region, comm_des.d_region]).sum()
    df_od = comm_des.stack(0).reset_index()
    df_od.drop(['o_region'], inplace=True, axis=1)
    df_od.columns = ['d_region', 'good', 'value']

    if min_rice == True:
        df_od = df_od.loc[~(df_od.good.isin(['max_rice', 'min_tons', 'max_tons']))]
    else:
        df_od = df_od.loc[~(df_od.good.isin(['min_rice', 'min_tons', 'max_tons']))]

    df_od['good'] = df_od.good.apply(lambda x: map_comm_ind(x))
    df_od['good'] = df_od.good.apply(lambda x: map_ind(x))

    df_od_ind = df_od.groupby(['d_region', 'good']).sum()
    df_od_ind = df_od_ind.unstack(1)
    df_od_ind.columns = df_od_ind.columns.get_level_values(1)

    return df_od_ind


if __name__ == "__main__":
    data_path = load_config()['paths']['data']

    input_file = os.path.join(
        data_path, 'Results', 'Failure_results',
        'single_edge_failures_totals_national_road_min.csv')
