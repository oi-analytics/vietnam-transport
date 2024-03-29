# -*- coding: utf-8 -*-
"""Run MRIO
"""
import os
import subprocess

import pandas as pd
from vtra.mrio.functions import (create_proxies, estimate_gva,
                                 load_output, load_provincial_stats)
from vtra.mrio.ras_method import ras_method
from vtra.utils import load_config


def run_mrio_disaggregate(notrade=False, min_rice=True, own_production_ratio=0.8):
    """This function will disaggregate the (single-region) national Input-Output table to a provincial multiregional Input-Output table

    Parameters
        - notrade - Boolean to specify whether we should include trade in the disaggregation. This should be set to **True** in the first step of the disaggregation. The default is set to **False**
        - min_rice - Boolean to determine whether you want to use the minimal rice value or the maximum rice value from the flow analysis. The default is set to **True**
        - own_production_ratio - Specify how much supply and demand is locally supplied and used, and how much is imported/exported. The default is set to **0.8**

    Outputs
        - .csv file containing the new multiregional Input-Output table.
        - pandas DataFrame with a multiregional Input-Output table

    """
    data_path = load_config()['paths']['data']

    # load provincial shapefile
    provinces = load_provincial_stats(data_path)
    provinces.name_eng = provinces.name_eng.apply(
        lambda x: x.replace(' ', '_').replace('-', '_'))

    # estimate gross value added
    provinces['raw_gva'] = estimate_gva(provinces, in_million=True)

    # prepare proxies for settings_trade
    create_proxies(data_path, notrade=notrade,
                   own_production_ratio=own_production_ratio, min_rice=min_rice)

    # run mrio_disaggregate
    if notrade == False:
        p = subprocess.Popen(['mrio_disaggregate', 'settings_trade.yml'],
                             cwd=os.path.join(data_path, 'IO_analysis', 'MRIO_TABLE'))
        p.wait()
    else:
        p = subprocess.Popen(['mrio_disaggregate', 'settings_notrade.yml'],
                             cwd=os.path.join(data_path, 'IO_analysis', 'MRIO_TABLE'))
        p.wait()

    # get reordered mrio with new region classification
    Xin = load_output(data_path, provinces, notrade=notrade)

    # convert to numpy matrix
    X0 = Xin.as_matrix()

    # get sum of rows and columns
    u = X0.sum(axis=1)
    v = X0.sum(axis=0)

    # and only keep T
    v[:(len(u)-3)] = u[:-3]

    # apply RAS method to rebalance the table
    X1 = ras_method(X0, u, v, eps=5e-5)

    # copy new balanced table into dataframe
    Xin.iloc[:, :] = X1

    # add indices to it
    index = list(Xin.index)
    index[567], index[568], index[569] = (
        'total', 'tax_sub'), ('total', 'import_'), ('total', 'valueA')

    Xin.index = pd.MultiIndex.from_tuples(index, names=['region', 'sector'])

    # save outpout
    if notrade == True:
        Xin.to_csv(os.path.join(data_path, 'IO_analysis', 'MRIO_TABLE', 'notrade_trade.csv'))
    else:
        Xin.to_csv(os.path.join(data_path, 'IO_analysis', 'MRIO_TABLE', 'IO_VIETNAM.csv'))

    if notrade == False:
        Region_sum = Xin.groupby(Xin.columns.get_level_values(0), axis='columns').sum().groupby(
            Xin.index.get_level_values(0), axis='index').sum()
        Region_sum.to_csv(os.path.join(data_path, 'IO_analysis',
                                       'MRIO_TABLE', 'region_trade.csv'))

    return Xin


def mrio_to_excel(Xin, min_rice=True):
    """Save the newly created multiregional Input-Output table to Excel, in the format required for the MRIA calculation.

    Parameters
        - Xin - pandas DataFrame of the new multiregional Input-Output table
        - min_rice - Boolean to determine whether you want to use the minimal rice value or the maximum rice value from the flow analysis. The default is set to **True**

    Outputs
        - .xlsx file with the multiregional Input-Output table

    """
    data_path = load_config()['paths']['data']

    Xnew = Xin.copy()

    # prepare exoort and finalD data
    Exports = pd.DataFrame(Xnew.iloc[:, Xnew.columns.get_level_values(
        1) == 'col3'].sum(axis=1), columns=['Exports'])
    Exports.columns = pd.MultiIndex.from_tuples(list(zip(['Total'], ['Export'])))
    FinalD_ToT = Xnew.iloc[:, ((Xnew.columns.get_level_values(1) == 'col1') | (
        Xnew.columns.get_level_values(1) == 'col2'))]
    FinalD_ToT = FinalD_ToT.groupby(level=0, axis=1).sum()
    FinalD_ToT.columns = pd.MultiIndex.from_tuples(
        list(zip(FinalD_ToT.columns, len(FinalD_ToT.columns)*['FinDem'])))

    Xnew.drop(['col1', 'col2', 'col3'], axis=1, level=1, inplace=True)

    Xnew = pd.concat([Xnew, FinalD_ToT, Exports], axis=1)

    # write to excel
    if min_rice == True:
        writer = pd.ExcelWriter(os.path.join(data_path, 'input_data', 'IO_VIETNAM_MIN.xlsx'))
    else:
        writer = pd.ExcelWriter(os.path.join(data_path, 'input_data', 'IO_VIETNAM_MAX.xlsx'))

    # write T
    df_T = Xnew.iloc[:567, :567]
    df_T.columns = df_T.columns.droplevel()
    df_labels_T = pd.DataFrame(df_T.reset_index()[['region', 'sector']])
    df_T.reset_index(inplace=True, drop=True)
    df_T.to_excel(writer, 'T', index=False, header=False)
    df_labels_T.to_excel(writer, 'labels_T', index=False, header=False)

    # write FD
    df_FD = Xnew.iloc[:567, 567:630]
    df_labels_FD = pd.DataFrame(list(df_FD.columns))
    df_FD.columns = df_FD.columns.droplevel()
    df_FD.reset_index(inplace=True, drop=True)
    df_FD.to_excel(writer, 'FD', index=False, header=False)
    df_labels_FD.to_excel(writer, 'labels_FD', index=False, header=False)

    # write ExpROW
    df_ExpROW = Exports[:567]
    df_labels_ExpROW = pd.DataFrame(list(df_ExpROW.columns.get_level_values(1)))
    df_ExpROW.reset_index(inplace=True, drop=True)
    df_ExpROW.columns = df_ExpROW.columns.droplevel()
    df_ExpROW.to_excel(writer, 'ExpROW', index=False, header=False)
    df_labels_ExpROW.reset_index(inplace=True, drop=True)
    df_labels_ExpROW.columns = ['Export']
    df_labels_ExpROW.to_excel(writer, 'labels_ExpROW', index=False, header=False)

    # write VA
    df_VA = pd.DataFrame(Xnew.iloc[567:, :].T[('total', 'tax_sub')] +
                         Xnew.iloc[567:, :].T[('total', 'valueA')], columns=['VA'])
    df_VA['imports'] = Xnew.iloc[567:, :].T[('total', 'import_')]
    df_VA.reset_index(inplace=True, drop=True)
    df_VA.to_excel(writer, 'VA', index=False, header=False)
    df_labels_VA = pd.DataFrame(['Import', 'VA']).T
    df_labels_VA.to_excel(writer, 'labels_VA', index=False, header=False)

    # save excel
    writer.save()


def main():
    run_mrio_disaggregate(notrade=True, min_rice=False)
    Xin = run_mrio_disaggregate(notrade=False, own_production_ratio=0.8, min_rice=False)

    mrio_to_excel(Xin, min_rice=False)

    return Xin


if __name__ == "__main__":
    Xin = main()
