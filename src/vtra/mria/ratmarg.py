# -*- coding: utf-8 -*-
"""Create output file for gams and run a quick gams module to estimate marginal values of
rationing demand
"""
import os
from shutil import copyfile

import pandas as pd
from gams import GamsParameter, GamsWorkspace
from vtra.utils import load_config


def load_db_IO(table_in):
    """Load the Input-Output data from the **io_basic** Class object and converts it to a GAMS .gdx file.

    Parameters
        - table_in - **io_basic** class object, containing all IO data

    Outputs
        - .gdx file of the IO data

    """

    data_path = load_config()['paths']['data']

    # CREATE GAMS WORKSPACE
    ws = GamsWorkspace(os.path.join(data_path, 'gams_runs'))

    # CREATE INPUT FILES GAMS GDX
    db = ws.add_database()

    # set regions
    reg = db.add_set("reg", 1, "Regions")
    for r in (table_in.countries):
        reg.add_record(r)

    # set rowcol
    rowcol = db.add_set("rowcol", 1, "All rows and columns")

    industries = list(table_in.sectors)
    final_demand = list(table_in.FD_labels['tfd'].unique())

    Import_lab = ['Import']
    Export_lab = ['Export']
    VA_lab = ['VA']

    rowcol_input = industries + final_demand + VA_lab + Import_lab + Export_lab
    for r in (rowcol_input):
        rowcol.add_record(r)

    # set row
    row = db.add_set("row", 1, "All rows")
    row_input = industries + VA_lab + Import_lab
    for r in (row_input):
        row.add_record(r)

    # set col
    col = db.add_set("col", 1, "All columns")
    col_input = industries + final_demand
    for r in (col_input):
        col.add_record(r)

    # set industries
    industries_ = db.add_set("S", 1, "Industries")
    for r in industries:
        industries_.add_record(r)

    # set FinalD
    fd_ = GamsParameter(db, "FinDem_ini", 4, "FinDem")
    for k, v in table_in.FinalD.items():
        fd_.add_record(k).value = v

    # set interaction matrix of intermediate demand
    z_m = db.add_parameter("Z_matrix_ini", 4, "Interaction matrix")
    for k, v in table_in.Z_matrix.items():
        z_m.add_record(k).value = v

    # set interaction matrix of intermediate demand
    a_m = db.add_parameter("A_matrix_ini", 4, "A matrix")
    for k, v in table_in.A_matrix.items():
        a_m.add_record(k).value = v

    # set Export ROW
    exp = db.add_parameter("ExpROW_ini", 3, "Exports to ROW")
    for k, v in table_in.ExpROW.items():
        exp.add_record(k).value = v

    # set ValueA
    val = db.add_parameter("ValueA_ini", 3, "Value Added")
    for k, v in table_in.ValueA.items():
        val.add_record(k).value = v

    # And save to GDX file
    db.export(os.path.join(data_path, "gams_runs", "{}.gdx".format(table_in.name)))

def ratmarg_IO(table_in):
    """Estimate marginal values of the rationing variable in GAMS. GAMS is required, as the marginal values of a variable are not returned in the free python solvers.

    Parameters
        - table_in - **io_basic** class object, containing all IO data

    Outputs
        - pandas DataFrame with the marginal values of the rationing variable

    """

    data_path = load_config()['paths']['data']

    table_in.prep_data()

    load_db_IO(table_in)

    # RUN SCRIPT WITH DISRUPTION
    setdir = os.path.join(data_path, 'gams_runs')
    ws = GamsWorkspace(setdir)
    ws.get_working_directory()

    gamsfile_in = os.path.join(data_path, "gams_runs", "obtain_marg_value.gms")
    gamsfile = os.path.join(data_path, "gams_runs",
                            "obtain_marg_value_{}.gms".format(table_in.name))
    copyfile(gamsfile_in, gamsfile)
    str_ctry = ','.join(table_in.countries)
    str_fd = ','.join(list(table_in.FD_labels['tfd'].unique()))

    with open(gamsfile, 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    gdx_file = "{}.gdx".format(table_in.name)
    data[26] = '$GDXIN '+gdx_file+'\n'

    str_ind = ','.join(table_in.sectors)
    data[32] = 'S(col) list of industries  /'+str_ind+'/\n'

    data[34] = '/'+str_ctry+'/\n'
    data[36] = '/'+str_ctry+'/\n'
    data[38] = '/'+str_fd+'/\n'

    with open(gamsfile, 'w') as file:
        file.writelines(data)

    gamsfile_run = gamsfile.replace("..\\..\\gams_runs\\", "")
    t1 = ws.add_job_from_file(gamsfile_run)

    t1.run()

    Ratmarg = []
    index_ = []
    for rec in t1.out_db["Ratmarg"]:
        index_.append((rec.keys[0], rec.keys[1]))
        Ratmarg.append(rec.get_value())

    index_ = pd.MultiIndex.from_tuples(index_, names=('CNTRY', 'IND'))
    Ratmarginal = pd.DataFrame(Ratmarg, index=index_).unstack()
    Ratmarginal.columns = Ratmarginal.columns.droplevel()

    Ratmarginal.to_csv(os.path.join(data_path, 'input_data',
                                    'Ratmarg_{}.csv'.format(table_in.name)))

    return Ratmarginal
