# -*- coding: utf-8 -*-
"""
Create output file for gams and run a quick gams module to estimate marginal values of rationing demand

@author: Elco Koks

@date: Nov, 2017
"""
import os
import sys
import pandas as pd
from gams import GamsWorkspace, GamsParameter
from shutil import copyfile

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
from scripts.utils import load_config


def load_db_SUT(table_in,RoW=None):
    
    data_path = load_config()['paths']['data']
    
    '''CREATE GAMS WORKSPACE'''
    
    ws = GamsWorkspace(os.path.join(data_path,'gams_runs'))    

    ''' CREATE INPUT FILES GAMS GDX '''
    
    db = ws.add_database()

    #set regions
    reg = db.add_set("reg",1,"Regions")
    for r in (table_in.countries):
        reg.add_record(r)        
   
    #set rowcol
    rowcol = db.add_set("rowcol",1,"All rows and columns")   
    industries = list(table_in.sectors)  
    products = list(table_in.products)
    final_demand = ['FinalD']

    Import_lab  = ['Import']
    Export_lab  = ['Export']
    VA_lab = ['VA']
    
    rowcol_input = industries + final_demand + Export_lab  + products + VA_lab 
    for r in (rowcol_input):
        rowcol.add_record(r) 
    
    #set row
    row = db.add_set("row",1,"All rows")    
    row_input = products + VA_lab + Import_lab
    for r in (row_input):
        row.add_record(r) 
      
    #set col
    col = db.add_set("col",1,"All columns")    
    col_input = industries + final_demand 
    for r in (col_input):
        col.add_record(r) 
        
    #set industries
    industries_ = db.add_set("ind",1,"Industries")
    for r in  industries:
        industries_.add_record(r)        
    
    #set Use table
    use_m = db.add_parameter("REG_USE2013", 4, "Interaction matrix")
    for k, v in table_in.Use.items():
        use_m.add_record(k).value = v 
    
    #set Supply table
    sup_m = db.add_parameter("REG_SUP2013", 4, "Interaction matrix")
    for k, v in table_in.Sup.items():
        sup_m.add_record(k).value = v 

    #set export ROW    
    exp = db.add_parameter("ExpROW_ini", 3, "Exports to ROW")
    for k, v in table_in.ExpROW.items():
        exp.add_record(k).value = v 

    #set export ROW    
    imp = db.add_parameter("ImpROW_ini", 3, "Imports from ROW")
    for k, v in table_in.ImpROW.items():
        imp.add_record(k).value = v 
 
    #set ValueA
    val = db.add_parameter("ValueA_ini", 3, "Value Added")
    for k, v in table_in.ValueA.items():
        val.add_record(k).value = v 
        
    # And save to GDX file
    db.export(os.path.join(data_path,"gams_runs","{}.gdx".format(table_in.name)))   
   
    
def load_db_IO(table_in,EORA=False,RoW=None):

    data_path = load_config()['paths']['data']

    
    '''CREATE GAMS WORKSPACE'''
    
    ws = GamsWorkspace(os.path.join(data_path,'gams_runs'))
    
    ''' CREATE INPUT FILES GAMS GDX '''
    
    db = ws.add_database()
    
    #set regions
    reg = db.add_set("reg",1,"Regions")
    if EORA is True:
        for r in (table_in.countries+['ROW']):
            reg.add_record(r)
    else:
        for r in (table_in.countries):
            reg.add_record(r)        
   
    #set rowcol
    rowcol = db.add_set("rowcol",1,"All rows and columns")   
    if EORA is True:
        industries = list(table_in.sectors)  + ['Total']
        final_demand = list(table_in.FD_labels['FD'].unique())

    else:
        industries = list(table_in.sectors)  
        final_demand = list(table_in.FD_labels['tfd'].unique())

    Import_lab  = ['Import']
    Export_lab  = ['Export']
    VA_lab = ['VA']
    
    rowcol_input = industries + final_demand + VA_lab + Import_lab + Export_lab
    for r in (rowcol_input):
        rowcol.add_record(r) 
    
    #set row
    row = db.add_set("row",1,"All rows")    
    row_input = industries + VA_lab + Import_lab
    for r in (row_input):
        row.add_record(r) 
      
    #set col
    col = db.add_set("col",1,"All columns")    
    col_input = industries + final_demand 
    for r in (col_input):
        col.add_record(r) 
        
    #set industries
    industries_ = db.add_set("S",1,"Industries")
    for r in  industries:
        industries_.add_record(r)    
    
    #set FinalD
    fd_ = GamsParameter(db,"FinDem_ini", 4, "FinDem")
    for k, v in table_in.FinalD.items():
        fd_.add_record(k).value = v    
    
    #set interaction matrix of intermediate demand
    z_m = db.add_parameter("Z_matrix_ini", 4, "Interaction matrix")
    for k, v in table_in.Z_matrix.items():
        z_m.add_record(k).value = v 
    
    #set interaction matrix of intermediate demand
    a_m = db.add_parameter("A_matrix_ini", 4, "A matrix")
    for k, v in table_in.A_matrix.items():
        a_m.add_record(k).value = v 

    if EORA is not True:
        #set Export ROW
        exp = db.add_parameter("ExpROW_ini", 3, "Exports to ROW")
        for k, v in table_in.ExpROW.items():
            exp.add_record(k).value = v 
    
    #set ValueA
    val = db.add_parameter("ValueA_ini", 3, "Value Added")
    for k, v in table_in.ValueA.items():
        val.add_record(k).value = v 
        
    # And save to GDX file
    db.export(os.path.join(data_path,"gams_runs","{}.gdx".format(table_in.name)))

def ratmarg_SUT(table_in,EORA=False):

    data_path = load_config()['paths']['data']
    
    table_in.prep_data()

    load_db_SUT(table_in)
    
    '''
    RUN SCRIPT WITH DISRUPTION
    '''
    setdir = os.path.join(data_path,'gams_runs')
    ws = GamsWorkspace(setdir)
    ws.get_working_directory()
    
    gamsfile_in =  os.path.join(setdir,"obtain_marg_value_SUT.gms")
    gamsfile = os.path.join(setdir,"obtain_marg_value_SUT_{}.gms".format(table_in.name))
    copyfile(gamsfile_in, gamsfile)
    str_ctry = ','.join(table_in.countries) 
    str_fd = 'FinalD' #','.join(list(table_in.FD_labels['FD'].unique()))


    with open(gamsfile, 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    gdx_file = "%s.gdx" % table_in.name
    data[30] = '$GDXIN '+gdx_file+'\n'

    str_ind = ','.join(table_in.sectors) 
    data[40] ='ind(col) list of industries  /'+str_ind+'/\n'

    data[37] = '/'+str_ctry+'/\n'
    data[39] = '/'+str_ctry+'/\n'
    data[43]  = '/'+str_fd+'/\n'

    with open(gamsfile, 'w') as file:
        file.writelines( data )
    
    gamsfile_run = gamsfile.replace("..\\..\\gams_runs\\", "")
    t1 = ws.add_job_from_file(gamsfile_run)
    
    t1.run()
    
    Ratmarg = []
    index_ = []
    for rec in t1.out_db["Ratmarg"]:
        index_.append((rec.keys[0],rec.keys[1]))
        Ratmarg.append(rec.get_value())     
        
    index_ = pd.MultiIndex.from_tuples(index_, names=('CNTRY', 'IND'))
    Ratmarginal = pd.DataFrame(Ratmarg, index=index_).unstack()
    Ratmarginal.columns = Ratmarginal.columns.droplevel()

    Ratmarginal.to_csv(os.path.join(data_path,'input_data','Ratmarg_{}.csv'.format(table_in.name)))

    return Ratmarginal    

def ratmarg_IO(table_in,EORA=False):

    data_path = load_config()['paths']['data']

    table_in.prep_data()

    if EORA is True:
        load_db_IO(table_in,EORA=True)
    else:
        load_db_IO(table_in)
    
    '''
    RUN SCRIPT WITH DISRUPTION
    '''
    setdir = os.path.join(data_path,'gams_runs')
    ws = GamsWorkspace(setdir)
    ws.get_working_directory()
    
    if EORA is False:
        gamsfile_in =  os.path.join(data_path,"gams_runs","obtain_marg_value.gms")
        gamsfile = os.path.join(data_path,"gams_runs","obtain_marg_value_{}.gms".format(table_in.name))
        copyfile(gamsfile_in, gamsfile)
        str_ctry = ','.join(table_in.countries)    
        str_fd = ','.join(list(table_in.FD_labels['tfd'].unique()))

    else:
        gamsfile_in =  os.path.join(data_path,"gams_runs","obtain_marg_value_EORA.gms")
        gamsfile = os.path.join(data_path,"gams_runs","obtain_marg_value_{}.gms".format(table_in.name))
        copyfile(gamsfile_in, gamsfile)
        str_ctry = ','.join(table_in.countries+['ROW']) 
        str_fd = ','.join(list(table_in.FD_labels['FD'].unique()))


    with open(gamsfile, 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    gdx_file = "{}.gdx".format(table_in.name)
    data[26] = '$GDXIN '+gdx_file+'\n'

    str_ind = ','.join(table_in.sectors) 
    data[32] ='S(col) list of industries  /'+str_ind+'/\n'


    data[34] = '/'+str_ctry+'/\n'
    data[36] = '/'+str_ctry+'/\n'
    data[38]  = '/'+str_fd+'/\n'

    with open(gamsfile, 'w') as file:
        file.writelines( data )
    
    gamsfile_run = gamsfile.replace("..\\..\\gams_runs\\", "")
    t1 = ws.add_job_from_file(gamsfile_run)
    
    t1.run()
    
    Ratmarg = []
    index_ = []
    for rec in t1.out_db["Ratmarg"]:
        index_.append((rec.keys[0],rec.keys[1]))
        Ratmarg.append(rec.get_value())     
        
    index_ = pd.MultiIndex.from_tuples(index_, names=('CNTRY', 'IND'))
    Ratmarginal = pd.DataFrame(Ratmarg, index=index_).unstack()
    Ratmarginal.columns = Ratmarginal.columns.droplevel()


    Ratmarginal.to_csv(os.path.join(data_path,'input_data','Ratmarg_{}.csv'.format(table_in.name)))


    return Ratmarginal    
