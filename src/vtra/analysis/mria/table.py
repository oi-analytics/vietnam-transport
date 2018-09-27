# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 15:26:37 2017

@author: Elco Koks

@date: Nov, 2017
"""

import pandas as pd
import numpy as np

class io_basic(object):

    """
    This is the class object 'EORA' which is used to set up the table.
    """
  
    def __init__(self, name,filepath,year,list_countries=None):
        
        self.year = year
        self.name = name
        self.file = filepath
        if list_countries is not None:
            self.countries = list_countries
            self.total_countries = len(list_countries)
        else:
            self.countries = []
            self.total_countries = 0


    def load_labels(self):

        """
        LOAD LABELS
        """

        if 'xls' in self.file:
            FD_labels = pd.read_excel(self.file,sheet_name="labels_FD",names=['reg','tfd'],header=None)
            Exp_labels = pd.read_excel(self.file,sheet_name="labels_ExpROW",names=['export'],header=None)
            T_labels = pd.read_excel(self.file,sheet_name="labels_T",header=None,names=['reg','ind'])
            VA_labels = pd.read_excel(self.file,sheet_name="labels_VA",names=['Import','ValueA'],header=None)


        if len(self.countries) == 0:
            self.countries = list(T_labels['reg'].unique())
            self.total_countries = len(self.countries)   
 
        self.FD_labels = FD_labels
        self.FD_cat = list(self.FD_labels['tfd'].unique())
        self.Exp_labels = Exp_labels
        self.T_labels = T_labels
        self.VA_labels = VA_labels
        self.sectors = list(T_labels['ind'].unique())


    def load_all_data(self):
        
        try: 
            self.FD_labels is None
        except:
            self.load_labels()

        """
        LOAD DATA
        """
        FD_data = pd.read_excel(self.file,sheet_name="FD",header=None)
        T_data  = pd.read_excel(self.file,sheet_name="T",header=None)
        VA_data = pd.read_excel(self.file,sheet_name="VA",header=None)
        ExpROW_data = pd.read_excel(self.file,sheet_name="ExpROW",header=None)

        """
        Add labels to the data from 'load_labels'
        """
        FD_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        ExpROW_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        T_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)

        reg_label = np.array(list(self.T_labels.values.T[0])+list(self.FD_labels.values.T[0])+['export'])
        ind_label = np.array(list(self.T_labels.values.T[1])+list(self.FD_labels.values.T[1])+['export'])
        va_index = np.vstack((reg_label,ind_label))
    
        VA_data.index = pd.MultiIndex.from_arrays(va_index)
        
        FD_data.columns = pd.MultiIndex.from_arrays(self.FD_labels.values.T)
        ExpROW_data.columns = pd.MultiIndex.from_arrays(self.Exp_labels.values.T)
        T_data.columns = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        VA_data.columns= pd.MultiIndex.from_arrays(self.VA_labels.values)
        
        """
        And return the data to the mother class
        """
        self.FD_data = FD_data
        self.T_data = T_data
        self.VA_data = pd.DataFrame(VA_data['VA'])
        self.ImpROW_data = pd.DataFrame(VA_data['Import'])
        self.ExpROW_data = ExpROW_data


    def prep_data(self):

        try: 
            self.FD_data is None
        except:
            self.load_all_data()
        
        self.sum_data = self.T_data.sum(axis=1)+self.FD_data.sum(axis=1)+self.ExpROW_data.sum(axis=1)

        self.A = self.T_data.divide(self.sum_data,axis=1)
        
        """
        Return all the parts of the dataset to the class again
        """
        
        self.Z_matrix = {r + k: v for r, kv in self.T_data.iterrows() for k,v in kv.to_dict().items()}
        self.A_matrix = {r + k: v for r, kv in self.A.iterrows() for k,v in kv.to_dict().items()}
        self.FinalD = {r + k: v for r, kv in self.FD_data.iterrows() for k,v in kv.to_dict().items()}
        self.ValueA = {r + k: v for r, kv in self.VA_data.iterrows() for k,v in kv.to_dict().items()}
        self.ImpROW = {r + k: v for r, kv in self.ImpROW_data.iterrows() for k,v in kv.to_dict().items()}
        self.ExpROW = {r + k: v for r, kv in self.ExpROW_data.iterrows() for k,v in kv.to_dict().items()}


class io_eora(object):

    """
    This is the class object 'EORA' which is used to set up the table.
    """
   
    def __init__(self, name, year,list_countries=None):
        
        self.year = year
        self.name = name
        if list_countries is not None:
            self.countries = list_countries
            self.total_countries = len(list_countries)
        else:
            self.countries = []
            self.total_countries = 0


    def load_labels(self):

        """
        LOAD LABELS
        """
    
        FD_labels = pd.read_csv('..\..\input_data\labels_FD.txt', sep='\t',index_col=False, 
                                names=['COUNTRY','CNTRY_CODE','Final Demand','FD'])
        T_labels = pd.read_csv('..\..\input_data\labels_T.txt', sep='\t',index_col=False,
                               names=['COUNTRY','CNTRY_CODE','Industries','IND'])
        VA_labels = pd.read_csv('..\..\input_data\labels_VA.txt', sep='\t',index_col=False,
                                names=['Value Added','VA'])
        

        if len(self.countries) == 0:
            self.countries = list(T_labels['COUNTRY'].unique())
            self.total_countries = len(self.countries)   

        """
        Convert labels to more usable format
        """

        # For Final Demand
        new_FD = ['P3h','P3n','P3g','P51','P52','P53']
        FD_to_replace = dict(zip(FD_labels['FD'].unique(),new_FD))
    
        FD_labels.replace(FD_to_replace,inplace=True)

        # For each of the sectors        
        new_T = ['i'+str(n+1) for n in range(26)]
        T_to_replace = dict(zip(T_labels['IND'].unique(),new_T))
    
        T_labels.replace(T_to_replace,inplace=True)

        self.FD_labels = FD_labels
        self.T_labels = T_labels
        self.VA_labels = VA_labels
        self.sectors = new_T


    def load_all_data(self,data_path=None):
        
        try: 
            self.FD_labels is None
        except:
            self.load_labels()

        """
        LOAD DATA
        """
        FD_data = pd.read_csv('..\..\input_data\Eora26_%s_bp_FD.txt' % self.year, sep='\t',index_col=False,header=None)/1000
        T_data = pd.read_csv('..\..\input_data\Eora26_%s_bp_T.txt' % self.year, sep='\t',index_col=False,header=None)/1000
        VA_data = pd.read_csv('..\..\input_data\Eora26_%s_bp_VA.txt' % self.year, sep='\t',index_col=False,header=None)/1000

        """
        Add labels to the data from 'load_labels'
        """
        FD_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T).droplevel(0).droplevel(1)
        T_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T).droplevel(0).droplevel(1)
        VA_data.index = pd.MultiIndex.from_arrays(self.VA_labels.values.T).droplevel(0)
        
        FD_data.columns = pd.MultiIndex.from_arrays(self.FD_labels.values.T).droplevel(0).droplevel(1)
        T_data.columns = pd.MultiIndex.from_arrays(self.T_labels.values.T).droplevel(0).droplevel(1)
        VA_data.columns= pd.MultiIndex.from_arrays(self.T_labels.values.T).droplevel(0).droplevel(1)
        
        """
        And return the data to the mother class
        """
        self.FD_data = FD_data
        self.T_data = T_data
        self.VA_data = VA_data


    def prep_data(self,curdir=None):
        
        try: 
            self.FD_data is None
        except:
            self.load_all_data(curdir)
        
        """
        Create subsets for the 'list_countries'
        """

        # create transaction matrix for the subset 
        subset_T = self.T_data[self.T_data.index.get_level_values(0).isin(self.countries)]
        subset_T = subset_T.iloc[:, subset_T.columns.get_level_values(0).isin(self.countries)]
        subset_T.index.names = ['CNTRY', 'IND']
        
        # create final demand matrix for the subset
        subset_FD = self.FD_data[self.FD_data.index.get_level_values(0).isin(self.countries)]                    
        subset_FD = subset_FD.iloc[:, subset_FD.columns.get_level_values(0).isin(self.countries)]
        subset_FD.index.names = ['CNTRY', 'IND']
        
        # create value added matrix for the subset
        subset_VA = self.VA_data.iloc[:, self.VA_data.columns.get_level_values(0).isin(self.countries)]
        
        # create aggregated transaction matrix for the rest of the world
        T_row = self.T_data.iloc[~self.T_data.columns.get_level_values(0).isin(self.countries),~self.T_data.index.get_level_values(0).isin(self.countries)]
        T_row = T_row.groupby(level=[1],axis=0).sum().groupby(level=[1],axis=1).sum()
        
        T_row['ROW'] = 'ROW'
        T_row = T_row.set_index('ROW',append=True,drop=True)
        T_row = T_row.reorder_levels([1,0], axis=0)
        T_row.columns = pd.MultiIndex.from_product([['ROW'],T_row.columns])
        T_row.index.names = ['CNTRY', 'IND']
        
        # create aggregated transaction import matrix for the rest of the world
        ImpTrade_ROW = self.T_data.iloc[~self.T_data.index.get_level_values(0).isin(self.countries),self.T_data.columns.get_level_values(0).isin(self.countries)]
        ImpTrade_ROW= ImpTrade_ROW.groupby(level=[1],axis=0).sum()
        
        ImpTrade_ROW[('ROW','ROW')] = 'ROW'
        ImpTrade_ROW = ImpTrade_ROW.set_index(('ROW','ROW'),append=True,drop=True)
        ImpTrade_ROW = ImpTrade_ROW.reorder_levels([1,0], axis=0)
        ImpTrade_ROW.index.names = ['CNTRY', 'IND']
        
        # create aggregated transaction export matrix for the rest of the world
        ExpTrade_ROW = self.T_data.iloc[self.T_data.index.get_level_values(0).isin(self.countries),~self.T_data.columns.get_level_values(0).isin(self.countries)]
        ExpTrade_ROW= ExpTrade_ROW.groupby(level=[1],axis=1).sum()
        ExpTrade_ROW.columns = pd.MultiIndex.from_product([['ROW'],ExpTrade_ROW.columns])
        ExpTrade_ROW.index.names = ['CNTRY', 'IND']
        
        # create aggregated final demand matrix for the rest of the world
        FD_row = self.FD_data.iloc[~self.FD_data.index.get_level_values(0).isin(self.countries), ~self.FD_data.columns.get_level_values(0).isin(self.countries)]
        FD_row = FD_row.groupby(level=[1],axis=0).sum().groupby(level=[1],axis=1).sum()
        
        FD_row['ROW'] = 'ROW'
        FD_row = FD_row.set_index('ROW',append=True,drop=True)
        FD_row = FD_row.reorder_levels([1,0], axis=0)
        FD_row.columns = pd.MultiIndex.from_product([['ROW'],FD_row.columns])
        FD_row.index.names = ['CNTRY', 'IND']
        
        # create aggregated final demand import matrix for the rest of the world
        ImpFD_ROW = self.FD_data.iloc[~self.FD_data.index.get_level_values(0).isin(self.countries),self.FD_data.columns.get_level_values(0).isin(self.countries)]
        ImpFD_ROW= ImpFD_ROW.groupby(level=[1],axis=0).sum()
        
        ImpFD_ROW[('ROW','ROW')] = 'ROW'
        ImpFD_ROW = ImpFD_ROW.set_index(('ROW','ROW'),append=True,drop=True)
        ImpFD_ROW = ImpFD_ROW.reorder_levels([1,0], axis=0)
        ImpFD_ROW.index.names = ['CNTRY', 'IND']
        
        # create aggregated final demand export matrix for the rest of the world
        ExpFD_ROW = self.FD_data.iloc[self.FD_data.index.get_level_values(0).isin(self.countries),~self.FD_data.columns.get_level_values(0).isin(self.countries)]
        ExpFD_ROW= ExpFD_ROW.groupby(level=[1],axis=1).sum()
        ExpFD_ROW.columns = pd.MultiIndex.from_product([['ROW'],ExpFD_ROW.columns])
        ExpFD_ROW.index.names = ['CNTRY', 'IND']
        
        # create aggregated value added matrix for the rest of the world
        VA_row = self.VA_data.iloc[:,~ self.VA_data.columns.get_level_values(0).isin(self.countries)]
        VA_row = VA_row.groupby(level=[1],axis=1).sum()
        VA_row.columns = pd.MultiIndex.from_product([['ROW'],VA_row.columns])
        
        # combine matrices to have a new transaction matrix
        df1 = ImpTrade_ROW.reset_index().join(T_row,on=['CNTRY', 'IND']).set_index(ImpTrade_ROW.index.names)
        df2 = pd.merge(subset_T.reset_index(), ExpTrade_ROW.reset_index(), on=['CNTRY', 'IND'], how='inner').set_index(['CNTRY', 'IND'])
        subset_T = pd.concat([df1,df2])
    
       
        # combine matrices to have a new final demand matrix
        df1 = ImpFD_ROW.reset_index().join(FD_row,on=['CNTRY', 'IND']).set_index(ImpFD_ROW.index.names)
        df2 = pd.merge(subset_FD.reset_index(), ExpFD_ROW.reset_index(), on=['CNTRY', 'IND'], how='inner').set_index(['CNTRY', 'IND'])
        subset_FD = pd.concat([df1,df2])
        
        # combine matrices to have a new Value Added matrix
        subset_VA = subset_VA.join(VA_row)
        sum_subset = subset_T.sum(axis=1) + subset_FD.sum(axis=1)             
        sum_outlays =  subset_VA.sum()   + subset_T.sum(axis=0)     
        sum_outlays.index.names = ['CNTRY', 'IND']
        
        VA_TZAetal = pd.DataFrame((subset_VA.sum() + (sum_subset-sum_outlays)).rename("VA"))            


        subset_T.columns.names = ['CNTRY', 'IND']
        self.A = subset_T.divide(sum_subset,axis=1)
        
        self.A = self.A.fillna(0)

        self.A[self.A == 0] = 1e-10
    
        """
        Return all the parts of the dataset to the class again
        """
        
        self.Z_matrix = {r + k: v for r, kv in subset_T.iterrows() for k,v in kv.to_dict().items()}
        self.A_matrix = {r + k: v for r, kv in self.A.iterrows() for k,v in kv.to_dict().items()}
        self.FinalD = {r + k: v for r, kv in subset_FD.iterrows() for k,v in kv.to_dict().items()}
        self.ValueA = {r + (k,): v for r, kv in VA_TZAetal.iterrows() for k,v in kv.to_dict().items()}

class io_oecd(object):

    """
    This is the class object 'EORA' which is used to set up the table.
    """
   
    def __init__(self, name, filepath, year,list_countries=None):
        
        self.year = year
        self.name = name
        self.file = filepath
        if list_countries is not None:
            self.countries = list_countries
            self.total_countries = len(list_countries)
        else:
            self.countries = []
            self.total_countries = 0

    def load_labels(self):

        """
        LOAD LABELS
        """

        if 'csv' in self.file:
            FD_labels = pd.read_excel('..\..\input_data\ICIO_2016_2011_labels.xlsx',sheetname="labels_FD",names=['reg','tfd'],header=None)
            Exp_labels = pd.read_excel('..\..\input_data\ICIO_2016_2011_labels.xlsx',sheetname="labels_ExpROW",names=['export'],header=None)
            T_labels = pd.read_excel('..\..\input_data\ICIO_2016_2011_labels.xlsx',sheetname="labels_T",header=None,names=['reg','ind','reg2'])
            VA_labels = pd.read_excel('..\..\input_data\ICIO_2016_2011_labels.xlsx',sheetname="labels_VA",names=['ValueA'],header=None)

        """
        Return labels to the mother class
        """

        if len(self.countries) == 0:
            self.countries = list(T_labels['reg'].unique())
            self.total_countries = len(self.countries)   

        
        self.FD_labels = FD_labels
        self.Exp_labels = Exp_labels
        self.T_labels = T_labels
        self.VA_labels = VA_labels
        self.sectors = list(T_labels['ind'].unique())
        

    def load_all_data(self):

        try: 
            self.FD_labels is None
        except:
            self.load_labels()

        
        oecd_raw = pd.read_csv(self.file,header=None)

        """
        LOAD DATA
        """
        FD_data = oecd_raw.iloc[:2414,2414:-1]
        T_data  = oecd_raw.iloc[:2414,:2414]
        VA_data = pd.DataFrame(oecd_raw.iloc[-2,:2414])
        ExpROW_data = pd.DataFrame(oecd_raw.iloc[:2414,-1])

        """
        Add labels to the data from 'load_labels'
        """
        FD_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        ExpROW_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        T_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        VA_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        
        FD_data.columns = pd.MultiIndex.from_arrays(self.FD_labels.values.T)
        ExpROW_data.columns = pd.MultiIndex.from_arrays(self.Exp_labels.values.T)
        T_data.columns = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        VA_data.columns= pd.MultiIndex.from_arrays(self.VA_labels.values)

        """
        Aggregate mexico and china, we don't care too much about the split ups
        """
        inb = T_data.groupby(level=[0,1],axis=0).sum()
        T_data = inb.groupby(level=[0,1],axis=1).sum()        

        FD_data = FD_data.groupby(level=[0,1],axis=0).sum()
        ExpROW_data = ExpROW_data.groupby(level=[0,1],axis=0).sum()
        VA_data = VA_data.groupby(level=[0,1],axis=0).sum() 

        """
        And return the data to the mother class
        """
        self.FD_data = FD_data
        self.T_data = T_data
        self.VA_data = pd.DataFrame(VA_data['VA'])
        self.ExpROW_data = ExpROW_data

                
    def prep_data(self):

        try: 
            self.FD_data is None
        except:
            self.load_all_data()
        
        self.sum_data = self.T_data.sum(axis=1)+self.FD_data.sum(axis=1)+self.ExpROW_data.sum(axis=1)

        self.A = self.T_data.divide(self.sum_data,axis=1)
        
        self.A = self.A.fillna(0)
    
        """
        Return all the parts of the dataset to the class again
        """
        
        self.Z_matrix = {r + k: v for r, kv in self.T_data.iterrows() for k,v in kv.to_dict().items()}
        self.A_matrix = {r + k: v for r, kv in self.A.iterrows() for k,v in kv.to_dict().items()}
        self.FinalD = {r + k: v for r, kv in self.FD_data.iterrows() for k,v in kv.to_dict().items()}
        self.ValueA = {r + (k,): v for r, kv in self.VA_data.iterrows() for k,v in kv.to_dict().items()}
        self.ImpROW = {r + ('Import',): 0 for r, kv in self.VA_data.iterrows() for k,v in kv.to_dict().items()}
        self.ExpROW = {r + (k,): v for r, kv in self.ExpROW_data.iterrows() for k,v in kv.to_dict().items()}


class sut_basic(object):

    """
    This is the class object 'EORA' which is used to set up the table.
    """
   
    def __init__(self, name,filepath,year,list_countries=None):
        
        self.year = year
        self.name = name
        self.file = filepath
        if list_countries is not None:
            self.countries = list_countries
            self.total_countries = len(list_countries)
        else:
            self.countries = []
            self.total_countries = 0



    def load_all_data(self):
        
        """
        LOAD DATA
        """
        self.Use_data = pd.read_excel(self.file,sheetname="USE",index_col=[0,1],header=[0,1])
        self.Sup_data  = pd.read_excel(self.file,sheetname="SUP",index_col=[0,1],header=[0,1])
        self.VA_data = pd.read_excel(self.file,sheetname="VA",index_col=[0,1],header=[0])
        self.ExpROW_data = pd.read_excel(self.file,sheetname="ExpROW",index_col=[0,1],header=[0])
        self.ImpROW_data = pd.read_excel(self.file,sheetname="ImpROW",index_col=[0,1],header=[0])

        """
        Extract indices
        """
        self.countries = list(set(list(self.Use_data.index.get_level_values(0))))
        self.sectors = list(set(list(self.Sup_data.index.get_level_values(1))))
        self.products = list(set(list(self.Use_data.index.get_level_values(1))))
        self.FD_cat = ['FinalD']

    def prep_data(self):

        try: 
            self.FD_data is None
        except:
            self.load_all_data()
   
        """
        Return all the parts of the dataset to the class again
        """
        self.Use = {r + k: v for r, kv in self.Use_data.iterrows() for k,v in kv.to_dict().items()}
        self.Sup = {r + k: v for r, kv in self.Sup_data.iterrows() for k,v in kv.to_dict().items()}
        self.ValueA = {r + (k,): v for r, kv in self.VA_data.iterrows() for k,v in kv.to_dict().items()}
        self.ImpROW = {r + (k,): v for r, kv in self.ImpROW_data.iterrows() for k,v in kv.to_dict().items()}
        self.ExpROW = {r + (k,): v for r, kv in self.ExpROW_data.iterrows() for k,v in kv.to_dict().items()}

        

