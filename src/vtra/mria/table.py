# -*- coding: utf-8 -*-
"""Create the economic tables required to run the MRIA model.
"""
import numpy as np
import pandas as pd


class io_basic(object):
    """io_basic is used to set up the table.
    """

    def __init__(self, name, filepath, list_regions):
        """Creation of a the object instance, specify the file path and sectors to include.

        Parameters
            - *self* - **io_basic** class object
            - name - string name for the **io_basic** class
            - filepath - string path name to location of IO table
            - list_regions - list of regions to include

        Output
            - *self*.name - string name of the model in the **io_basic** class
            - *self*.file - filepath in the **MRIA_IO** class
            - *self*.regions - list of regions in the **MRIA_IO** class
            - *self*.total_regions - Integer of total amount of regions in the **io_basic** class

        """
        self.name = name
        self.file = filepath
        self.regions = list_regions
        self.total_regions = len(list_regions)

    def load_labels(self):
        """Load all labels for the **io_basic** class.

        Parameters
            - *self* - **io_basic** class object

        Output
            - *self*.FD_labels - labels for Final Demand columns in the **io_basic** class
            - *self*.FD_cat - labels for Final Demand categories in the **io_basic** class
            - *self*.Exp_labels - labels for Export columns in the **io_basic** class
            - *self*.T_labels - region and sector labels for Z-matrix  in the **io_basic** class
            - *self*.VA_labels - labels for Value Added in the **io_basic** class
            - *self*.sectors - labels for the sectors in the **io_basic** class

        """

        if 'xls' in self.file:
            FD_labels = pd.read_excel(self.file, sheet_name="labels_FD",
                                      names=['reg', 'tfd'], header=None)
            Exp_labels = pd.read_excel(self.file, sheet_name="labels_ExpROW", names=[
                                       'export'], header=None)
            T_labels = pd.read_excel(self.file, sheet_name="labels_T",
                                     header=None, names=['reg', 'ind'])
            VA_labels = pd.read_excel(self.file, sheet_name="labels_VA", names=[
                                      'Import', 'ValueA'], header=None)

        if len(self.regions) == 0:
            self.regions = list(T_labels['reg'].unique())
            self.total_regions = len(self.regions)

        self.FD_labels = FD_labels
        self.FD_cat = list(self.FD_labels['tfd'].unique())
        self.Exp_labels = Exp_labels
        self.T_labels = T_labels
        self.VA_labels = VA_labels
        self.sectors = list(T_labels['ind'].unique())

    def load_all_data(self):
        """Load all data for the **io_basic** class.

        Parameters
            - *self* - **io_basic** class object

        Output
            - *self*.FD_data - pandas Dataframe of Final Demand in the **io_basic** class
            - *self*.T_data - pandas Dataframe of Z matrix in the **io_basic** class
            - *self*.VA_data - pandas Dataframe of Value Added in the **io_basic** class
            - *self*.ImpROW_data - pandas Dataframe of import from the Rest of the World in the **io_basic** class
            - *self*.ExpROW_data - pandas Dataframe of exports to the Rest of The World in the **io_basic** class

        """
        try:
            self.FD_labels is None
        except:
            self.load_labels()

        #LOAD DATA
        FD_data = pd.read_excel(self.file, sheet_name="FD", header=None)
        T_data = pd.read_excel(self.file, sheet_name="T", header=None)
        VA_data = pd.read_excel(self.file, sheet_name="VA", header=None)
        ExpROW_data = pd.read_excel(self.file, sheet_name="ExpROW", header=None)

        # Add labels to the data from 'load_labels'
        FD_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        ExpROW_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        T_data.index = pd.MultiIndex.from_arrays(self.T_labels.values.T)

        reg_label = np.array(
            list(self.T_labels.values.T[0])+list(self.FD_labels.values.T[0])+['export'])
        ind_label = np.array(
            list(self.T_labels.values.T[1])+list(self.FD_labels.values.T[1])+['export'])
        va_index = np.vstack((reg_label, ind_label))

        VA_data.index = pd.MultiIndex.from_arrays(va_index)

        FD_data.columns = pd.MultiIndex.from_arrays(self.FD_labels.values.T)
        ExpROW_data.columns = pd.MultiIndex.from_arrays(self.Exp_labels.values.T)
        T_data.columns = pd.MultiIndex.from_arrays(self.T_labels.values.T)
        VA_data.columns = pd.MultiIndex.from_arrays(self.VA_labels.values)


        # And return the data to the mother class
        self.FD_data = FD_data
        self.T_data = T_data
        self.VA_data = pd.DataFrame(VA_data['VA'])
        self.ImpROW_data = pd.DataFrame(VA_data['Import'])
        self.ExpROW_data = ExpROW_data

    def prep_data(self):
        """Transform the dataframes into dictionaries, ready to be used in the **MRIA_IO** class instance.

        Parameters
            - *self* - **io_basic** class object

        Output
            - *self*.FinalD - dictionary of Final Demand in the **io_basic** class
            - *self*.A_matrix - dictionary of A matrix in the **io_basic** class
            - *self*.Z_matrix - dictionary of Z matrix in the **io_basic** class
            - *self*.ValueA - dictionary of Value Added in the **io_basic** class
            - *self*.ImpROW - dictionary of import from the Rest of the World in the **io_basic** class
            - *self*.ExpROW - dictionary of exports to the Rest of The World in the **io_basic** class

        """
        try:
            self.FD_data is None
        except:
            self.load_all_data()

        self.sum_data = self.T_data.sum(
            axis=1)+self.FD_data.sum(axis=1)+self.ExpROW_data.sum(axis=1)

        self.A = self.T_data.divide(self.sum_data, axis=1)


        #Return all the parts of the dataset to the class again

        self.Z_matrix = {r + k: v for r, kv in self.T_data.iterrows()
                         for k, v in kv.to_dict().items()}
        self.A_matrix = {r + k: v for r, kv in self.A.iterrows()
                         for k, v in kv.to_dict().items()}
        self.FinalD = {r + k: v for r, kv in self.FD_data.iterrows()
                       for k, v in kv.to_dict().items()}
        self.ValueA = {r + k: v for r, kv in self.VA_data.iterrows()
                       for k, v in kv.to_dict().items()}
        self.ImpROW = {r + k: v for r, kv in self.ImpROW_data.iterrows()
                       for k, v in kv.to_dict().items()}
        self.ExpROW = {r + k: v for r, kv in self.ExpROW_data.iterrows()
                       for k, v in kv.to_dict().items()}
