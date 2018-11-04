# -*- coding: utf-8 -*-
"""MRIA Model

Purpose
-------

The Multiregional Impact Assessment (MRIA) Model allows for estimating a new post-disaster economic situation in equilibrium, given a set of disruptions.

References
----------

1) Koks, E. E., & Thissen, M. (2016). A multiregional impact assessment model for disaster analysis. Economic Systems Research, 28(4), 429-449.

"""
import os

import numpy as np
import pandas as pd
from pyomo.environ import (ConcreteModel, Constraint, Objective, Param, Set,
                           SetOf, Var, minimize)
from pyomo.opt import SolverFactory
from vtra.mria.ratmarg import ratmarg_IO
from vtra.utils import load_config


class MRIA_IO(object):
    """MRIA_IO sets up the modelling framework.

    In this class we define the type of model, sets, set up the core variables and specify the
    constraints and objectives for different model setups.

    Attributes
    ----------
    name : string
        name of the model
    m : Pyomo.ConcreteModel
        model instance
    regions : list
        model regions
    total_regions : int
        number of regions
    sectors : list
        model sectors
    fd_cat : list
        final demand categories
    """

    def __init__(self, name, list_regions, list_sectors, list_fd_cats=[]):
        """Create a Concrete Model, specify the regions and sectors to include.

        Parameters
        ----------
        name : str
            name of the model
        list_regions : list
            regions to include in the model calculations
        list_sectors : list
            sectors to include in the model calculations
        list_fd_cats : list
            the final demand categories in the table. This will be aggregated to one column.

        """
        self.name = name
        self.m = ConcreteModel()
        self.regions = list_regions
        self.total_regions = len(list_regions)
        self.sectors = list_sectors
        self.fd_cat = list_fd_cats


    def create_sets(self, FD_SET=None, VA_SET=None):
        """Create of the various sets, allowing for specification of set inputs.

        Parameters
        ----------
        FD_SET : list or None, optional
            final demand categories, default to an **empty list**
        VA_SET : list or None, optional
            value added categories, default to an **empty list**

        """

        if FD_SET is None:
            FD_SET = []

        if VA_SET is None:
            VA_SET = []

        self.m.S = Set(initialize=self.sectors, doc='sectors')
        self.m.rROW = Set(initialize=self.regions, ordered=True,
                          doc='regions including export')
        self.m.R = Set(initialize=self.regions, ordered=True, doc='regions')
        self.m.fdemand = Set(initialize=self.fd_cat, doc='Final Demand')
        self.m.VA = Set(initialize=VA_SET, doc='value added')

    def create_alias(self):
        """Set aliases.

        Parameters
        ----------
        list_regions : list[str]
            regions to include in the model calculations
        list_sectors : list[str]
            sectors to include in the model calculations
        list_fd_cats : list[str]
            the final demand categories in the table. This will be aggregated to one column.

        """
        self.m.Rb = SetOf(self.m.R)  # an alias of region R
        self.m.r = SetOf(self.m.R)  # an alias of region R
        self.m.Sb = SetOf(self.m.S)  # an alias of sector S

    def create_A_mat(self, A_mat_in):
        """Creation of the A-matrix for the optimization model.

        Parameters
        ----------
        A_mat_in
            A-matrix dictionary from the **io_basic** class object

        """
        model = self.m

        def A_matrix_init(model, R, S, Rb, Sb):
            return A_mat_in[R, S, Rb, Sb]

        model.A_matrix = Param(model.R, model.S, model.R, model.Sb,
                               initialize=A_matrix_init, doc='A matrix')

        self.A_matrix = model.A_matrix


    def create_FD(self, FinalD, disr_dict_fd):
        """Specify Final Demand and Local Final Demand.

        Parameters
        ----------
        FinalD
            Final Demand dictionary from the **io_basic** class object
        disr_dict_fd
            dictionary containing the disruptions in final demand

        Create *self*.ttfd - Pyomo Parameter instance for the total final demand
        and *self*.fd - Pyomo Parameter instance for the final demand

        """

        disrupted_des = list(np.unique([x[1] for x in disr_dict_fd]))

        model = self.m

        model.Rdes = Set(initialize=disrupted_des, doc='Final Demand')

        def tfd_init(model, R, S, Rb):
            if (R, Rb, S) in list(disr_dict_fd.keys()):
                return sum(FinalD[R, S, Rb, fdemand] for fdemand in model.fdemand)*disr_dict_fd[R, Rb, S]
            else:
                return sum(FinalD[R, S, Rb, fdemand] for fdemand in model.fdemand)

        def fd_init(model, R, S):
            return sum(model.tfd[R, S, Rb] for Rb in model.Rb)

        model.tfd = Param(model.R, model.S, model.Rb, initialize=tfd_init, doc='Final Demand')

        model.fd = Param(model.R, model.S, initialize=fd_init, doc='Final Demand')

        self.ttfd = model.tfd
        self.fd = model.fd


    def create_LFD(self, FinalD):
        """Specify local final demand

        Parameters
        ----------
        FinalD
            Final Demand dictionary from the **io_basic** class object

        Create *self*.lfd - Pyomo Parameter instance for the local final demand

        """
        model = self.m

        def lfd_init(model, R, S):
            return sum(FinalD[R, S, R, fdemand] for fdemand in model.fdemand)
        model.lfd = Param(model.R, model.S, initialize=lfd_init, doc='Final Demand')

        self.lfd = model.lfd


    def create_ExpImp(self, ExpROW, ImpROW):
        """Specify export and import to the rest of the world

        Parameters
        ----------
        ExpROW
            Exports to the Rest of the World dictionary from the **io_basic** class object
        ImpROW
            Imports from the Rest of the World dictionary from the **io_basic** class object

        Create *self*.ExpROW - Pyomo Parameter instance for the Exports to the Rest of the World
        and *self*.ImpROW - Pyomo Parameter instance for the Imports from the Rest of the World

        """
        model = self.m

        # Specify Export ROW
        def ExpROW_ini(m, R, S):
            return (ExpROW[R, S, 'Export'])
        model.ExpROW = Param(model.R, model.S, initialize=ExpROW_ini,
                             doc='Exports to the rest of the world')

        # Specify Import ROW
        def ImpROW_init(m, R, S):
            return (ImpROW[R, S, 'Import'])
        model.ImpROW = Param(model.R, model.S, initialize=ImpROW_init,
                             doc='Imports from the rest of the world')

        self.ExpROW = model.ExpROW
        self.ImpROW = model.ImpROW

    def create_X_up(self, disr_dict, Regmaxcap=0.98):
        """Specify upper bound of total production **X**.

        Parameters
        ----------
        disr_dict : dict
            dictionary containing the reduction in production capacity
        Regmaxcap : float, optional
            maximum regional capacity. The default value is set to **0.98**

        Create *self*.X_up - Pyomo Parameter instance for the upper bound of total production **X**

        """
        model = self.m

        def shock_init(model, R, S):
            if (R, S) in list(disr_dict.keys()):
                return disr_dict[R, S]
            else:
                return 1.05

        model.X_up = Param(model.R, model.S, initialize=shock_init,
                           doc='Maximum production capacity')
        self.X_up = model.X_up

    def create_Xbase(self, Z_matrix, disr_dict, FinalD=None):
        """Specify Baseline value of total production **X**

        Parameters
        ----------
        Z_matrix
            Z-matrix dictionary from the **io_basic** class object
        disr_dict
            dictionary containing the reduction in production capacity
        FinalD
            Final Demand dictionary from the **io_basic** class object

        Create *self*.X_up - Pyomo Parameter instance for the Baseline value of total production **X**

        """
        model = self.m

        if self.fd.active is not True:
            self.create_FD(FinalD, disr_dict)

        if self.ExpROW.active is not True:
            self.create_ExpImp(Z_matrix)

        def x_init_base(model, R, S):
            return(sum(Z_matrix[R, S, Rb, Sb] for Rb in model.Rb for Sb in
                       model.Sb) + self.fd[R, S] + self.ExpROW[R, S])

        model.Xbase = Param(model.R, model.S, initialize=x_init_base,
                            doc='Total Production baseline')
        self.Xbase = model.Xbase

    def create_X(self, disr_dict, Regmaxcap=0.98, A_matrix_ini=None, Z_matrix=None,
                 FinalD=None, Xbase=None, fd=None, ExpROW=None):
        """Create the total production **X** variable

        Parameters
        ----------
        disr_dict
            dictionary containing the reduction in production capacity
        Regmaxcap
            maximum regional capacity. The default value is set to **0.98**
        A_matrix_ini
             A-matrix dictionary from the **io_basic** class object
        Z_matrix
            Z-matrix dictionary from the **io_basic** class object
        FinalD
            Final Demand dictionary from the **io_basic** class object
        Xbase
            Total Production **X** parameter from the **MRIA** class object
        fd
            Final Demand parameter from the **MRIA** class object
        ExpROW
            Export to the Rest of the World parameter from the **MRIA** class object

        Create *self*.X_up - Pyomo Variable instance of total production **X**

        """

        model = self.m

        if self.Xbase.active is not True:
            self.create_Xbase(Z_matrix, FinalD)

        if self.A_matrix.active is not True:
            self.create_A_mat(A_matrix_ini)

        def X_bounds(model, R, S):
            if (R, S) in list(disr_dict.keys()):
                return (0.0, (1/Regmaxcap*self.Xbase[R, S])*disr_dict[R, S])
            else:
                return (0.0, (1/Regmaxcap*self.Xbase[R, S])*1.1)

        def x_init(model, R, S):
            return(sum(self.A_matrix[R, S, Rb, Sb]*self.Xbase[Rb, Sb] for Rb in model.Rb
                       for Sb in model.Sb) + self.fd[R, S] + self.ExpROW[R, S])

        model.X = Var(model.R, model.S, bounds=X_bounds,
                      initialize=x_init, doc='Total Production')

        self.X = model.X


    def create_VA(self, ValueA):
        """Specify Value Added

        Parameters
        ----------
        ValueA
            Value Added dictionary from the **io_basic** class object

        Create *self*.ValueA - Pyomo Parameter instance for the total Value Added

        """
        model = self.m

        def va_init(model, R, S):
            return ValueA[R, S, 'VA']

        model.ValueA = Param(model.R, model.S, initialize=va_init, doc='Value Added')

        self.ValueA = model.ValueA

    def create_Z_mat(self):
        """Specify Trade between regions

        Create *self*.Z_matrix - Pyomo Parameter instance for the total trade matrix

        """
        model = self.m

        def Z_matrix_init(model, R, S, Rb, Sb):
            return self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb]

        model.Z_matrix = Param(model.R, model.S, model.R, model.Sb,
                               initialize=Z_matrix_init, doc='Z matrix')
        self.Z_matrix = model.Z_matrix

    def create_Trade(self, FinalD, Z_matrix=None):
        """Create Trade Matrix

        Parameters
        ----------
        FinalD
            Final Demand dictionary from the **io_basic** class object
        Z_matrix
            Z-matrix dictionary from the **io_basic** class object

        Create *self*.trade - Pyomo Parameter instance for the trade matrix between regions

        """
        model = self.m

        def Trade_init(model, R, Rb, S):
            while R != Rb:
                return sum(self.Z_matrix[Rb, S, R, i] for i in model.Sb) + sum(FinalD[Rb, S, R, i] for i in model.fdemand)

        model.trade = Param(model.R, model.Rb, model.S, initialize=Trade_init, doc='Trade')
        self.trade = model.trade


    def create_TotExp(self):
        """Estimate Total Export

        Create *self*.TotExp - Pyomo Parameter instance for the total export

        """
        model = self.m

        def totexp_init(model, R, S):
            return sum(self.trade[Rb, R, S] for Rb in model.Rb if (R != Rb))

        model.TotExp = Param(model.R, model.S, initialize=totexp_init,
                             doc='Total exports between regions')
        self.TotExp = model.TotExp


    def create_TotImp(self):
        """Estimate Total Import

        Create *self*.TotExp - Pyomo Parameter instance for the total import

        """
        model = self.m

        def totimp_init(model, R, S):
            return sum(self.trade[R, Rb, S] for Rb in model.Rb if (R != Rb))

        model.TotImp = Param(model.R, model.S, initialize=totimp_init,
                             doc='Total imports between regions')
        self.TotImp = model.TotImp


    def create_ImpShares(self):
        """Estimate Import shares and Import share DisImp

        Create *self*.ImportShare - Pyomo Parameter instance for the total import shares

        """
        model = self.m

        def impsh_init(model, R, Rb, S):
            # & ((sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb] for Sb in model.Sb) + self.fd[Rb, S]) != None):
            while (self.trade[Rb, R, S] != None):
                try:
                    return self.trade[Rb, R, S]/(sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb] for Sb in model.Sb) + self.fd[Rb, S])
                except ZeroDivisionError:
                    return 0

        def impshdis_init(model, R, Rb, S):
            return (sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb] for Sb in model.Sb) + self.fd[Rb, S])

        model.ImportShare = Param(model.R, model.Rb, model.S,
                                  initialize=impsh_init, doc='Importshare of each region')
        model.ImportShareDisImp = Param(
            model.R, model.Rb, model.S, initialize=impshdis_init, doc='Importshare DisImp of each region')

        self.ImportShare = model.ImportShare
        self.ImportShareDisImp = model.ImportShareDisImp


    def create_Rdem(self):
        """Create reconstruction demand variable.

        Create *self*.Rdem - Pyomo Parameter instance for the total reconstruction demand

        """
        model = self.m
        model.Rdem = Param(model.R, model.S, initialize=0, doc='Reconstruction demand')
        self.Rdem = model.Rdem

    def create_Rat(self, FinalD=None, Z_matrix=None):
        """Create rationing variable

        Parameters
        ----------
        FinalD
            Final Demand dictionary from the **io_basic** class object
        Z_matrix
            Z-matrix dictionary from the **io_basic** class object

        Set *self*.Rat - Pyomo Variable instance for rationing

        """

        model = self.m

        if self.lfd.active is not True:
            self.create_LFD(FinalD)

        if self.ExpROW.active is not True:
            self.create_ExpImp(Z_matrix)

        def Rat_bounds(model, R, S):
            return (0, abs(self.lfd[R, S]+self.ExpROW[R, S]))

        model.Rat = Var(model.R, model.S, bounds=Rat_bounds, initialize=0, doc='Rationing')
        self.Rat = model.Rat

    def create_Ratmarg(self, Table):
        """Estimate the marginal values of the rationing variable

        Parameters
        ----------
        Table
            the **io_basic** class object

        Set *self*.Ratmarg - Pyomo Parameter instance for the marginal value of rationing

        """
        model = self.m

        try:
            data_path = load_config()['paths']['data']
            RatMarg = pd.read_csv(os.path.join(data_path, 'input_data',
                                               'Ratmarg_{}.csv'.format(self.name)), index_col=[0], header=0)

            if (set(list(RatMarg.index.values)) != set(list(self.regions))):
                RatMarg = ratmarg_IO(Table)
        except:
            RatMarg = ratmarg_IO(Table)

        Ratmarginal = {(r, k): v for r, kv in RatMarg.iterrows()
                       for k, v in kv.to_dict().items()}

        model.Ratmarg = Param(model.R, model.S, initialize=Ratmarginal,
                              doc='Rationing marginal', mutable=True)
        self.Ratmarg = model.Ratmarg

    def create_DisImp(self, disr_dict, Regmaxcap=0.98):
        """Create disaster import variable.

        Parameters
        ----------
        disr_dict
            dictionary containing the reduction in production capacity
        Regmaxcap
            maximum regional capacity. The default value is set to **0.98**

        Set *self*.DisImp - Pyomo Variable instance for disaster imports

        """
        model = self.m

        disrupted_ctry = list(np.unique([x[0] for x in disr_dict]))

        # problem regions
        dimp_ctry = ['KEN', 'UGA']
        dimp_ind = ['i3']

        def Dis_bounds(model, R, S):
            if R in dimp_ctry and S in dimp_ind:
                return (0, 0)
            elif (model.X_up[R, S] < (1.05) or R in disrupted_ctry):
                return (0, None)
            else:
                return (0, None)

        model.DisImp = Var(model.R, model.S, bounds=Dis_bounds,
                           initialize=0, doc='Disaster Imports')
        self.DisImp = model.DisImp


    def create_demand(self):
        """Specify demand function

        Create *self*.Demand - Pyomo Variable instance for total demand

        """
        model = self.m

        def demand_init(model, R, S):
            return (
                sum(self.A_matrix[R, S, R, Sb]*self.X[R, Sb]
                    for Sb in model.Sb) + self.lfd[R, S] + self.Rdem[R, S] - self.Rat[R, S]
                + sum(self.ImportShare[R, Rb, S]*(sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb]
                                                      for Sb in model.Sb) + self.fd[Rb, S] + self.Rdem[Rb, S] - self.Rat[R, S]) for Rb in model.Rb if (R != Rb))
                + sum(self.ImportShare[R, Rb, S]*(self.DisImp[Rb, S])
                      for Rb in model.Rb if (R != Rb))
                + self.ExpROW[R, S]
            )

        model.Demand = Var(model.R, model.S, bounds=(0.0, None), initialize=demand_init)
        self.Demand = model.Demand

    #
    # Create baseline dataset to use in model
    #

    def baseline_data(self, Table, disr_dict_sup, disr_dict_fd):
        """Set up all the baseline variables for the MRIA model

        Parameters
        ----------
        Table
            the **io_basic** class object
        disr_dict_sup
            dictionary containing the reduction in production capacity
        disr_dict_fd
            dictionary containing the disruptions in final demand

        Set all required parameters and variables for the **MRIA_IO** class and the **MRIA**
        model.

        """

        self.create_ExpImp(Table.ExpROW, Table.ImpROW)

        self.create_A_mat(Table.A_matrix)
        self.create_FD(Table.FinalD, disr_dict_fd)
        self.create_LFD(Table.FinalD)
        self.create_Xbase(Table.Z_matrix, Table.FinalD, disr_dict_fd)
        self.create_X(disr_dict_sup, Z_matrix=Table.Z_matrix, FinalD=Table.FinalD)
        self.create_VA(Table.ValueA)
        self.create_Z_mat()
        self.create_Trade(Table.FinalD)
        self.create_TotExp()
        self.create_TotImp()
        self.create_ImpShares()

    def impact_data(self, Table, disr_dict_sup, disr_dict_fd, Regmaxcap=0.98):
        """Create additional parameters and variables required for impact analysis

        Parameters
        ----------
        Table
            the **io_basic** class object
        disr_dict_sup
            dictionary containing the reduction in production capacity
        disr_dict_fd
            dictionary containing the disruptions in final demand
        Regmaxcap
            maximum regional capacity. The default value is set to **0.98**

        Set all additional parameters and variables required for the **MRIA_IO** class and the
        **MRIA** model to do an impact analysis.
        """

        self.create_X_up(disr_dict_sup)
        self.create_Rdem()
        self.create_Rat(Table.FinalD, Table.Z_matrix)
        self.create_Ratmarg(Table)
        self.create_DisImp(disr_dict_sup)
        self.create_demand()

    #
    # Set up baseline model
    #

    def run_basemodel(self, solver=None):
        """Run the baseline model of the MRIA model.

        This should return the baseline situation (i.e. no changes between input matrix and
        output matrix).

        Parameters
        ----------
        solver : str
            Specify the solver to be used with Pyomo. The Default value is set to
            **None**. If set to **None**, the ipopt solver will be used

        Write out the output of an optimized **MRIA_IO** class and the **MRIA** model

        """
        model = self.m

        if solver is None:
            solver = 'ipopt'

        def demSup(model, R, S):
            return (self.X[R, S] >=
                    sum(self.A_matrix[R, S, R, Sb]*self.X[R, Sb]
                        for Sb in model.Sb) + self.lfd[R, S]
                    + sum(self.ImportShare[R, Rb, S]*(sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb]
                                                          for Sb in model.Sb) + self.fd[Rb, S]) for Rb in model.Rb if (R != Rb))
                    + self.ExpROW[R, S])

        model.demSup = Constraint(model.R, model.S, rule=demSup, doc='Satisfy demand')

        def objective_base(model):
            return sum(self.X[R, S] for R in model.R for S in model.S)

        model.objective = Objective(rule=objective_base, sense=minimize,
                                    doc='Define objective function')

        opt = SolverFactory(solver)
        if solver is 'ipopt':
            opt.options['warm_start_init_point'] = 'yes'
            opt.options['warm_start_bound_push'] = 1e-6
            opt.options['warm_start_mult_bound_push'] = 1e-6
            opt.options['mu_init'] = 1e-6
        results = opt.solve(model, tee=True)
        # sends results to stdout
        results.write()

    def run_impactmodel(self, solver=None, output=None, tol=1e-6, DisWeight=1.75, RatWeight=2):
        """Run the **MRIA** model with disruptions. This will return an economy with a new
        equilibrium, based on the new production and demand values.

        Parameters
        ----------
        solver
            Specify the solver to be used with Pyomo. The Default value is set to **None**. If set to **None**, the ipopt solver will be used
        output
            Specify whether you want the solver to print its progress.The default value is set to **None**
        tol
            the tolerance value that determines whether the outcome of the model is feasible. The default value is set to **1e-6**
        DisWeight
            the weight that determines the penalty set to let the model allow for additional imports. A higher penalty value will result in less imports. The default value is set to **1.75**
        RatWeight
            the weight that determines the penalty set to let the model allow to ration goods. A higher penalty value will result in less rationing. The default value is set to **2**

        Optionally write out the output of an optimized **MRIA_IO** class and **MRIA** model

        """
        model = self.m

        if solver is None:
            solver = 'ipopt'

        if DisWeight is None:
            DisWeight = 1.75

        if RatWeight is None:
            RatWeight = 2

        def demDisRat(model, R, S):
            return (
                self.Demand[R, S] == (sum(self.A_matrix[R, S, R, Sb]*self.X[R, Sb] for Sb in model.Sb) + self.lfd[R, S] + self.Rdem[R, S] - self.Rat[R, S]
                                      + sum(self.ImportShare[R, Rb, S]*(sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb] for Sb in model.Sb) +
                                                                        self.fd[Rb, S] + self.Rdem[Rb, S] - self.Rat[Rb, S]) for Rb in model.Rb if (R != Rb))
                                      + sum(self.ImportShare[R, Rb, S]*(self.DisImp[Rb, S])
                                            for Rb in model.Rb if (R != Rb))
                                      + self.ExpROW[R, S])
            )

        model.demDisRat = Constraint(model.R, model.S, rule=demDisRat, doc='Satisfy demand')

        def demsupDis(model, R, S):
            return (self.DisImp[R, S]+self.X[R, S]) >= self.Demand[R, S]

        model.demsupDis = Constraint(model.R, model.S, rule=demsupDis, doc='Satisfy demand')

        def DisImpA(model, R, S):
            return (self.DisImp[R, S]*(self.DisImp[R, S] + (self.Xbase[R, S]*self.X_up[R, S]) - self.Demand[R, S])) == 0

        model.DisImpA = Constraint(model.R, model.S, rule=DisImpA, doc='Satisfy demand')

        def DisImpEq(model, R, S):
            #    return m.DisImp[R, S] >=  (m.Demand[R, S] - (m.X[R, S]))
            return self.DisImp[R, S] >= (self.Demand[R, S] - (self.X[R, S]*self.X_up[R, S]))

#        model.DisImpEq = Constraint(model.R, model.S, rule=DisImpEq, doc='Satisfy demand')

        def ObjectiveDis2(model):
            return (
                sum(self.X[R, S] for S in model.S for R in model.R)
                + DisWeight*sum((self.Ratmarg[R, S]*self.DisImp[R, S])
                                for R in model.R for S in model.S)
                + RatWeight*sum((self.Ratmarg[R, S]*self.Rat[R, S])
                                for R in model.R for S in model.S)
                + sum((sum(self.ImportShare[R, Rb, S]*(sum(self.A_matrix[R, S, Rb, Sb]*self.X[Rb, Sb] for Sb in model.Sb) + self.fd[Rb, S] + self.Rdem[Rb, S] - self.Rat[Rb, S]) for Rb in model.Rb if (R != Rb))
                       + sum(self.ImportShare[R, Rb, S]*(self.DisImp[Rb, S]) for Rb in model.Rb if (R != Rb))) for R in model.R for S in model.S)
            )

        model.objective = Objective(rule=ObjectiveDis2, sense=minimize,
                                    doc='Define objective function')

        opt = SolverFactory(solver)
        if solver is 'ipopt':
            opt.options['max_iter'] = 5000
            opt.options['warm_start_init_point'] = 'yes'
            opt.options['warm_start_bound_push'] = 1e-6
            opt.options['warm_start_mult_bound_push'] = 1e-6
            opt.options['mu_init'] = 1e-6
            if tol != 1e-6:
                opt.options['tol'] = tol

        if output is None:
            opt.solve(model, tee=False)
        else:
            results = opt.solve(model, tee=True)
            # sends results to stdout
            results.write()
