# -*- coding: utf-8 -*-
"""
This script builds the MRIA model

@author: Elco Koks

@date: Nov, 2017

"""

from pyomo.environ import ConcreteModel,Set,SetOf,Param,Var,Constraint,Objective,minimize
import pandas as pd
import numpy as np
from pyomo.opt import SolverFactory
from vtra.analysis.mria.ratmarg import ratmarg_IO,ratmarg_SUT
import os
import sys


from vtra.utils import load_config

class MRIA_IO(object):

    """
    This is the class object 'MRIA' which is used to set up the modelling framework.

    We define the type of model, sets, set up the core variables and specify the
    constraints and objectives for different model setups.
    """

    def __init__(self, name, list_countries,list_sectors,EORA=False,list_fd_cats=[]):

        """
        Creation of a Concrete Model, specify the countries and sectors
        to include.
        """
        self.name = name
        self.m = ConcreteModel()
        self.countries = list_countries
        self.total_countries = len(list_countries)
        self.sectors = list_sectors
        self.fd_cat = list_fd_cats
        if EORA is True:
            self.EORA = True
        else:
            self.EORA = False

    def create_sets(self,FD_SET=[],VA_SET=[]):

        """
        Creation of the various sets. First step in future-proofing by allowing
        for own specification of set inputs
        """

        self.m.S = Set(initialize=self.sectors, doc='sectors')

        if self.EORA is True:
            self.m.rROW = Set(initialize=self.countries+['ROW'],ordered=True, doc='regions including export')
            self.m.R = Set(initialize=self.countries+['ROW'],ordered=True, doc='regions')
        else:
            self.m.rROW = Set(initialize=self.countries,ordered=True, doc='regions including export')
            self.m.R = Set(initialize=self.countries,ordered=True, doc='regions')

        if self.EORA is True:
            self.m.fdemand = Set(initialize=['P3h', 'P3n','P3g', 'P51','P52','P53'], doc='Final Demand')
        else:
            self.m.fdemand = Set(initialize=self.fd_cat, doc='Final Demand')

        if self.EORA is True:
            self.m.VA = Set(initialize=['VA'], doc='value added')
        else:
            self.m.VA = Set(initialize=VA_SET, doc='value added')

    def create_alias(self):
        """
        Set aliases
        """
        self.m.Rb   = SetOf(self.m.R)  # an alias of region R
        self.m.r   = SetOf(self.m.R)  # an alias of region R
        self.m.Sb   = SetOf(self.m.S)  # an alias of sector S


    """
    This part focuses on tables, parameters and variables
    """
    def create_A_mat(self,A_mat_in):
        model = self.m
        def A_matrix_init(model,R,S,Rb,Sb):
                    return A_mat_in[R,S,Rb,Sb]

        model.A_matrix = Param(model.R,model.S,model.R,model.Sb,initialize=A_matrix_init,doc = 'A matrix')

        self.A_matrix = model.A_matrix

    ''' Specify Final Demand and Local Final Demand'''
    def create_FD(self,FinalD,disr_dict_fd):

        disrupted_des = list(np.unique([x[1] for x in disr_dict_fd]))

        model = self.m

        model.Rdes = Set(initialize=disrupted_des, doc='Final Demand')

        def tfd_init(model,R,S,Rb):
            if (R,Rb,S) in list(disr_dict_fd.keys()):
                return sum(FinalD[R,S,Rb,fdemand] for fdemand in model.fdemand)*disr_dict_fd[R,Rb,S]
            else:
                return sum(FinalD[R,S,Rb,fdemand] for fdemand in model.fdemand)

        def fd_init(model,R,S):
                return sum(model.tfd[R,S,Rb] for Rb in model.Rb)


        model.tfd = Param(model.R, model.S,model.Rb, initialize=tfd_init, doc='Final Demand')

        model.fd = Param(model.R, model.S, initialize=fd_init, doc='Final Demand')

        self.ttfd = model.tfd
        self.fd = model.fd

    ''' Specify local final demand '''
    def create_LFD(self,FinalD):
        model = self.m
        def lfd_init(model,R,S):
            return sum(FinalD[R,S,R,fdemand] for fdemand in model.fdemand)
        model.lfd = Param(model.R, model.S, initialize=lfd_init, doc='Final Demand')

        self.lfd = model.lfd

    ''' Specify export and import to the rest of the world '''
    def create_ExpImp(self,ExpROW,ImpROW):
        model = self.m

        # Specify Export ROW
        def ExpROW_ini(m,R,S):
            return (ExpROW[R,S,'Export'])
        model.ExpROW = Param(model.R, model.S, initialize=ExpROW_ini, doc='Exports to the rest of the world')

        # Specify Import ROW
        def ImpROW_init(m,R,S):
            return (ImpROW[R,S,'Import'])
        model.ImpROW = Param(model.R, model.S, initialize=ImpROW_init, doc='Imports from the rest of the world')

        self.ExpROW = model.ExpROW
        self.ImpROW = model.ImpROW

    def create_ExpImp_EORA(self,Z_matrix):
        model = self.m

        # Specify Export ROW
        def ExpROW_ini(m,R,S):
            return (Z_matrix[R,S,'ROW','Total'])
        model.ExpROW = Param(model.R, model.S, initialize=ExpROW_ini, doc='Exports to the rest of the world')

        # Specify Import ROW
        def ImpROW_init(m,R,S):
            return (Z_matrix['ROW','Total',R,S])
        model.ImpROW = Param(model.R, model.S, initialize=ImpROW_init, doc='Imports from the rest of the world')

        self.ExpROW = model.ExpROW
        self.ImpROW = model.ImpROW


    """ Specify X variables """
    def create_X_up(self,disr_dict,Regmaxcap=0.98):
        model = self.m

        def shock_init(model, R,S):
            if (R,S) in list(disr_dict.keys()):
                return disr_dict[R,S]
            else:
                return 1.05

        model.X_up = Param(model.R,model.S,initialize=shock_init, doc='Maximum production capacity')
        self.X_up = model.X_up

    '''create Xbase'''
    def create_Xbase(self,Z_matrix,disr_dict,FinalD=None):
        model = self.m

        if self.fd.active is not True:
            self.create_FD(FinalD,disr_dict)

        if self.ExpROW.active is not True:
            self.create_ExpImp(Z_matrix)

        def x_init_base(model,R,S):
            return( sum(Z_matrix[R,S,Rb,Sb] for Rb in model.Rb for Sb in
                        model.Sb) + self.fd[R,S] + self.ExpROW[R,S])

        model.Xbase = Param(model.R, model.S,initialize=x_init_base,doc='Total Production baseline')
        self.Xbase = model.Xbase

    '''create X'''
    def create_X(self,disr_dict,Regmaxcap=0.98,
                 A_matrix_ini=None,Z_matrix=None,FinalD=None,Xbase=None,fd=None,ExpROW=None):

        model = self.m

        if self.Xbase.active is not True:
            self.create_Xbase(Z_matrix,FinalD)

        if self.A_matrix.active is not True:
            self.create_A_mat(A_matrix_ini)

        def X_bounds(model, R,S):
            if (R,S) in list(disr_dict.keys()):
                return (0.0, (1/Regmaxcap*self.Xbase[R,S])*disr_dict[R,S])
            else:
                return (0.0, (1/Regmaxcap*self.Xbase[R,S])*1.1)

        def x_init(model,R,S):
            return( sum(self.A_matrix[R,S,Rb,Sb]*self.Xbase[Rb,Sb] for Rb in model.Rb
                        for Sb in model.Sb) + self.fd[R,S] + self.ExpROW[R,S])

        model.X = Var(model.R, model.S, bounds=X_bounds,initialize=x_init, doc='Total Production')

        self.X = model.X

    """ Specify trade and value added """

    ''' Specify Value Added '''
    def create_VA(self,ValueA):
        model = self.m

        def va_init(model,R,S):
            return ValueA[R,S,'VA']

        model.ValueA = Param(model.R, model.S, initialize=va_init, doc='Value Added')

        self.ValueA = model.ValueA

    ''' Specify Trade between regions '''
    def create_Z_mat(self):
        model = self.m
        def Z_matrix_init(model,R,S,Rb,Sb):
            return self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb]

        model.Z_matrix = Param(model.R,model.S,model.R,model.Sb,initialize=Z_matrix_init,doc = 'Z matrix')
        self.Z_matrix = model.Z_matrix

    def create_Trade(self,FinalD,Z_matrix=None):
        model = self.m

        def Trade_init(model,R,Rb,S):
            while R != Rb:
                return sum(self.Z_matrix[Rb,S,R,i] for i in model.Sb)  + sum(FinalD[Rb,S,R,i] for i in model.fdemand)


        model.trade = Param(model.R,model.Rb, model.S, initialize=Trade_init, doc='Trade')
        self.trade = model.trade

    '''Estimate Total Export'''
    def create_TotExp(self):
        model = self.m
        def totexp_init(model, R, S):
            return sum(self.trade[Rb,R,S] for Rb in model.Rb if (R != Rb))

        model.TotExp = Param(model.R, model.S,initialize=totexp_init,doc='Total exports between regions')
        self.TotExp = model.TotExp

    '''Estimate Total Import'''
    def create_TotImp(self):
        model = self.m
        def totimp_init(model, R, S):
            return sum(self.trade[R,Rb,S] for Rb in model.Rb if (R != Rb))

        model.TotImp = Param(model.R, model.S,initialize=totimp_init,doc='Total imports between regions')
        self.TotImp = model.TotImp

    '''Estimate Import shares and Import share DisImp'''
    def create_ImpShares(self):
        model = self.m
        def impsh_init(model, R, Rb, S):
            while (self.trade[Rb,R,S] != None): # & ((sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S]) != None):
                try:
                    return self.trade[Rb,R,S]/(sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S])
                except ZeroDivisionError:
                    return 0

        def impshdis_init(model, R, Rb, S):
            return (sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S])

        model.ImportShare = Param(model.R, model.Rb, model.S,initialize=impsh_init,doc='Importshare of each region')
        model.ImportShareDisImp = Param(model.R, model.Rb,model.S,initialize=impshdis_init,doc='Importshare DisImp of each region')

        self.ImportShare = model.ImportShare
        self.ImportShareDisImp = model.ImportShareDisImp

    """ Specify specific variables for impact analysis """

    '''Reconstruction demand variable'''
    def create_Rdem(self):
        model = self.m
        model.Rdem = Param(model.R, model.S,initialize=0, doc='Reconstruction demand')
        self.Rdem = model.Rdem

    '''Rationing variable'''
    def create_Rat(self,FinalD=None,Z_matrix=None):
        model = self.m

        if self.lfd.active is not True:
            self.create_LFD(FinalD)

        if self.ExpROW.active is not True:
            self.create_ExpImp(Z_matrix)

        def Rat_bounds(model,R,S):
            return (0,abs(self.lfd[R,S]+self.ExpROW[R,S]))

        model.Rat = Var(model.R, model.S,bounds=Rat_bounds,initialize=0,doc='Rationing')
        self.Rat = model.Rat

    def create_Ratmarg(self,Table):
        model = self.m

        try:
            data_path = load_config()['paths']['data']
            RatMarg = pd.read_csv(os.path.join(data_path,'input_data','Ratmarg_{}.csv'.format(self.name)), index_col =[0],header=0)
            if self.EORA is True and (set(list(RatMarg.index.values)) != set(list(self.countries+['ROW']))):
                RatMarg = ratmarg_IO(Table,self.EORA)
            elif (set(list(RatMarg.index.values)) != set(list(self.countries))):
                RatMarg = ratmarg_IO(Table,self.EORA)
        except:
            RatMarg = ratmarg_IO(Table,self.EORA)

        Ratmarginal = {(r,k): v for r, kv in RatMarg.iterrows() for k,v in kv.to_dict().items()}

        model.Ratmarg = Param(model.R, model.S,initialize=Ratmarginal, doc='Rationing marginal',mutable=True)
        self.Ratmarg = model.Ratmarg

    '''Disaster import variable'''
    def create_DisImp(self,disr_dict,Regmaxcap=0.98):
        model = self.m

        disrupted_ctry = list(np.unique([x[0] for x in disr_dict]))

        #problem regions
        dimp_ctry = ['KEN','UGA']
        dimp_ind = ['i3']

        def Dis_bounds(model,R,S):
            if R in dimp_ctry and S in dimp_ind:
                return (0,0)
            elif (model.X_up[R,S] < (1.05) or R in disrupted_ctry):
                return (0,None)
            else:
                return (0,None)

        model.DisImp = Var(model.R, model.S,bounds=Dis_bounds,initialize=0, doc='Disaster Imports')
        self.DisImp = model.DisImp

    '''Specify demand function'''
    def create_demand(self):
        model = self.m

        def demand_init(model,R,S):
            return  (
            sum(self.A_matrix[R,S,R,Sb]*self.X[R,Sb] for Sb in model.Sb) + self.lfd[R,S] + self.Rdem[R,S] - self.Rat[R,S]
            +  sum(self.ImportShare[R,Rb,S]*(sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S] + self.Rdem[Rb,S]- self.Rat[R,S]) for Rb in model.Rb if (R != Rb))
            +  sum(self.ImportShare[R,Rb,S]*(self.DisImp[Rb,S]) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,S]
            )

        model.Demand = Var(model.R, model.S, bounds=(0.0,None),initialize=demand_init)
        self.Demand = model.Demand

    """ Create baseline dataset to use in model """
    def baseline_data(self,Table,disr_dict_sup,disr_dict_fd,EORA=None):

        if self.EORA is True:
            self.create_ExpImp_EORA(Table.Z_matrix)
        else:
            self.create_ExpImp(Table.ExpROW,Table.ImpROW)

        self.create_A_mat(Table.A_matrix)
        self.create_FD(Table.FinalD,disr_dict_fd)
        self.create_LFD(Table.FinalD)
        self.create_Xbase(Table.Z_matrix,Table.FinalD,disr_dict_fd)
        self.create_X(disr_dict_sup,Z_matrix=Table.Z_matrix,FinalD = Table.FinalD)
        self.create_VA(Table.ValueA)
        self.create_Z_mat()
        self.create_Trade(Table.FinalD)
        self.create_TotExp()
        self.create_TotImp()
        self.create_ImpShares()

    """ Create additional parameters and variables required for impact
    analysis """


    def impact_data(self,Table,disr_dict_sup,disr_dict_fd,Regmaxcap=0.98):


        self.create_X_up(disr_dict_sup)
        self.create_Rdem()
        self.create_Rat(Table.FinalD,Table.Z_matrix)
        self.create_Ratmarg(Table)
        self.create_DisImp(disr_dict_sup)
        self.create_demand()

    """
    Set up baseline model
    """
    def run_basemodel(self,solver=None):
        model = self.m

        if solver is None:
            solver = 'ipopt'

        def demSup(model, R,S):
            return  (self.X[R,S] >=
                     sum(self.A_matrix[R,S,R,Sb]*self.X[R,Sb] for Sb in model.Sb) + self.lfd[R,S]
            +  sum(self.ImportShare[R,Rb,S]*(sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S]) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,S])

        model.demSup = Constraint(model.R,model.S, rule=demSup, doc='Satisfy demand')

        def objective_base(model):
            return sum(self.X[R,S] for R in model.R for S in model.S)

        model.objective = Objective(rule=objective_base, sense=minimize, doc='Define objective function')

        opt = SolverFactory(solver)
        if solver is 'ipopt':
            opt.options['warm_start_init_point'] = 'yes'
            opt.options['warm_start_bound_push'] = 1e-6
            opt.options['warm_start_mult_bound_push'] = 1e-6
            opt.options['mu_init'] = 1e-6
        results = opt.solve(model,tee=True)
        #sends results to stdout
        results.write()

    def run_impactmodel(self,solver=None,output=None,tol=1e-6,DisWeight=1.75,RatWeight=2):
        model = self.m

        if solver is None:
            solver = 'ipopt'

        if DisWeight is None:
            DisWeight = 1.75

        if RatWeight is None:
            RatWeight = 2

        def demDisRat(model, R,S):
            return  (
            self.Demand[R,S] ==  (sum(self.A_matrix[R,S,R,Sb]*self.X[R,Sb] for Sb in model.Sb) + self.lfd[R,S] + self.Rdem[R,S] - self.Rat[R,S]
            +  sum(self.ImportShare[R,Rb,S]*(sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S] + self.Rdem[Rb,S]- self.Rat[Rb,S]) for Rb in model.Rb if (R != Rb))
            +  sum(self.ImportShare[R,Rb,S]*(self.DisImp[Rb,S]) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,S]            )
            )

        model.demDisRat = Constraint(model.R,model.S, rule=demDisRat, doc='Satisfy demand')

        def demsupDis(model,R,S):
             return (self.DisImp[R,S]+self.X[R,S]) >= self.Demand[R,S]

        model.demsupDis = Constraint(model.R,model.S, rule=demsupDis, doc='Satisfy demand')

        def DisImpA(model,R,S):
            return (self.DisImp[R,S]*(self.DisImp[R,S] + (self.Xbase[R,S]*self.X_up[R,S]) - self.Demand[R,S])) == 0

        model.DisImpA = Constraint(model.R,model.S, rule=DisImpA, doc='Satisfy demand')

        def DisImpEq(model,R,S):
        #    return m.DisImp[R,S] >=  (m.Demand[R,S] - (m.X[R,S]))
            return self.DisImp[R,S] >=  (self.Demand[R,S] - (self.X[R,S]*self.X_up[R,S]))

#        model.DisImpEq = Constraint(model.R,model.S, rule=DisImpEq, doc='Satisfy demand')


        def ObjectiveDis2(model):
            return (
                sum(self.X[R,S] for S in model.S for R in model.R)
                + DisWeight*sum((self.Ratmarg[R,S]*self.DisImp[R,S]) for R in model.R for S in model.S)
                + RatWeight*sum((self.Ratmarg[R,S]*self.Rat[R,S]) for R in model.R for S in model.S)
                + sum((sum(self.ImportShare[R,Rb,S]*(sum(self.A_matrix[R,S,Rb,Sb]*self.X[Rb,Sb] for Sb in model.Sb) + self.fd[Rb,S] + self.Rdem[Rb,S] - self.Rat[Rb,S]) for Rb in model.Rb if (R != Rb))
                +  sum(self.ImportShare[R,Rb,S]*(self.DisImp[Rb,S]) for Rb in model.Rb if (R != Rb))) for R in model.R for S in model.S)
                )

        model.objective = Objective(rule=ObjectiveDis2, sense=minimize, doc='Define objective function')

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
            opt.solve(model,tee=False)
        else:
            results = opt.solve(model,tee=True)
            #sends results to stdout
            results.write()


class MRIA_SUT(object):

    """
    This is the class object 'MRIA' which is used to set up the modelling framework.

    We define the type of model, sets, set up the core variables and specify the
    constraints and objectives for different model setups.
    """

    def __init__(self, name, list_countries,list_sectors,list_products,EORA=False):

        """
        Creation of a Concrete Model, specify the countries and sectors
        to include.
        """
        self.name = name
        self.m = ConcreteModel()
        self.countries = list_countries
        self.total_countries = len(list_countries)
        self.sectors = list_sectors
        self.products = list_products

        if EORA is True:
            self.EORA = True
        else:
            self.EORA = False

    def create_sets(self,FD_SET=['FinalD'],VA_SET=['VA']):

        """
        Creation of the various sets. First step in future-proofing by allowing
        for own specification of set inputs
        """

        self.m.S = Set(initialize=self.sectors, doc='sectors')
        self.m.P = Set(initialize=self.products, doc='sectors')
        self.m.row = Set(initialize=self.products, doc='products')
        self.m.col = Set(initialize=self.sectors+['FinalD'], doc='sectors and final demand')

        self.m.rROW = Set(initialize=self.countries,ordered=True, doc='regions including export')
        self.m.R = Set(initialize=self.countries,ordered=True, doc='regions')

        self.m.fdemand = Set(initialize=FD_SET, doc='Final Demand')

        self.m.VA = Set(initialize=VA_SET, doc='value added')

    def create_alias(self):
        """
        Set aliases
        """
        self.m.Rb   = SetOf(self.m.R)  # an alias of region R
        self.m.r   = SetOf(self.m.R)  # an alias of region R
        self.m.Sb   = SetOf(self.m.S)  # an alias of sector S


    """
    This part focuses on tables, parameters and variables
    """
    def create_UseAbs(self,REG_USE):

        model = self.m

        def usetable_init(model,R,P,Rb,col):
            return REG_USE[R,P,Rb,col]

        model.UseAbs = Param(model.R,model.P,model.Rb,model.col,initialize=usetable_init,doc='Absolute use table')

        self.UseAbs = model.UseAbs

    def create_SupAbs(self,REG_SUP):
        model = self.m

        def suptable_init(model,R,S,Rb,P):
            return REG_SUP[R,S,Rb,P]

        model.SupAbs = Param(model.R,model.S,model.Rb,model.P,initialize=suptable_init,doc='Absolute sup table')

        self.SupAbs = model.SupAbs

    def create_TotSup(self):
        model = self.m
        def totsup_init(model,R,S):
            return sum(self.SupAbs[Rb,S,R,P] for Rb in model.Rb for P in model.P)

        model.TotSup = Param(model.R,model.S,initialize=totsup_init)
        self.TotSup = model.TotSup

    def create_Sup(self):
        model = self.m
        def sup_init(model,R,S,P):
            return sum(self.SupAbs[R,S,Rb,P] for Rb in model.Rb)/self.TotSup[R,S]

        model.Sup = Param(model.R,model.S,model.P,initialize=sup_init)
        self.Sup = model.Sup

    '''create Xbase'''
    def create_Xbase(self):
        model = self.m

        def xbase_init(model,R,S):
            return sum(self.SupAbs[Rb,S,R,P] for Rb in model.Rb for P in model.P)

        model.Xbase = Param(model.R,model.S,initialize=xbase_init)
        self.Xbase = model.Xbase

    def create_Use(self):
        model = self.m
        def use_init(model,R,P,S):
            return sum(self.UseAbs[Rb,P,R,S] for Rb in model.Rb)/self.Xbase[R,S]

        model.Use = Param(model.R,model.P,model.S,initialize=use_init)
        self.Use = model.Use


    """ Specify X variables """
    def create_X_up(self,disruption=1.1,disrupted_ctry=[],disrupted_sctr=[],Regmaxcap=0.98):
        model = self.m

        def shock_init(model, R,S):
            if R in disrupted_ctry and S in disrupted_sctr:
                return 1/Regmaxcap*disruption
            else:
                return 1/Regmaxcap*1.1

        model.X_up = Param(model.R,model.S,initialize=shock_init, doc='Maximum production capacity')
        self.X_up = model.X_up


    '''create X'''
    def create_X(self,disruption=1.1,disrupted_ctry=[],disrupted_sctr=[],Regmaxcap=0.98):

        model = self.m

        def X_bounds(model, R,S):
            if R in disrupted_ctry and S in disrupted_sctr:
                return (0.0, (1/Regmaxcap*self.Xbase[R,S])*disruption)
            else:
                return (0.0, (1/Regmaxcap*self.Xbase[R,S])*1.1)

        def x_init(model,R,S):
            return sum(self.SupAbs[Rb,S,R,P] for Rb in model.Rb for P in model.P)

        model.X = Var(model.R, model.S, bounds=X_bounds,initialize=x_init, doc='Total Production')

        self.X = model.X

    def create_fd(self,REG_USE):

        model = self.m

        def findem_init(model,R,P):
            return sum(REG_USE[Rb,P,R,fdemand] for Rb in model.Rb for fdemand in model.fdemand)

        model.fd = Param(model.R,model.P,initialize=findem_init)

        self.fd = model.fd

    ''' Specify Trade between regions '''
    def create_Trade(self,REG_USE):
        model = self.m

        def Trade_init(model,R,Rb,P):
            return sum(REG_USE[R,P,Rb,col] for col in model.col if (R != Rb))

        model.trade = Param(model.R,model.Rb,model.P,initialize=Trade_init)

        self.trade = model.trade

    '''Estimate Total Export'''
    def create_TotExp(self):
        model = self.m
        def totexp_init(model, R, P):
            return sum(self.trade[R,Rb,P] for Rb in model.Rb if (R != Rb))

        model.TotExp = Param(model.R, model.P,initialize=totexp_init,doc='Total exports between regions')
        self.TotExp = model.TotExp

    '''Estimate Export ROW'''
    def create_ExpImp(self,ExpROW_in):

        model = self.m
        # Specify Export ROW
        def ExpROW_ini(m,R,P):
            return (ExpROW_in[R,P,'Exports'])

        def ImpROW_ini(m,R,P):
            return (ExpROW_in[R,P,'Exports'])*0

        model.ExpROW = Param(model.R, model.P, initialize=ExpROW_ini, doc='Exports to the rest of the world')
        model.ImpROW = Param(model.R, model.P, initialize=ImpROW_ini, doc='Imports from the rest of the world')

        self.ExpROW = model.ExpROW
        self.ImpROW = model.ImpROW


    '''Estimate Total Import'''
    def create_TotImp(self):
        model = self.m
        def totimp_init(model, R, P):
            return sum(self.trade[Rb,R,P] for Rb in model.Rb if (R != Rb))

        model.TotImp = Param(model.R, model.P,initialize=totimp_init,doc='Total imports between regions')
        self.TotImp = model.TotImp

    '''Estimate Import shares and Import share DisImp'''
    def create_ImpShares(self):
        model = self.m
        def impsh_init(model, R, Rb, P):
            while self.trade[Rb,R,P] != None:
                return self.trade[Rb,R,P]/self.TotImp[R,P]

        model.ImportShare = Param(model.R, model.Rb, model.P,initialize=impsh_init,doc='Importshare of each region')
        model.ImportShareDisImp = Param(model.R, model.Rb,model.P,initialize=impsh_init,doc='Importshare DisImp of each region')

        self.ImportShare = model.ImportShare
        self.ImportShareDisImp = model.ImportShareDisImp

        '''Estimate Importratio'''
    def create_ImportRatio(self):
        model = self.m

        def imprat_init(model, R, P):
                return (self.TotImp[R,P] + self.ImpROW[R,P])/(sum(self.Xbase[R,S]*self.Use[R,P,S] for S in model.S) + self.fd[R,P])

        model.Importratio = Param(model.R, model.P,initialize=imprat_init,doc='Importratio of each region')

        self.Importratio = model.Importratio


    """ Specify specific variables for impact analysis """

    '''Reconstruction demand variable'''
    def create_Rdem(self):
        model = self.m
        model.Rdem = Param(model.R, model.P,initialize=0, doc='Reconstruction demand')
        self.Rdem = model.Rdem

    '''Rationing variable'''
    def create_Rat(self):
        model = self.m

        def Rat_bounds(model,R,P):
            return (0,abs(self.fd[R,P]+self.ExpROW[R,P]))

        model.Rat = Var(model.R, model.P,bounds=Rat_bounds,initialize=0,doc='Rationing')
        self.Rat = model.Rat

    def create_Ratmarg(self,Table):
        model = self.m

        try:
            RatMarg = pd.read_csv('..\..\input_data\Ratmarg_Vale_SuT.csv', index_col =[0],header=0)
            if (set(list(RatMarg.index.values)) != set(list(self.countries))):
                RatMarg = ratmarg_SUT(Table)
        except:
            RatMarg = ratmarg_SUT(Table)

#        RatMarg = pd.read_csv('..\..\input_data\Ratmarg_Vale_SuT.csv', index_col =[0],header=0)
        Ratmarginal = {(r,k): v for r, kv in RatMarg.iterrows() for k,v in kv.to_dict().items()}

        model.Ratmarg = Param(model.R, model.P,initialize=Ratmarginal, doc='Rationing marginal',mutable=True)
        self.Ratmarg = model.Ratmarg

    '''Disaster import variable'''
    def create_DisImp(self,disrupted_ctry={},Regmaxcap=0.98):
        model = self.m

        model.DisImp = Var(model.R, model.P,bounds=(0,None),initialize=0, doc='Disaster Imports')
        self.DisImp = model.DisImp

    '''Specify demand function'''
    def create_demand(self):
        model = self.m

        def demand_init(model,R,P):
            return  (
            (sum(self.X[R,S]*self.Use[R,P,S] for S in model.S) + self.fd[R,P] - self.Rat[R,P])*(1-self.Importratio[R,P])
            +  sum(self.ImportShare[Rb,R,P]*((sum(self.X[Rb,S]*self.Use[Rb,P,S] for S in model.S) + self.fd[Rb,P]- self.Rat[R,P])*(self.Importratio[Rb,P])) for Rb in model.Rb if (R != Rb))
            +  sum(self.ImportShare[R,Rb,P]*(self.DisImp[Rb,P]) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,P]
            )


        model.Demand = Var(model.R, model.P, bounds=(0.0,None),initialize=demand_init)
        self.Demand = model.Demand

    """
    Set up baseline model
    """

    """ Create baseline dataset to use in model """
    def baseline_data(self,Table,disruption=1.1,disrupted_ctry=[],disrupted_sctr=[],EORA=None):

        self.create_UseAbs(Table.Use)
        self.create_SupAbs(Table.Sup)
        self.create_TotSup()
        self.create_Sup()
        self.create_Xbase()
        self.create_Use()
        self.create_X(disruption,disrupted_ctry,disrupted_sctr)
        self.create_fd(Table.Use)
        self.create_Trade(Table.Use)
        self.create_TotExp()
        self.create_TotImp()
        self.create_ExpImp(Table.ExpROW)
        self.create_ImpShares()
        self.create_ImportRatio()

    """ Create additional parameters and variables required for impact
    analysis """
    def impact_data(self,Table,disruption=1.1,disrupted_ctry=[],disrupted_sctr=[],Regmaxcap=0.98):
        self.create_X_up(disruption,disrupted_ctry,disrupted_sctr)
        self.create_Rdem()
        self.create_Rat()
        self.create_Ratmarg(Table)
        self.create_DisImp(disrupted_ctry)
        self.create_demand()

    def run_basemodel(self,solver=None):
        model = self.m

        def test_init(m,R,P):
            return ((sum(self.X[R,S]*self.Use[R,P,S] for S in model.S) + self.fd[R,P])*(1-self.Importratio[R,P])
            +  sum(self.ImportShare[Rb,R,P]*((sum(self.X[Rb,S]*self.Use[Rb,P,S] for S in model.S) + self.fd[Rb,P])*(self.Importratio[Rb,P])) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,P])

        model.test = Param(model.R,model.P,initialize=test_init)

        self.test = model.test

        if solver is None:
            solver = 'ipopt'

        #Demand equals supply, baseline equation
        def demSup(m, R,P,S):
            return  (
            sum(self.X[R,S]*self.Sup[R,S,P] for S in model.S) >= (sum(self.X[R,S]*self.Use[R,P,S] for S in model.S) + self.fd[R,P])*(1-self.Importratio[R,P])
            +  sum(self.ImportShare[Rb,R,P]*((sum(self.X[Rb,S]*self.Use[Rb,P,S] for S in model.S) + self.fd[Rb,P])*(self.Importratio[Rb,P])) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,P]
            )

        model.demSup = Constraint(model.R,model.P,model.S, rule=demSup, doc='Satisfy demand')

        def objective_base(model):
            return sum(self.X[R,S] for R in model.R for S in model.S)

        model.objective = Objective(rule=objective_base, sense=minimize, doc='Define objective function')

        opt = SolverFactory(solver)
        if solver is 'ipopt':
            opt.options['warm_start_init_point'] = 'yes'
            opt.options['warm_start_bound_push'] = 1e-6
            opt.options['warm_start_mult_bound_push'] = 1e-6
            opt.options['mu_init'] = 1e-6
        results = opt.solve(model,tee=True)
        #sends results to stdout
        results.write()

    def run_impactmodel(self,solver=None,tol=None,output=None,DisWeight=None,RatWeight=None,Regmaxcap=0.98):
        model = self.m

        if solver is None:
            solver = 'ipopt'

        if DisWeight is None:
            DisWeight = 0.75

        if RatWeight is None:
            RatWeight = 2

        def demDisRat(model, R,P):
            return  (
            self.Demand[R,P] == ((sum(self.X[R,S]*self.Use[R,P,S] for S in model.S) + self.fd[R,P] - self.Rat[R,P])*(1-self.Importratio[R,P])
            +  sum(self.ImportShare[Rb,R,P]*((sum(self.X[Rb,S]*self.Use[Rb,P,S] for S in model.S) + self.fd[Rb,P]- self.Rat[R,P])*(self.Importratio[Rb,P])) for Rb in model.Rb if (R != Rb))
            +  sum(self.ImportShare[R,Rb,P]*(self.DisImp[Rb,P]) for Rb in model.Rb if (R != Rb))
            + self.ExpROW[R,P])
            )

        model.demDisRat = Constraint(model.R,model.P, rule=demDisRat, doc='Satisfy demand')

        def demsupDis(model,R,P):
             return (self.DisImp[R,P]+ sum(self.X[R,S]*self.Sup[R,S,P] for S in model.S)) >= self.Demand[R,P]

        model.demsupDis = Constraint(model.R,model.P, rule=demsupDis, doc='Satisfy demand')

#        def DisImpA(model,R,P):
#             return self.DisImp[R,P]*(self.DisImp[R,P] - Regmaxcap*(sum((self.X[R,S]*self.X_up[R,S])*self.Sup[R,S,P] for S in model.S) - self.Demand[R,P])) == 0
##             return self.DisImp[R,P]*(self.DisImp[R,P] - (sum((self.X[R,S])*self.Sup[R,S,P] for S in model.S) - self.Demand[R,P])) == 0
#
#        model.DisImpA = Constraint(model.R,model.P,rule=DisImpA, doc='Satisfy demand')

        def DisImpEq(model,R,P):
             return (self.DisImp[R,P] >=  (self.Demand[R,P] - Regmaxcap*sum((self.X[R,S]*self.X_up[R,S])*self.Sup[R,S,P] for S in model.S)))

        model.DisImpEq = Constraint(model.R,model.P,rule=DisImpEq, doc='Satisfy demand')

        def ObjectiveDis2(model):
            return (
                sum(sum(self.X[R,S]*self.Sup[R,S,P] for S in model.S) for R in model.R for P in model.P)
                + DisWeight*sum((self.Ratmarg[R,P]*self.DisImp[R,P]) for R in model.R for P in model.P)
                + RatWeight*sum((self.Ratmarg[R,P]*self.Rat[R,P]) for R in model.R for P in model.P)
                + sum((sum(self.ImportShare[Rb,R,P]*((sum(self.X[Rb,S]*self.Use[Rb,P,S] for S in model.S) + self.fd[Rb,P]+ self.Rdem[Rb,P] - self.Rat[Rb,P])*(self.Importratio[Rb,P])) for Rb in model.Rb if (R != Rb))
                +  sum(self.ImportShareDisImp[Rb,R,P]*(self.TotImp[Rb,P]/(self.TotImp[Rb,P]+self.ImpROW[Rb,P])*self.DisImp[Rb,P]) for Rb in model.Rb if (Rb != R))) for R in model.R for P in model.P)
                )

        model.objective = Objective(rule=ObjectiveDis2, sense=minimize, doc='Define objective function')

        opt = SolverFactory(solver)
        if solver is 'ipopt':
            opt.options['max_iter'] = 5000
            opt.options['warm_start_init_point'] = 'yes'
            opt.options['warm_start_bound_push'] = 1e-6
            opt.options['warm_start_mult_bound_push'] = 1e-6
            opt.options['mu_init'] = 1e-6
            if tol is None:
                opt.options['tol'] = 1e-6
            else:
                opt.options['tol'] = tol

        if output is None:
            opt.solve(model,tee=False)
        else:
            results = opt.solve(model,tee=True)
            #sends results to stdout
            results.write()
