# -*- coding: utf-8 -*-
"""
PURPOSE: estimate a new matrix X with exogenously given row and column
totals that is a close as possible to a given original matrix X0 using
the Generalized RAS (GRAS) approach
-------------------------------------------------------------------------
USAGE: X = gras(X0,u,v) OR [X,r,s] = gras(X0,u,v) with or without eps
included as the fourth argument, where
INPUT:
-> X0 = benchmark (base) matrix, not necessarily square
-> u = column vector of (new) row totals
-> v = column vector of (new) column totals
-> eps = convergence tolerance level; if empty, the default threshold
is 0.1e-5 (=0.000001)
OUTPUT:
-> X = estimated/adjusted/updated matrix
-> r = substitution effects (row multipliers)
-> s = fabrication effects (column multipliers)
-------------------------------------------------------------------------
REFERENCES: 1) Junius T. and J. Oosterhaven (2003), The solution of
updating or regionalizing a matrix with both positive and negative
entries, Economic Systems Research, 15, pp. 87-96.
2) Lenzen M., R. Wood and B. Gallego (2007), Some comments on the GRAS
method, Economic Systems Research, 19, pp. 461-465.
3) Temurshoev, U., R.E. Miller and M.C. Bouwmeester (2013), A note on the
GRAS method, Economic Systems Research, 25, pp. 361-367.
-------------------------------------------------------------------------
NOTE FROM THE AUTHOR: If you use this program and publish the results in
the form of working/discussion papers, journal articles etc., you are
kindly asked to cite the third paper mentioned above (as this code is the
online Appendix to that paper).
-------------------------------------------------------------------------
Written by:   Umed Temurshoev, 07/10/2010 with later adjustments
              Current e-mail: umed.temurshoev@ec.europa.eu
"""

import numpy as np
import pandas as pd

from prepare_table import load_config,load_sectors,load_table,aggregate_table,load_provincial_stats,estimate_gva
from read_table import load_output

def invd(x):
    invd = 1./x
    invd[x==0] = 1
    return np.diag(invd)   

def ras_method(X0, u, v, eps=1e-5):
    
    m, n = np.shape(X0)
    N = np.zeros((m,n))
    N[X0<0] = -X0[X0<0]
    P = X0+N
    
    #initial guess for r (suggested by J&O, 2003)
    r = np.ones((m))      
    pr = np.dot(P.T, r)
    nr = N.T.dot(invd(r)).dot(np.ones((m)))
    s1 = np.dot(invd(2*pr),(v+np.sqrt((np.square(v)+4*(pr.dot(nr))))))
    ss = -invd(v).dot(nr)
    s1[pr==0] = ss[pr==0]      

    ps = np.dot(P, s1)
    ns = N.dot(invd(s1)).dot(np.ones((n)))
    r = np.dot(invd(2*ps),(u+np.sqrt((np.square(u)+4*(ps.dot(ns))))))
    rr = - invd(u).dot(ns)
    r[ps==0] = rr[ps==0]
 
    pr = np.dot(P.T,r)
    nr = N.T.dot(invd(r)).dot(np.ones((m)))
    
    #%second step s
    s2 = np.dot(invd(2*pr),v+np.sqrt((np.square(v)+4*(pr.dot(nr)))))      
    ss = -invd(v).dot(nr)
    s2[pr==0] = ss[pr==0]      
    dif = s2-s1

    M = np.max(abs(dif))
    i = 1                #first iteration
    while (M > eps):
        print(M)
        s1 = s2
        ps = P.dot(s1)
        ns = N.dot(invd(s1)).dot(np.ones((n)))
        r = np.dot(invd(2*ps),(u+np.sqrt((np.square(u)+4*(ps.dot(ns))))))
        rr = -invd(u).dot(ns)
        r[ps==0] = rr[ps==0]
        pr = P.T.dot(r)
        nr = N.T.dot(invd(r)).dot(np.ones((m)))
        s2 = np.dot(invd(2*pr),v+np.sqrt((np.square(v)+4*(pr.dot(nr)))))    
        ss = -invd(v).dot(nr)
        s2[pr==0] = ss[pr==0]      
        dif = s2-s1
        i = i+1
        M = np.max(abs(dif))
        if i == 10000:
            break

    #%final step s        
    s = s2                                        
    ps = P.dot(s)
    ns = N.dot(invd(s)).dot(np.ones((n)))
    r = np.dot(invd(2*ps),(u+np.sqrt((np.square(u)+4*(ps.dot(ns))))))
    rr = -invd(u).dot(ns)
    r[ps==0] = rr[ps==0]
    return np.diag(r).dot(P).dot(np.diag(s))-invd(r).dot(N).dot(invd(s))      #%updated matrix

