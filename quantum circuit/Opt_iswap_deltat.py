#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 22:20:06 2017

@author: chen
"""

from time import clock
from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from gate_evolution import *
from multiprocessing import Pool
import os

def Opt_iswap_deltat(deltat):
    gateset = gatesetting()
    Operator = ['iswap1'  ,  'I' ]
    psi0 = tensor(basis(N,n) , basis(level,1) ,  basis(level,0))
    target = tensor(basis(N,n) , basis(level,0) ,  basis(level,1))
    
    gateset.iswap_deltat = deltat
    result , tlist = gate_evolution(psi0 , Operator,gateset = gateset)
    fid0=fidelity(ptrace(result.states[-1],1), ptrace(target,1))
    fid1=fidelity(ptrace(result.states[-1],2), ptrace(target,2))
    
    return([fid0,fid1,deltat])

if __name__ == '__main__':
    starttime=clock()
    g = linspace(20,50,5)
#    g = [0.033358578617]
    p = Pool(4)
    A = p.map(Opt_iswap_deltat,g)
    p.close()
    p.join()
    fid0 =  np.array([x[0] for x in A])
    fid1 =  np.array([x[1] for x in A])
    iswap_deltat = np.array([x[2] for x in A])
    opt0 = iswap_deltat[np.where(fid0== max(fid0))]
    opt1 = iswap_deltat[np.where(fid1== max(fid1))]
    print(opt0[0]  , max(fid0))
    print(opt1[0] , max(fid1))

    
    finishtime=clock()
    print( 'Time used: ', (finishtime-starttime), 's')