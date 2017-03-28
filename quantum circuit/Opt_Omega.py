#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 21:54:08 2017

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

#def Opt_Omega(Omega_list):
#    gateset = gatesetting()
#    Operator = ['X'  ,  'X' ]
#    fid = []
#    target = sxm0 * psi0
#    for omega in Omega_list:
#        gateset.omega = omega
#        result = gate_evolution(Operator,gateset = gateset)
#        fid.append(fidelity(result.states[-1],target))
#        
#    opt = Omega_list[np.where(max(fid))]
#    return(opt,max(fid))
#
#if __name__ == '__main__':
#    starttime=clock()
#    g = linspace(0.02,0.04,10)
#
#    opt,Max = Opt_Omega(g)
#    print(opt[0])
#    print(Max)
#    
#    finishtime=clock()
#    print( 'Time used: ', (finishtime-starttime), 's')


def Opt_Omega(omega):
    gateset = gatesetting()
    Operator = ['X'  ,  'X' ]
    psi0 = tensor(basis(N,n) , basis(level,1) ,  (basis(level,0)+1*basis(level,1)).unit())
    target = sxm1*sxm0 * psi0   
    gateset.omega = omega
    result , tlist = gate_evolution(psi0,Operator,gateset = gateset)
    fid0=fidelity(ptrace(result.states[-1],1), ptrace(target,1))
    fid1=fidelity(ptrace(result.states[-1],2), ptrace(target,2))
    
    return([fid0,fid1,omega])

if __name__ == '__main__':
    starttime=clock()
    g = linspace(0.0333,0.0334,10)
#    g = [0.033358578617]
    p = Pool(4)
    A = p.map(Opt_Omega,g)
    p.close()
    p.join()
    fid0 =  np.array([x[0] for x in A])
    fid1 =  np.array([x[1] for x in A])
    omega = np.array([x[2] for x in A])
    opt0 = omega[np.where(fid0== max(fid0))]
    opt1 = omega[np.where(fid1== max(fid1))]
    print(opt0[0]  , max(fid0))
    print(opt1[0] , max(fid1))

    
    finishtime=clock()
    print( 'Time used: ', (finishtime-starttime), 's')