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
    
    
    quset = qusetting()
    Operator = ['iswap1'  ,  'I' ]
    qtype = 1
    quset.qtype = qtype
    if quset.qtype == 1:
        a,sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,w01,w02 = initial(quset)
        psi0 = tensor(basis(quset.N,n) , basis(3,1) ,  basis(3,0))
        target = tensor(basis(quset.N,n) , basis(3,0) ,  basis(3,1))
        quset.iswap_deltat = deltat
        result , tlist = gate_evolution(psi0,Operator,setting = quset)
        
        rf01 = np.exp(1j*(w02[0])*tlist[-1])*basis(3,2)*basis(3,2).dag()
        rf02 = np.exp(1j*(w01[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf0 = basis(3,0)*basis(3,0).dag()+rf01+rf02
        rf11 = np.exp(1j*(w02[1])*tlist[-1])*basis(3,2)*basis(3,2).dag()
        rf12 = np.exp(1j*(w01[1])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf1 = basis(3,0)*basis(3,0).dag()+rf11+rf12
        U = tensor(qeye(quset.N),rf0,rf1)
    elif quset.qtype == 2:
        sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,w01,w02 = initial(quset)
        psi0 = tensor(basis(3,1) ,  basis(3,0))
        target = tensor(basis(3,0) ,  basis(3,1))
        quset.iswap_deltat = deltat
        result , tlist = gate_evolution(psi0,Operator,setting = quset)
        
        rf01 = np.exp(1j*(w02[0])*tlist[-1])*basis(3,2)*basis(3,2).dag()
        rf02 = np.exp(1j*(w01[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf0 = basis(3,0)*basis(3,0).dag()+rf01+rf02
        rf11 = np.exp(1j*(w02[1])*tlist[-1])*basis(3,2)*basis(3,2).dag()
        rf12 = np.exp(1j*(w01[1])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf1 = basis(3,0)*basis(3,0).dag()+rf11+rf12
        U = tensor(rf0,rf1)
    
    fid0=fidelity(ptrace(U*result.states[-1],1), ptrace(target,1))
    fid1=fidelity(ptrace(U*result.states[-1],2), ptrace(target,2))
    
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