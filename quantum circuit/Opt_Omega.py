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
    '''
    对X,Y门进行参数标定
    输入X,Y的可变参量omega
    返回保真度和对应的参量
    '''
    quset = qusetting()
    Operator = ['X'  ,  'X' ]
    qtype = 1
    quset.qtype = qtype
    if quset.qtype == 1:
        a,sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,En = initial(quset)
        psi0 = tensor(basis(quset.N,n) , basis(3,1) ,  (basis(3,0)+1*basis(3,1)).unit())
        #初态
        target = sxm[1]*sxm[0] * psi0   
        #目标状态
        quset.omega = omega
        result , tlist = gate_evolution(psi0,Operator,setting = quset)
        

        rf0 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(En[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(En[1])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        U = tensor(qeye(quset.N),rf0,rf1)

        fid=fidelity(U*result.states[-1]*result.states[-1].dag()*U.dag(), target)
    elif quset.qtype == 2:
        sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,En = initial(quset)
        psi0 = tensor(basis(3,1) ,  (basis(3,0)+1*basis(3,1)).unit())
        target = sxm[1]*sxm[0] * psi0   
        quset.omega = omega
        result , tlist = gate_evolution(psi0,Operator,setting = quset)

        rf0 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(En[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(En[1])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        U = tensor(rf0,rf1)

        fid=fidelity(U*result.states[-1]*result.states[-1].dag()*U.dag(), target)
    
    return([fid,omega])

if __name__ == '__main__':
    starttime=clock()
    g = linspace(0.0333,0.0334,5)
    p = Pool(4)
    A = p.map(Opt_Omega,g)
    p.close()
    p.join()
    fid =  np.array([x[0] for x in A])#取输出中的第一列
    omega = np.array([x[1] for x in A])#取输出中的第二列
    opt = omega[np.where(fid== max(fid))]
    
    print(opt[0]  , max(fid))

    
    finishtime=clock()
    print( 'Time used: ', (finishtime-starttime), 's')