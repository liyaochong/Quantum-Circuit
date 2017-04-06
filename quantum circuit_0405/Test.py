#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 21:36:26 2017

@author: chen
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from gate_evolution import *


from time import clock
starttime=clock()

quset = qusetting()
Operator = ['CZ1','I']
#Operator = ['X'  ,  'X' ]
#Operator = ['iswap1'  ,  'I' ]
#psi0 = tensor(basis(N,n) , basis(3,1) ,  (basis(3,0)+basis(3,1)).unit())
#psi0 = tensor(basis(N,n) , (basis(level,0)+basis(level,1)).unit() ,  (basis(level,0)+basis(level,1)).unit())
#psi0 = tensor(basis(N,n) , basis(level,1) ,  basis(level,0))
psi0 = tensor(basis(3,1) ,  (basis(3,0)+basis(3,1)).unit())
quset.DRAG = True
#quset.omega = 0.033243043043
quset.CZ_deltat = 400
#quset.iswap_deltat = 20
quset.qtype = 2 

result , tlist = gate_evolution(psi0 , Operator , setting = quset)
evolutionplot(0 , result , tlist , setting = quset)
evolutionplot(1 , result , tlist , setting = quset)


#rf01 =np.exp(1j*(w02[0])*tlist[-1])*basis(3,2)*basis(3,2).dag()
#rf02 = np.exp(1j*(w01[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
#rf0 = basis(3,0)*basis(3,0).dag()+rf01+rf02
#rf11 = np.exp(1j*(w02[1])*tlist[-1])*basis(3,2)*basis(3,2).dag()
#rf12 = np.exp(1j*(w01[1])*tlist[-1])*basis(3,1)*basis(3,1).dag()
#rf1 = basis(3,0)*basis(3,0).dag()+rf11+rf12
#U = tensor(qeye(N),rf0,rf1)
#
##target = sxm0*sxm1*psi0
#target = tensor(basis(N,n) , basis(3,1) ,  (basis(level,0)-basis(level,1)).unit())
#fid0=fidelity(ptrace(U*result.states[-1],1), ptrace(target,1))
#fid1=fidelity(ptrace(U*result.states[-1],2), ptrace(target,2))
#print(fid0,fid1)


finishtime=clock()
print ('Time used: ', (finishtime-starttime), 's')