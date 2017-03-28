#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 20:44:16 2017

@author: chen
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *




#==============================================================================    
'''
qubit's parameters
'''
w_c = 7.0  * 2 * np.pi  # cavity frequency
w_q = np.array([ 5.0 , 6.0 ]) * 2 * np.pi
g = np.array([0.2 , 0.2]) * 2 * np.pi
eta_q = np.array([0.8 , 0.8]) * 2 * np.pi
N = 3              # number of cavity fock states
n= 0
level = 3  #能级数
#==============================================================================
#psi0 = tensor(basis(N,n) , basis(level,0) ,  (basis(level,0)+1j*basis(level,1)).unit())
#tlist = np.linspace(0,41,201)
#==============================================================================

#==============================================================================
'''
Operators
'''
a = tensor(destroy(N),qeye(level),qeye(level))
sm0 = tensor(qeye(N),destroy(level),qeye(level))
sm1 = tensor(qeye(N),qeye(level),destroy(level))

E_uc0 = tensor(qeye(N),basis(3,2)*basis(3,2).dag(),qeye(level))
E_uc1 = tensor(qeye(N),qeye(level), basis(3,2)*basis(3,2).dag())#用以表征非简谐性的对角线最后一项(非计算能级)
#E_uc1 = tensor(qeye(N),qeye(level), Qobj([[0,0],[0,1]]))

E_e0 = tensor(qeye(N),basis(3,1)*basis(3,1).dag(),qeye(level))
E_e1 = tensor(qeye(N),qeye(level),basis(3,1)*basis(3,1).dag())  #激发态

E_g0 = tensor(qeye(N),basis(3,0)*basis(3,0).dag(),qeye(level))
E_g1 = tensor(qeye(N),qeye(level),basis(3,0)*basis(3,0).dag())  #基态

sn0 = sm0.dag()*sm0
sn1 = sm1.dag()*sm1

sx0 = sm0.dag()+sm0;sxm0 = tensor(qeye(N),Qobj([[0,1,0],[1,0,0],[0,0,0]]),qeye(level))
sx1 = sm1.dag()+sm1;sxm1 = tensor(qeye(N),qeye(level),Qobj([[0,1,0],[1,0,0],[0,0,0]]))

sy0 = -1j*(sm0.dag()-sm0);sym0 = tensor(qeye(N),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]),qeye(level))
sy1 = -1j*(sm1.dag()-sm1);sym1 = tensor(qeye(N),qeye(level),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]))

sz0 = E_g0 - E_e0
sz1 = E_g1 - E_e1

#==============================================================================
HCoupling = g[0]*(a*sm0.dag() + sm0*a.dag()) + g[1]*(a*sm1.dag() + sm1*a.dag())
Hc = w_c * a.dag() * a 
H_eta = -eta_q[0] * E_uc0 - eta_q[1] * E_uc1
#Caculate the energy of qubit
HT = w_q[0]*sn0 + w_q[1]*sn1 + Hc + H_eta + HCoupling
w = HT.eigenenergies()
#w0 = (ptrace(HT,1)/(N*level)).eigenenergies();w1 = (ptrace(HT,2)/(N*level)).eigenenergies()
w01 = [(w[1]-w[0]),(w[2]-w[0])]
w02 = [(w[4]-w[0]),(w[6]-w[0])]