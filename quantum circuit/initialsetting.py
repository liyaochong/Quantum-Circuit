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
class qusetting(object):
    def __init__(self):
        '''
        gate
        '''
        self.DRAG = True
        self.omega = 0.033348628532133037
        self.delta = 0.069641896319
        self.iswap_deltat = 361
        self.CZ_deltat = 267.517517518
#==============================================================================
        '''
        initial
        '''
        self.qtype = 1
#==============================================================================
        '''
        qubit
        '''
        self.w_c = 7.0  * 2 * np.pi  # cavity frequency
        self.w_q = np.array([ 5.0 , 6.0]) * 2 * np.pi
        self.g = np.array([0.06 , 0.06]) * 2 * np.pi
        self.eta_q = np.array([0.2 , 0.2]) * 2 * np.pi
        self.N = 3              # number of cavity fock states
        self.n= 0
        self.level = 3  #能级数
#==============================================================================
        '''
        Plot
        '''
        self.RF = True
#==============================================================================



def initial(setting = qusetting()):
    if setting.qtype == 1:#two qubits,coupled with a resonator

        
#==============================================================================
        '''
        Operators
        '''
        a = tensor(destroy(setting.N),qeye(3),qeye(3))
        sm = np.array([tensor(qeye(setting.N),destroy(3),qeye(3)) , tensor(qeye(setting.N),qeye(3),destroy(3))])
        
        E_uc = np.array([tensor(qeye(setting.N),basis(3,2)*basis(3,2).dag(),qeye(3)) , tensor(qeye(setting.N),qeye(3), basis(3,2)*basis(3,2).dag())])
        #用以表征非简谐性的对角线最后一项(非计算能级)
        #E_uc1 = tensor(qeye(N),qeye(3), Qobj([[0,0],[0,1]]))
        
        E_e = np.array([tensor(qeye(setting.N),basis(3,1)*basis(3,1).dag(),qeye(3)),tensor(qeye(setting.N),qeye(3),basis(3,1)*basis(3,1).dag())])
        #激发态
        
        E_g = np.array([tensor(qeye(setting.N),basis(3,0)*basis(3,0).dag(),qeye(3)) , tensor(qeye(setting.N),qeye(3),basis(3,0)*basis(3,0).dag())])
        #基态
        
        sn = np.array([sm[0].dag()*sm[0] , sm[1].dag()*sm[1]])
        
        sx = np.array([sm[0].dag()+sm[0],sm[1].dag()+sm[1]]);
        sxm = np.array([tensor(qeye(setting.N),Qobj([[0,1,0],[1,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(setting.N),qeye(3),Qobj([[0,1,0],[1,0,0],[0,0,0]]))])
        
        
        sy = np.array([-1j*(sm[0].dag()-sm[0]) , -1j*(sm[1].dag()-sm[1])]);
        sym = np.array([tensor(qeye(setting.N),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(setting.N),qeye(3),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]))])
        
        sz = np.array([E_g[0] - E_e[0] , E_g[1] - E_e[1]])
        
#==============================================================================
        HCoupling = setting.g[0]*(a+a.dag())*(sm[0]+sm[0].dag()) + setting.g[1]*(a+a.dag())*(sm[1]+sm[1].dag())
        Hc = setting.w_c * a.dag() * a 
        H_eta = setting.eta_q[0] * E_uc[0] + setting.eta_q[1] * E_uc[1]
        #Caculate the energy of qubit
        HT = setting.w_q[0]*sn[0] + setting.w_q[1]*sn[1] + Hc + H_eta + HCoupling
        #HT = w_q[0]*sn0 + w_q[1]*sn1 + Hc  + HCoupling
        w = HT.eigenenergies()
        #w0 = (ptrace(HT,1)/(setting.N*3)).eigenenergies();w1 = (ptrace(HT,2)/(setting.N*3)).eigenenergies()
        w01 = [(w[1]-w[0]),(w[2]-w[0])]
        w02 = [(w[4]-w[0]),(w[6]-w[0])]
        #print(w[2]-w[0],w[5]-w[1])
        
        return(a,sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,w01,w02)
        
    elif setting.qtype == 2:#two qubits,coupled directly

#==============================================================================
        '''
        Operators
        '''
        sm = np.array([tensor(destroy(3),qeye(3)) , tensor(qeye(3),destroy(3))])
        
        E_uc = np.array([tensor(basis(3,2)*basis(3,2).dag(),qeye(3)) , tensor(qeye(3), basis(3,2)*basis(3,2).dag())])
        #用以表征非简谐性的对角线最后一项(非计算能级)
        #E_uc1 = tensor(qeye(3), Qobj([[0,0],[0,1]]))
        
        E_e = np.array([tensor(basis(3,1)*basis(3,1).dag(),qeye(3)),tensor(qeye(3),basis(3,1)*basis(3,1).dag())])
        #激发态
        
        E_g = np.array([tensor(basis(3,0)*basis(3,0).dag(),qeye(3)) , tensor(qeye(3),basis(3,0)*basis(3,0).dag())])
        #基态
        
        sn = np.array([sm[0].dag()*sm[0] , sm[1].dag()*sm[1]])
        
        sx = np.array([sm[0].dag()+sm[0],sm[1].dag()+sm[1]]);
        sxm = np.array([tensor(Qobj([[0,1,0],[1,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(3),Qobj([[0,1,0],[1,0,0],[0,0,0]]))])
        
        
        sy = np.array([-1j*(sm[0].dag()-sm[0]) , -1j*(sm[1].dag()-sm[1])]);
        sym = np.array([tensor(Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(3),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]))])
        
        sz = np.array([E_g[0] - E_e[0] , E_g[1] - E_e[1]])
        
#==============================================================================
        Delta1 = abs(setting.w_q[0]-setting.w_c)
        Delta2 = abs(setting.w_q[1]-setting.w_c)
        g_effect = 0.5*setting.g[0]*setting.g[1]*(Delta1+Delta2)/(Delta1*Delta2)
        HCoupling = g_effect*(sm[0]+sm[0].dag())*(sm[1]+sm[1].dag())
        H_eta = setting.eta_q[0] * E_uc[0] + setting.eta_q[1] * E_uc[1]
        #Caculate the energy of qubit
        HT = setting.w_q[0]*sn[0] + setting.w_q[1]*sn[1] + H_eta + HCoupling
        #HT = w_q[0]*sn0 + w_q[1]*sn1 + Hc  + HCoupling
        w = HT.eigenenergies()
        #w0 = (ptrace(HT,1)/(N*level)).eigenenergies();w1 = (ptrace(HT,2)/(N*level)).eigenenergies()
        w01 = [(w[1]-w[0]),(w[2]-w[0])]
        w02 = [(w[4]-w[0]),(w[6]-w[0])]
        #print(w[2]-w[0],w[5]-w[1])
    
        
        return(sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,w01,w02)
        
        
        
        
        
    