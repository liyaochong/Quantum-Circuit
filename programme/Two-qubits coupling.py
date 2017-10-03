# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 15:40:22 2017

@author: Chen
"""
'''
两个qubit通过resonator进行耦合
'''
'''
事例是两个qubit有相同的频率，和resonator有个detuning，初态都在+态，观察是否有能量交换和相位积累
'''
from qutip import *
import numpy as np


from time import clock
starttime=clock()
'''
parameters
'''
wc = 7.0  * 2 * pi  # cavity frequency
wq = np.array([5.0,6.0]) * 2 * np.pi  # atom frequency
g  = np.array([0.03,0.03]) * 2 * np.pi  # coupling strength
delta = np.array([wq[0]-wc,wq[1]-wc])   #detuning
eta_q = np.array([-0.25 , -0.25]) * 2 * np.pi   #非简谐项
'''
Operators
'''
a = tensor(destroy(3),qeye(3),qeye(3))
sm = np.array([tensor(qeye(3),destroy(3),qeye(3)) , tensor(qeye(3),qeye(3),destroy(3))])
        
E_uc = np.array([tensor(qeye(3),basis(3,2)*basis(3,2).dag(),qeye(3)) , tensor(qeye(3),qeye(3), basis(3,2)*basis(3,2).dag())])
#用以表征非简谐性的对角线最后一项(非计算能级)
#E_uc1 = tensor(qeye(N),qeye(3), Qobj([[0,0],[0,1]]))

E_e = np.array([tensor(qeye(3),basis(3,1)*basis(3,1).dag(),qeye(3)),tensor(qeye(3),qeye(3),basis(3,1)*basis(3,1).dag())])
#激发态

E_g = np.array([tensor(qeye(3),basis(3,0)*basis(3,0).dag(),qeye(3)) , tensor(qeye(3),qeye(3),basis(3,0)*basis(3,0).dag())])
#基态

sn = np.array([sm[0].dag()*sm[0] , sm[1].dag()*sm[1]])

sx = np.array([sm[0].dag()+sm[0],sm[1].dag()+sm[1]]);
sxm = np.array([tensor(qeye(3),Qobj([[0,1,0],[1,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(3),qeye(3),Qobj([[0,1,0],[1,0,0],[0,0,0]]))])


sy = np.array([1j*(sm[0].dag()-sm[0]) , 1j*(sm[1].dag()-sm[1])]);
sym = np.array([tensor(qeye(3),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(3),qeye(3),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]))])

sz = np.array([E_g[0] - E_e[0] , E_g[1] - E_e[1]])

'''
Hamilton
'''
HCoupling = g[0]*(a+a.dag())*(sm[0]+sm[0].dag()) + g[1]*(a+a.dag())*(sm[1]+sm[1].dag())
Hc = w_c * a.dag() * a 
H_eta = eta_q[0] * E_uc[0] + eta_q[1] * E_uc[1]
Hq = w_q[0]*sn[0] + w_q[1]*sn[1]
H = Hq + H_eta + Hc + HCoupling











finishtime=clock()
print ('Time used: ', (finishtime-starttime), 's')
