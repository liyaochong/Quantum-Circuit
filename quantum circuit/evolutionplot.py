#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 21:00:37 2017

@author: chen
"""
from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from initialsetting import *


#==============================================================================
class plotsetting(object):
    def __init__(self):
        self.RF = True
#==============================================================================

#==============================================================================
def evolutionplot(target , result , tlist , setting = plotsetting()):
    
    n_x = [] ; n_y = [] ; n_z = [];
    for t in range(0,len(tlist)):
        if setting.RF:
#            -(g[1]*g[1]/(w_q[1]-w_q[0]))
            rf01 = np.exp(1j*(w02[0])*tlist[t])*basis(3,2)*basis(3,2).dag()
            rf02 = np.exp(1j*(w01[0])*tlist[t])*basis(3,1)*basis(3,1).dag()
            rf0 = basis(3,0)*basis(3,0).dag()+rf01+rf02
            rf11 = np.exp(1j*(w02[1])*tlist[t])*basis(3,2)*basis(3,2).dag()
            rf12 = np.exp(1j*(w01[1])*tlist[t])*basis(3,1)*basis(3,1).dag()
            rf1 = basis(3,0)*basis(3,0).dag()+rf11+rf12
            U = tensor(qeye(N),rf0,rf1)
            opx = U.dag()*eval('sxm'+str(eval('target')))*U
            opy = U.dag()*eval('sym'+str(eval('target')))*U
            opz = U.dag()*eval('sz'+str(eval('target')))*U
            
        else:
            opx = eval('sxm'+str(eval('target')))
            opy = eval('sym'+str(eval('target')))
            opz = eval('sz'+str(eval('target')))
            
    
        n_x.append(expect(opx,result.states[t]))
        n_y.append(expect(opy,result.states[t]))
        n_z.append(expect(opz,result.states[t]))
    fig, axes = plt.subplots(1, 3, figsize=(10,6))
    
    axes[0].plot(tlist, n_x, label='X')
    axes[1].plot(tlist, n_y, label='Y')
    axes[2].plot(tlist, n_z, label='Z')
    sphere = Bloch()
    sphere.add_points([n_x , n_y , n_z])
    sphere.add_vectors([n_x[-1],n_y[-1],n_z[-1]])
    sphere.make_sphere() 
    plt.show()
#==============================================================================
