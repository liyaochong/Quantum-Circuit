#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 14:27:44 2017

@author: chen
"""

from scipy.special import jn,  jn_zeros
from scipy.optimize import *
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from pylab import *
import time

def Z_wave(t,args):
    A = args['A']
    omega = args['omega']
    w = A*np.cos(omega*t)
    return(w)

def optA(A):
    
    tlist= np.linspace(0,1000,1001)
    
    

    A = A
    print(A/2/np.pi,omega/2/np.pi,A/omega)
    H=[H0,[sm0.dag()*sm0 , Z_wave]]
    
    psi0=tensor(basis(N,1),basis(N,0))
    
    options=Options()
    args = {'A':A , 'omega':omega}
    output = mesolve(H,psi0,tlist,[],[],args = args,options = options)
    
    exp = [expect(sm0.dag() * sm0,output.states) , expect(sm1.dag() * sm1,output.states)]
    fig, axes = plt.subplots(2, 1, figsize=(10,8))
    labels=['Q0','Q1']
    for ii in range(0,2):   
        n_q = exp[ii]
        
        axes[ii].plot(tlist, n_q, label=labels[ii])
        
        axes[ii].set_ylim([-0.1,1.1])
        axes[ii].legend(loc=0)
        axes[ii].set_xlabel('Time')
        axes[ii].set_ylabel('P')
    
    Qmin = min(exp[0])
    Qmax = max(exp[0])
    print(A/2/np.pi , g*jn(0,A/omega)/2/np.pi,Qmax-Qmin)
    return(Qmax-Qmin)


if __name__=='__main__':
    starttime=time.time()
    
    N = 3
    
    g = 0.009 * 2 * np.pi
    wq= np.array([5.1 , 5.1 ]) * 2 * np.pi
    eta_q=  np.array([-0.250 , -0.250]) * 2 * np.pi
    
    

    sm0=tensor(destroy(N),qeye(N))
    sm1=tensor(qeye(N),destroy(N))
    E_uc0 = tensor(basis(3,2)*basis(3,2).dag() , qeye(3)) 
    E_uc1 = tensor(qeye(3) , basis(3,2)*basis(3,2).dag())
    
    
    H0= (wq[0]) * sm0.dag()*sm0 + (wq[1]) * sm1.dag()*sm1 + eta_q[0]*E_uc0 + eta_q[1]*E_uc1 + g * (sm0.dag()+sm0) * (sm1.dag()+sm1)
#    H0= (wq[0]) * sm0.dag()*sm0 + (wq[1]) * sm1.dag()*sm1 + g * (sm0.dag()+sm0) * (sm1.dag()+sm1)
    

    

    
    omega = 0.4*2*np.pi
    A = jn_zeros(0,1)*omega
                
#    x0 = [jn_zeros(0,1)*omega ]
#    result = minimize(optA, x0, method="Nelder-Mead",options={'disp': True})
#    print(result.x/2/np.pi)

    optA(A[0])
    
    
    finishtime=time.time()
    print( 'Time used: ', (finishtime-starttime), 's')