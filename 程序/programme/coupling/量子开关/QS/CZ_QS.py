#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 17:01:13 2017

@author: chen
CZ gate using quantum switch
"""



import time 
import csv
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from qutip import *
from scipy.optimize import *
from scipy.integrate import *
from scipy import interpolate 
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import *
from multiprocessing import Pool
from decimal import *
from math import *
import gc 
import sys



def CZ0(t,args):
    tp = args['tp']
    delta = args['delta']
    tlist = np.linspace(0,tp,resolution)
#    print(t,tlist[-1])
    w = interpolate.interp1d(tlist,delta*th,'slinear')
    if t<=tp and t>=0:
        w = w(t)
    else:
        w = 0
        
    return(w)
def plotCZ(t,args):
    tp = args['tp']
    delta = args['delta']
    tlist = np.linspace(0,tp,resolution)
    w = interpolate.interp1d(tlist,delta*th,'slinear')
    w = w(t)
    return(w)

def Z_wave(t,args):
    A = args['A']
    omega = args['omega']
    w = A*np.cos(omega*t)
    return(w)
    
    
    
def CZgate(P):
    tp  = P[0]
    delta = P[1]
    A = 0.05*2*1.25*2*np.pi
    omega = 0.05*2*np.pi
    args = {'tp':tp,'delta':delta,'A':A,'omega':omega}
    
    
    
    HCoupling = g[0]*(a+a.dag())*(sm[0]+sm[0].dag()) + g[1]*(a+a.dag())*(sm[1]+sm[1].dag())
    Hc = w_c * a.dag() * a 
    H_eta = -eta_q[0] * E_uc[0] - eta_q[1] * E_uc[1]
    Hq = w_q[0]*sn[0] + w_q[1]*sn[1]
    H0 = Hq + H_eta + Hc + HCoupling
    
    Hd0 = [sn[0],CZ0]
    Hd1 = [sn[1],Z_wave]
    H = [H0,Hd1]
    
#    geff = 0.5*g[1]*g[0]*(1/abs(w_q[1]+eta_q[0]-w_c)+1/abs(w_q[1]-w_c))
    
    [E,S] = H0.eigenstates()
    
#    tp = np.pi/np.sqrt(2)/geff
    tlist = np.linspace(0,tp,tp+1)
#    psi0 = (S[2]+S[5]).unit()
#    psi0 = tensor(basis(N,0),basis(3,1),(basis(3,0)+basis(3,1)).unit())
    psi0 = tensor(basis(N,0),basis(3,1),(basis(3,0)).unit())
    
    
    '''evolution'''
    
    options=Options()
    options.atol=1e-8
    options.rtol=1e-6
    options.first_step=0.01
    options.num_cpus= 4
    options.nsteps=1e6
    options.gui='True'
    options.ntraj=1000
    options.rhs_reuse=False
    start = time.time()
    
    result = mesolve(H,psi0,tlist,[],[],args = args,options = options)
    
    
    
    
    
    
    n_x0 = [] ; n_y0 = [] ; n_z0 = [];
    n_x1 = [] ; n_y1 = [] ; n_z1 = [];
    for t in range(0,len(tlist)):
        rf0 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E[3]-E[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E[2]-E[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        rf2 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E[1]-E[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
        U = tensor(qeye(3),rf1,rf2)
#        U = (1j*H0*tlist[t]).expm()
        
        opx0 = U.dag()*sx[0]*U
        opy0 = U.dag()*sy[0]*U
        opz0 = sz[0]
        opx1 = U.dag()*sx[1]*U
        opy1 = U.dag()*sy[1]*U
        opz1 = sz[1]
        n_x0.append(expect(opx0,result.states[t]))
        n_y0.append(expect(opy0,result.states[t]))
        n_z0.append(expect(opz0,result.states[t]))
        n_x1.append(expect(opx1,result.states[t]))
        n_y1.append(expect(opy1,result.states[t]))
        n_z1.append(expect(opz1,result.states[t]))

    
    fig, axes = plt.subplots(3, 1, figsize=(10,8))
    labels=['Q1X','Q1Y','Q1Z']  
    
    axes[0].plot(tlist, n_x0, label=labels[0]);axes[0].set_ylim([-1.1,1.1])
    axes[0].legend(loc=0);axes[0].set_xlabel('Time');axes[0].set_ylabel('P')
    axes[1].plot(tlist, n_y0, label=labels[1]);axes[1].set_ylim([-1.1,1.1])
    axes[1].legend(loc=0);axes[1].set_xlabel('Time');axes[1].set_ylabel('P')
    axes[2].plot(tlist, n_z0, label=labels[2]);axes[2].set_ylim([-1.1,1.1])
    axes[2].legend(loc=0);axes[2].set_xlabel('Time');axes[2].set_ylabel('P')
    
    fig, axes = plt.subplots(3, 1, figsize=(10,8))
    labels=['Q2X','Q2Y','Q2Z']  
    
    axes[0].plot(tlist, n_x1, label=labels[0]);axes[0].set_ylim([-1.1,1.1])
    axes[0].legend(loc=0);axes[0].set_xlabel('Time');axes[0].set_ylabel('P')
    axes[1].plot(tlist, n_y1, label=labels[1]);axes[1].set_ylim([-1.1,1.1])
    axes[1].legend(loc=0);axes[1].set_xlabel('Time');axes[1].set_ylabel('P')
    axes[2].plot(tlist, n_z1, label=labels[2]);axes[2].set_ylim([-1.1,1.1])
    axes[2].legend(loc=0);axes[2].set_xlabel('Time');axes[2].set_ylabel('P')
    
    
#    sphere = Bloch()
#    sphere.add_points([n_x0 , n_y0 , n_z0])
#    sphere.add_vectors([n_x0[-1],n_y0[-1],n_z0[-1]])
#    sphere.make_sphere() 
#    sphere = Bloch()
#    sphere.add_points([n_x1 , n_y1 , n_z1])
#    sphere.add_vectors([n_x1[-1],n_y1[-1],n_z1[-1]])
#    sphere.make_sphere() 
#    plt.show()
    
    rf0 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E[3]-E[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
    rf1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E[2]-E[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
    rf2 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E[1]-E[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
    UT = tensor(qeye(3),rf1,rf2)
#    target = (S[2]-S[5]).unit()
    target = tensor(basis(N,0),basis(3,1),(basis(3,0)-basis(3,1)).unit())


    fid = fidelity(UT*result.states[-1]*result.states[-1].dag()*UT.dag(),target)
    end = time.time()
    print(fid,P,end-start)
    
    return(1-fid)
    
    
def opt(gt):
    g[0] = gt;g[1] = gt
    geff = 0.5*g[1]*g[0]*(1/abs(w_q[1]+eta_q[0]-w_c)+1/abs(w_q[1]-w_c))
    x0 = [np.pi/np.sqrt(2)/geff+15,(w_q[1]-w_q[0]+eta_q[0])]
    result = minimize(CZgate, x0, method="Nelder-Mead",options={'disp': True})
    fid = 1-CZgate([result.x[0],result.x[1]])
    print(fid,gt,result.x[0],result.x[1])
    return([fid,[result.x[0],result.x[1]]])
        

if __name__=='__main__':
    thf = 0.55*np.pi/2;
    thi = 0.05;
    lam2 = -0.18;
    lam3 = 0.044;
    resolution = 1024
    ti=np.linspace(0,1,resolution)
    han2 = np.vectorize(lambda ti:(1-lam3)*(1-cos(2*pi*ti))+lam2*(1-cos(4*pi*ti))+lam3*(1-cos(6*pi*ti)))
    han2 = han2(ti)
    thsl=thi+(thf-thi)*han2/max(han2)
    x = 1/np.tan(thsl);
    x = x-x[0];
    
    tlu = np.cumsum(np.sin(thsl))*ti[1]
    tlu=tlu-tlu[0]
    ti=np.linspace(0, tlu[-1], resolution)
    th=interpolate.interp1d(tlu,thsl,'slinear')
    th = th(ti)
    th=1/np.tan(th)
    th=th-th[0]
    th=th/min(th)
    
    
    #==============================================================================
    global g
    ## Cavity Frequency
    w_c= 7.0 * 2 * np.pi
    
    ## Qubits frequency
    w_q = np.array([ 4.7 , 4.70 , 6 , 5.9 , 5.8 , 5.7 , 5.6 , 5.5 , 5.4 , 5.3]) * 2 * np.pi
    
    ## Coupling Strength
    g = np.array([0.1, 0.1, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]) * 2 * np.pi
    
    ## Qubits Anharmonicity
    
    eta_q=  np.array([0.250, 0.250, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]) * 2 * np.pi
    N = 3
    
    
    
    #==============================================================================
    
    a = tensor(destroy(N),qeye(3),qeye(3))
    sm = np.array([tensor(qeye(N),destroy(3),qeye(3)) , tensor(qeye(N),qeye(3),destroy(3))])
    
    E_uc = np.array([tensor(qeye(N),basis(3,2)*basis(3,2).dag(),qeye(3)) , tensor(qeye(N),qeye(3), basis(3,2)*basis(3,2).dag())])
    #用以表征非简谐性的对角线最后一�?非计算能�?
    #E_uc1 = tensor(qeye(N),qeye(3), Qobj([[0,0],[0,1]]))
    
    E_e = np.array([tensor(qeye(N),basis(3,1)*basis(3,1).dag(),qeye(3)),tensor(qeye(N),qeye(3),basis(3,1)*basis(3,1).dag())])
    #激发�?    
    E_g = np.array([tensor(qeye(N),basis(3,0)*basis(3,0).dag(),qeye(3)) , tensor(qeye(N),qeye(3),basis(3,0)*basis(3,0).dag())])
    #基�?    
    sn = np.array([sm[0].dag()*sm[0] , sm[1].dag()*sm[1]])
    
    sx = np.array([sm[0].dag()+sm[0],sm[1].dag()+sm[1]]);
    sxm = np.array([tensor(qeye(N),Qobj([[0,1,0],[1,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(N),qeye(3),Qobj([[0,1,0],[1,0,0],[0,0,0]]))])
    
    
    sy = np.array([1j*(sm[0].dag()-sm[0]) , 1j*(sm[1].dag()-sm[1])]);
    sym = np.array([tensor(qeye(N),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(N),qeye(3),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]))])
    
    sz = np.array([E_g[0] - E_e[0] , E_g[1] - E_e[1]])
    #==============================================================================
    starttime  = time.time()
    
#    geff = 0.5*g[1]*g[0]*(1/abs(w_q[1]+eta_q[0]-w_c)+1/abs(w_q[1]-w_c))
#    x0 = [np.pi/np.sqrt(2)/geff+15,(w_q[1]-w_q[0]+eta_q[0])]
#    result = minimize(CZgate, x0, method="Nelder-Mead",options={'disp': True})
#    print(result.x[0],result.x[1])


#    gt = np.linspace(0.08,0.012,41)
#    
#    p = Pool(16)
#    
#    A = p.map(opt,gt)
#    p.close()
#    p.join()
#    z =  np.array([x[0] for x in A])
#    para = np.array([x[1] for x in A])
#    figure()
#    plot(gt,z);xlable('g')
    
#    tp  = result.x[0]
#    delta = result.x[1]
#    tp  = 49
#    delta = 1.272
#    args = {'tp':tp,'delta':delta}
#    tlist = np.linspace(0,tp,tp+1)
#    D = plotCZ(tlist,args)/2/np.pi;
##    figure();
#    plot(tlist,D)
#    title(str(w_c/2/np.pi))

    CZgate([121.8612,-2.2147])
    
    
    endtime  = time.time()

    print('used time:',endtime-starttime,'s')