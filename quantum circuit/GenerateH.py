#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 20:18:52 2017

@author: chen
"""
'''
Generator a Hamilton
'''

from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from gate import *
from initialsetting import *

def SingleOpTime(Operator,setting = gatesetting()):
    lenc = len(Operator)
    time = 0
    for inxc in range(0,lenc):
        if Operator[inxc][0]=='X':
            time = np.max([time , 40])
        elif Operator[inxc][0]=='Y':
            time = np.max([time , 40])
        elif Operator[inxc][0]=='Z':
            time = np.max([time , 40])
        elif Operator[inxc][0]=='i':
            time = np.max([time , setting.iswap_deltat+3])
        elif Operator[inxc][0]=='C':
            if Operator[inxc][1]=='Z':
                time = np.max([time , setting.CZ_deltat+100])
        elif Operator[inxc][0]=='I':
            if len(Operator[inxc]) == 1:
                time = np.max([time , 0])
            else:
                time = np.max([time , int(Operator[inxc][1:])])
    return(time)
                
                
                
                
def GenerateH(Operator , setting = gatesetting()):
    t0 = 0 ;  t1 = SingleOpTime(Operator,setting)
    lenc = len(Operator)
    HCoupling = g[0]*(a*sm0.dag() + sm0*a.dag()) + g[1]*(a*sm1.dag() + sm1*a.dag())
    Hc = w_c * a.dag() * a 
    H_eta = -eta_q[0] * E_uc0 - eta_q[1] * E_uc1
    H= Hc + H_eta + HCoupling
    
    lenc = len(Operator)
    args = {}
    w_t=[]
    D_t=[]
    for inxc in range(0,lenc):# import drive
        if Operator[inxc][0]=='X':
            w_t_i,D_t_i,args_i=Gate_rx(inxc,t0,t1,setting = setting)
        elif Operator[inxc][0]=='Y':
            w_t_i,D_t_i,args_i=Gate_ry(inxc,t0,t1,setting = setting)
        elif Operator[inxc][0]=='Z':
            w_t_i,D_t_i,args_i=Gate_rz(inxc,t0,t1)
        elif Operator[inxc][0]=='i':
            w_t_i,D_t_i,args_i=Gate_iSWAP(inxc,int(Operator[inxc][5:]),t0,t1,setting = setting)
        elif Operator[inxc][0]=='C':
            if Operator[inxc][1]=='Z':
                w_t_i,D_t_i,args_i=Gate_CZ(inxc,int(Operator[inxc][2:]),t0,t1,setting = setting)
        elif Operator[inxc][0]=='I':
                w_t_i,D_t_i,args_i=Gate_i(inxc,t0,t1)
        w_t.append(w_t_i)
        D_t.append(D_t_i)
        args=dict(args,**args_i)##在args里加入args_i
    H_q0 = [sn0 , w_t[0]]
    H_q1 = [sn1 , w_t[1]]
    H_d0 = [sx0 , D_t[0]] 
    H_d1 = [sx1 , D_t[1]]
    H = [H , H_q0 , H_q1 , H_d0 , H_d1]
    tlist = np.linspace(0,t1,2*(t1-0)+1)
    return H , args , tlist 