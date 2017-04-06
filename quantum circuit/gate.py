#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 19:48:46 2017

@author: chen
"""
from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from initialsetting import *
#==============================================================================
'''
single qubit gate
'''
def Gate_rx(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):
#    print(setting.DRAG)
    args_i={}
    w_t_i='(w_t'+str(inxc)+')*('+str(t0)+'<t<='+str(t1)+')'    
    if setting.DRAG:
#        print(setting.DRAG)
        D_t_i='(Omega'+str(inxc)+'*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+')-(t-20-'+str(t0)+')/2.0/width'+str(inxc)+'**2/'+str(setting.eta_q[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2.0)))*('+str(t0)+'<t<='+str(min(t1,t0+500))+')'
#        D_t_i = 'Omega' + str(inxc) + '*np.cos(t*f'+str(inxc)+')'
    else:
        D_t_i='(Omega'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    args_i['f'+str(inxc)]=setting.w_q[inxc]
    args_i['Omega'+str(inxc)]=setting.omega*2*phi
    args_i['width'+str(inxc)]=width
    return w_t_i,D_t_i,args_i     

def Gate_ry(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):
    args_i={}
    w_t_i='(w_t'+str(inxc)+')*('+str(t0)+'<t<='+str(t1)+')'
    if setting.DRAG:
        D_t_i='(Omega'+str(inxc)+'*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2)-(t-20-'+str(t0)+')/2.0/width'+str(inxc)+'**2/'+str(setting.eta_q[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi)))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    else:
        D_t_i='(Omega'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    args_i['f'+str(inxc)]=setting.w_q[inxc]
    args_i['Omega'+str(inxc)]=setting.omega*2*phi
    args_i['width'+str(inxc)]=width
    return w_t_i,D_t_i,args_i   
def Gate_rz(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):

    args_i={}
    w_t_i='(w_t'+str(inxc)+'+delta'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2))*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    args_i['width'+str(inxc)]=width
    args_i['delta'+str(inxc)]=setting.delta*phi
    return w_t_i,D_t_i,args_i

def Gate_i(inxc,t0,t1,setting = qusetting()):
    args_i={}
    w_t_i='(w_t'+str(inxc)+')*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    return w_t_i,D_t_i,args_i
#==============================================================================

#==============================================================================
'''
two qubits gate
'''
def Gate_iSWAP(inxc,inxt,t0,t1,phi=np.pi,width = 0.5,setting = qusetting()):

    args_i={}
    w_t_i='(w_t'+str(inxc)+'+delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'0)/width'+str(inxc)+')) -delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'1)/width'+str(inxc)+')))*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    args_i['width'+str(inxc)]=width
    args_i['t'+str(inxc)+'0']=2
    args_i['t'+str(inxc)+'1']=setting.iswap_deltat
    args_i['delta'+str(inxc)]=setting.w_q[inxt]-setting.w_q[inxc]
#    print(setting.iswap_deltat)
    return w_t_i,D_t_i,args_i

def Gate_CZ(inxc,inxt,t0,t1,phi=np.pi,width = 10,setting = qusetting()):

    args_i={}
    w_t_i='(w_t'+str(inxc)+'+delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'0)/width'+str(inxc)+')) -delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'1)/width'+str(inxc)+')))*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    args_i['width'+str(inxc)]=width
    args_i['t'+str(inxc)+'0']=100
    args_i['t'+str(inxc)+'1']=setting.CZ_deltat
    args_i['delta'+str(inxc)]=(setting.w_q[inxt]-setting.w_q[inxc]+setting.eta_q[inxt])
    
#    print(w_t_i)
#    print(args_i)
    return w_t_i,D_t_i,args_i


#==============================================================================