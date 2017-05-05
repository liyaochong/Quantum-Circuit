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
from scipy.special import *
#==============================================================================
'''
single qubit gate
'''
def Gate_rx(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):#X门，通过改Phi可以更改旋转角度
    args_i={}
    w_t_i='(0)*('+str(t0)+'<t<='+str(t1)+')'    
    if setting.DRAG:
        D_t_i='(Omega'+str(inxc)+'*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+')'
        D_t_i += '-(t-20-'+str(t0)+')/2.0/width'+str(inxc)+'**2/'+str(setting.eta_q[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2.0)))'
        D_t_i+='*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    else:
        D_t_i='(Omega'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    args_i['f'+str(inxc)]=setting.En[inxc]
    args_i['Omega'+str(inxc)]=setting.omega*2*phi
    args_i['width'+str(inxc)]=width
    return w_t_i,D_t_i,args_i     

def Gate_ry(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):#Y门，通过改Phi可以更改旋转角度
    args_i={}
    w_t_i='(0)*('+str(t0)+'<t<='+str(t1)+')'
    if setting.DRAG:
        D_t_i='(Omega'+str(inxc)+'*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2)'
        D_t_i+='-(t-20-'+str(t0)+')/2.0/width'+str(inxc)+'**2/'+str(setting.eta_q[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi)))'
        D_t_i+='*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    else:
        D_t_i='(Omega'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    args_i['f'+str(inxc)]=setting.En[inxc]
    args_i['Omega'+str(inxc)]=setting.omega*2*phi
    args_i['width'+str(inxc)]=width
    return w_t_i,D_t_i,args_i   
def Gate_rz(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):#Z门，通过改Phi可以更改旋转角度

    args_i={}
    w_t_i='(delta'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2))*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['w_t'+str(inxc)]=setting.w_q[inxc]
    args_i['width'+str(inxc)]=width
    args_i['delta'+str(inxc)]=setting.delta*phi
    return w_t_i,D_t_i,args_i

def Gate_i(inxc,t0,t1,setting = qusetting()):#I门
    args_i={}
    w_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    return w_t_i,D_t_i,args_i

def Gate_H(inxc,t0,t1,phi=np.pi,width = 6,setting = qusetting()):#H门
    args_i={}
    w_t_i='(0)*('+str(t0)+'<t<='+str(t1)+')'    
    if setting.DRAG:
        D_t_i='(-Omega'+str(inxc)+'/2*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2)'
        D_t_i += '-(t-20-'+str(t0)+')/2.0/width'+str(inxc)+'**2/'+str(setting.eta_q[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi)))'
        D_t_i+='*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
        D_t_i+='+(Omega'+str(inxc)+'*(np.exp(-(t-20-'+str(t0)+'-40)**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+')'
        D_t_i += '-(t-20-'+str(t0)+'-40)/2.0/width'+str(inxc)+'**2/'+str(setting.eta_q[inxc])+'*np.exp(-(t-20-'+str(t0)+'-40)**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2.0)))'
        D_t_i+='*('+str(t0+40)+'<t<='+str(min(t1,t0+80))+')'
    else:
        D_t_i='(-Omega'+str(inxc)+'/2*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'+np.pi/2))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
        D_t_i+='+(Omega'+str(inxc)+'*np.exp(-(t-20-'+str(t0)+'-40)**2/2.0/width'+str(inxc)+'**2)*np.cos(t*f'+str(inxc)+'))*('+str(t0+40)+'<t<='+str(min(t1,t0+80))+')'
    args_i['f'+str(inxc)]=setting.En[inxc]
    args_i['Omega'+str(inxc)]=setting.omega*2*phi
    args_i['width'+str(inxc)]=width

    return w_t_i,D_t_i,args_i     
#==============================================================================

#==============================================================================
'''
two qubits gate
'''
def Gate_iSWAP(inxc,inxt,t0,t1,phi=np.pi,width = 0.5,setting = qusetting()):#iswap门，还未标定

    args_i={}
    w_t_i='(delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'0)/width'+str(inxc)+')) -delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'1)/width'+str(inxc)+')))*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['width'+str(inxc)]=width
    args_i['t'+str(inxc)+'0']=2
    args_i['t'+str(inxc)+'1']=setting.iswap_deltat
    args_i['delta'+str(inxc)]=setting.En[inxt]-setting.En[inxc]
    return w_t_i,D_t_i,args_i

def Gate_CZ(inxc,inxt,t0,t1,phi=np.pi,width = 10,setting = qusetting()):#CZ门

    args_i={}
    w_t_i='(delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'0)/width'+str(inxc)+')) -delta'+str(inxc)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'1)/width'+str(inxc)+')))*('+str(t0)+'<t<='+str(t1)+')'
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    args_i['width'+str(inxc)]=width
    args_i['t'+str(inxc)+'0']=50
    args_i['t'+str(inxc)+'1']=setting.CZ_deltat
    args_i['delta'+str(inxc)]=(setting.En[inxt]-setting.En[inxc]+setting.eta_q[inxt])
    return w_t_i,D_t_i,args_i

#def Gate_CZ(inxc,inxt,t0,t1,phi=np.pi,setting = qusetting()):#Eef波形的CZ门
#
#    args_i={}
#    w_t_i='delta'+str(inxc)+'/2*(erf((t-'+str(t0)+'-t'+str(inxc)+'ramp/2)/np.sqrt(2)/sigma'+str(inxc)+')-erf((t-'+str(t0)+'-t'+str(inxc)+'gate+t'+str(inxc)+'ramp/2)/np.sqrt(2)/sigma'+str(inxc)+'))*('+str(t0)+'<t<='+str(t1)+')'
#    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['t'+str(inxc)+'ramp']=setting.ramp
#    args_i['sigma'+str(inxc)]=setting.ramp/4/np.sqrt(2)
#    args_i['t'+str(inxc)+'gate']=t1-t0      
#    args_i['delta'+str(inxc)]=(setting.En[inxt]-setting.En[inxc]+setting.eta_q[inxt])
#    
#    return w_t_i,D_t_i,args_i


#==============================================================================