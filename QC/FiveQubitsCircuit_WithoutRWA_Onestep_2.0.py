#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 04:04:14 2017

@author: qubits0
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 22:29:10 2017

@author: qubits0
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 21:14:06 2017

@author: qubits0
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 01 14:56:12 2016

@author: Logmay
"""

import time
import csv
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from qutip import *
from scipy.optimize import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import *
from decimal import *
from multiprocessing import Pool


def CreateBasicOperator(APRType=0):
    print('APRType:'+str(APRType))
    if APRType==0:
        cmdstr=''
        for II in range(0,Num_Q):
            cmdstr+='qeye(3),'
        a= eval('tensor(destroy(N),'+cmdstr+')')
    
        sm=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='destroy(3),'
                else:
                    cmdstr+='qeye(3),'
            sm.append(eval('tensor(qeye(N),'+cmdstr+')'))
    
        smm=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='basis(3,2)*basis(3,2).dag(),'
                else:
                    cmdstr+='qeye(3),'
            smm.append(eval('tensor(qeye(N),'+cmdstr+')'))
        
        E_e=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='basis(3,1)*basis(3,1).dag(),'
                else:
                    cmdstr+='qeye(3),'
            E_e.append(eval('tensor(qeye(N),'+cmdstr+')'))
        
        E_g=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='basis(3,0)*basis(3,0).dag(),'
                else:
                    cmdstr+='qeye(3),'
            E_g.append(eval('tensor(qeye(N),'+cmdstr+')'))
        
    elif APRType==1:
        cmdstr=''
        for II in range(0,Num_Q):
            cmdstr+='qeye(3),'
        a= eval('tensor('+cmdstr+')')
    
        sm=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='destroy(3),'
                else:
                    cmdstr+='qeye(3),'
            sm.append(eval('tensor('+cmdstr+')'))
    
        smm=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='basis(3,2)*basis(3,2).dag(),'
                else:
                    cmdstr+='qeye(3),'
            smm.append(eval('tensor('+cmdstr+')'))
        
        E_e=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='basis(3,1)*basis(3,1).dag(),'
                else:
                    cmdstr+='qeye(3),'
            E_e.append(eval('tensor('+cmdstr+')'))
        
        E_g=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='basis(3,0)*basis(3,0).dag(),'
                else:
                    cmdstr+='qeye(3),'
            E_g.append(eval('tensor('+cmdstr+')'))
        
    elif APRType==2:
        cmdstr=''
        for II in range(0,Num_Q):
            cmdstr+='qeye(2),'
        a= eval('tensor(destroy(N),'+cmdstr+')')
    
        sm=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='destroy(2),'
                else:
                    cmdstr+='qeye(2),'
            sm.append(eval('tensor(qeye(N),'+cmdstr+')'))
        
        smm=[]
        
        E_e=[]
        for II in range(0,Num_Q):
            E_e.append(sm[II].dag()*sm[II])
        E_g=[]
        for II in range(0,Num_Q):
            E_g.append(sm[II]*sm[II].dag())
    
    elif APRType==3:
        cmdstr=''
        for II in range(0,Num_Q):
            cmdstr+='qeye(2),'
        a= eval('tensor('+cmdstr+')')
    
        sm=[]
        for II in range(0,Num_Q):
            cmdstr=''
            for JJ in range(0,Num_Q):
                if II==JJ:
                    cmdstr+='destroy(2),'
                else:
                    cmdstr+='qeye(2),'
            sm.append(eval('tensor('+cmdstr+')'))
        
        smm=[]
        
        E_e=[]
        for II in range(0,Num_Q):
            E_e.append(sm[II].dag()*sm[II])
        E_g=[]
        for II in range(0,Num_Q):
            E_g.append(sm[II]*sm[II].dag())
    
    Sz=[]
    Sx=[]
    Sy=[]
    for II in range(0,Num_Q):
        Sz.append(sm[II].dag()*sm[II])
        Sx.append(sm[II].dag()+sm[II])
        Sy.append(1j*(sm[II].dag()-sm[II]))   
        
    return a, sm, smm, E_e, E_g, Sz, Sx, Sy

def GenerateDissipation(Dis=False):
    c_op_list = []
    
    if Dis==True:
        for II in range(0,Num_Q):
            c_op_list.append(np.sqrt(gamma[II] * (1+n_th)) * sm[II])
            c_op_list.append(np.sqrt(gamma[II] * n_th) * sm[II].dag())
            c_op_list.append(np.sqrt(gamma_phi[II]) * sm[II].dag()*sm[II])
            
    return c_op_list

def Gate_rx(inxc,t0,t1,index,phi=np.pi,omega=0.03332):
    if TG:
        omega=fx
    args_i={}
    w_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['w_t'+str(inxc)+'_'+str(index)]=0
    if DRAG:
        D_t_i='(Omega'+str(inxc)+'_'+str(index)+'*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2)*np.cos(t*f'+str(inxc)+'_'+str(index)+')+(t-20-'+str(t0)+')/2/width'+str(inxc)+'_'+str(index)+'**2/'+str(w_qa[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2)*np.cos(t*f'+str(inxc)+'_'+str(index)+'+np.pi/2)))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    else:
        D_t_i='(Omega'+str(inxc)+'_'+str(index)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2)*np.cos(t*f'+str(inxc)+'_'+str(index)+'))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    args_i['f'+str(inxc)+'_'+str(index)]=E_q[inxc]
    args_i['Omega'+str(inxc)+'_'+str(index)]=omega*2*phi
    args_i['width'+str(inxc)+'_'+str(index)]=6
    return w_t_i,D_t_i,args_i     
         
def Gate_ry(inxc,t0,t1,index,phi=np.pi,omega=0.03332):
    if TG:
        omega=fx
    args_i={}
    w_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['w_t'+str(inxc)+'_'+str(index)]=0
    if DRAG:
        D_t_i='(Omega'+str(inxc)+'_'+str(index)+'*(np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2)*np.cos(t*f'+str(inxc)+'_'+str(index)+'-np.pi/2)+(t-20-'+str(t0)+')/2/width'+str(inxc)+'_'+str(index)+'**2/'+str(w_qa[inxc])+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2)*np.cos(t*f'+str(inxc)+'_'+str(index)+')))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    else:
        D_t_i='(Omega'+str(inxc)+'_'+str(index)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2)*np.cos(t*f'+str(inxc)+'_'+str(index)+'-np.pi/2))*('+str(t0)+'<t<='+str(min(t1,t0+40))+')'
    args_i['f'+str(inxc)+'_'+str(index)]=E_q[inxc]
    args_i['Omega'+str(inxc)+'_'+str(index)]=omega*2*phi
    args_i['width'+str(inxc)+'_'+str(index)]=6 
    return w_t_i,D_t_i,args_i
def Gate_H(inxc,t0,t1,index):
    w_t_i1,D_t_i1,args_i1=Gate_ry(inxc,t0,t0+40,index+100000,phi=np.pi/2)
    w_t_i2,D_t_i2,args_i2=Gate_rx(inxc,t0+40,t1,index+110000,phi=np.pi)
    args_i=dict(args_i1,**args_i2)
    
    w_t_i=w_t_i1+'+'+w_t_i2
    D_t_i=D_t_i1+'+'+D_t_i2
    return w_t_i,D_t_i,args_i
#==============================================================================
# def Gate_rz(inxc,t0,t1,index,phi=np.pi):
#     args_i={}
#     w_t_i='(w_t'+str(inxc)+'_'+str(index)+'+delta'+str(inxc)+'_'+str(index)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'0)/width'+str(inxc)+'_'+str(index)+')) -delta'+str(inxc)+'_'+str(index)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'1)/width'+str(inxc)+'_'+str(index)+')))*('+str(t0)+'<t<='+str(t1)+')'
#     args_i['w_t'+str(inxc)+'_'+str(index)]=w_q[inxc]
#     args_i['width'+str(inxc)+'_'+str(index)]=0.5
#     args_i['t'+str(inxc)+'_'+str(index)+'0']=2
#     args_i['t'+str(inxc)+'_'+str(index)+'1']=38
#     args_i['delta'+str(inxc)+'_'+str(index)]=0.0293*phi
#     D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#     return w_t_i,D_t_i,args_i
#==============================================================================

def Gate_rz(inxc,t0,t1,index,phi=np.pi,delta=0.0688001 ):
    if TG:
        delta=fx
    args_i={}
    w_t_i='(delta'+str(inxc)+'_'+str(index)+'*np.exp(-(t-20-'+str(t0)+')**2/2.0/width'+str(inxc)+'_'+str(index)+'**2))*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['w_t'+str(inxc)+'_'+str(index)]=w_q[inxc]
    args_i['width'+str(inxc)+'_'+str(index)]=6
    args_i['delta'+str(inxc)+'_'+str(index)]=-delta*phi
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    return w_t_i,D_t_i,args_i

def Gate_iSWAP(inxc,inxt,t0,t1,index,phi=np.pi,deltat=500):
    if TG:
        deltat=fx
    args_i={}
    w_t_i='(delta'+str(inxc)+'_'+str(index)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'0)/width'+str(inxc)+'_'+str(index)+')) -delta'+str(inxc)+'_'+str(index)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'1)/width'+str(inxc)+'_'+str(index)+')))*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['w_t'+str(inxc)+'_'+str(index)]=w_q[inxc]
    args_i['width'+str(inxc)+'_'+str(index)]=0.5
    args_i['t'+str(inxc)+'_'+str(index)+'0']=2
    args_i['t'+str(inxc)+'_'+str(index)+'1']=deltat
    args_i['delta'+str(inxc)+'_'+str(index)]=E_q[inxt]-E_q[inxc]
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    return w_t_i,D_t_i,args_i

def Gate_CZ(inxc,inxt,t0,t1,index,phi=np.pi,deltat=751):
    if TG:
        deltat=fx
    args_i={}
    w_t_i='(delta'+str(inxc)+'_'+str(index)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'0)/width'+str(inxc)+'_'+str(index)+')) -delta'+str(inxc)+'_'+str(index)+'/(1 + np.exp(-(t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'1)/width'+str(inxc)+'_'+str(index)+')))*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['w_t'+str(inxc)+'_'+str(index)]=w_q[inxc]
    args_i['width'+str(inxc)+'_'+str(index)]=10
    args_i['t'+str(inxc)+'_'+str(index)+'0']=deltat
    args_i['t'+str(inxc)+'_'+str(index)+'1']=850-deltat
    args_i['delta'+str(inxc)+'_'+str(index)]=(E_q[inxt]-E_q[inxc]-w_qa[inxt])
#    print(E_q[inxt]-E_q[inxc]-w_qa[inxt])
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#    print(w_t_i)
#    print(args_i)
    return w_t_i,D_t_i,args_i


#def Gate_CZ(inxc,inxt,t0,t1,index,phi=np.pi,ramp=139):
#    if TG:
#        ramp=fx
#    args_i={}
#    w_t_i='delta'+str(inxc)+'_'+str(index)+'/2*(erf((t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'ramp/2)/np.sqrt(2)/sigma'+str(inxc)+'_'+str(index)+')-erf((t-'+str(t0)+'-t'+str(inxc)+'_'+str(index)+'gate+t'+str(inxc)+'_'+str(index)+'ramp/2)/np.sqrt(2)/sigma'+str(inxc)+'_'+str(index)+'))*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['t'+str(inxc)+'_'+str(index)+'ramp']=ramp
#    args_i['sigma'+str(inxc)+'_'+str(index)]=ramp/4/np.sqrt(2)
#    args_i['t'+str(inxc)+'_'+str(index)+'gate']=t1-t0      
#    args_i['delta'+str(inxc)+'_'+str(index)]=(E_q[inxt]-E_q[inxc]-w_qa[inxt])
#    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#    return w_t_i,D_t_i,args_i


def Gate_i(inxc,t0,t1,index):
    args_i={}
    w_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
#    args_i['w_t'+str(inxc)+'_'+str(index)]=0
    D_t_i='0*('+str(t0)+'<t<='+str(t1)+')'
    return w_t_i,D_t_i,args_i

def GetLen(Operator):
    lenc=len(Operator)
    LenSingleOperator=0
    for inxc in range(0,lenc):
        if Operator[inxc][0]=='X':
            LenSingleOperator=np.max([LenSingleOperator,40])
        elif Operator[inxc][0]=='Y':
            LenSingleOperator=np.max([LenSingleOperator,40])
        elif Operator[inxc][0]=='Z':
            LenSingleOperator=np.max([LenSingleOperator,40])
        elif Operator[inxc][0]=='H':
            LenSingleOperator=np.max([LenSingleOperator,80])    
        elif Operator[inxc][0]=='i':
            LenSingleOperator=np.max([LenSingleOperator,500])
        elif Operator[inxc][0]=='C':
            if Operator[inxc][1]=='Z':
                LenSingleOperator=np.max([LenSingleOperator,850])
        elif Operator[inxc][0]=='I':
            if len(Operator[inxc])==1:
                LenSingleOperator=np.max([LenSingleOperator,0])
            else:
                LenSingleOperator=np.max([LenSingleOperator,int(Operator[inxc][1:])])
            
    return LenSingleOperator

def GenerateH(Operators, APRType=0):
    global E_q
    
    if APRType==0 or APRType==2:
        HCoupling=0
        for II in range(0,Num_Q):
#            HCoupling+= g[II] * (a * sm[II].dag() + a.dag() * sm[II]) 
            HCoupling += g[II]*(a+a.dag())*(sm[II]+sm[II].dag())
        Hc=w_c * a.dag() * a 
        H=HCoupling+Hc
    elif APRType==1 or APRType==3:
        HCoupling=0
        for II in range(0,Num_Q):
            for JJ in range(II,Num_Q):
#                HCoupling+= g[II]*g[JJ]/(w_q[II]-w_c)/(w_q[JJ]-w_c)*(w_q[II]-w_c+w_q[JJ]-w_c) * (sm[II].dag() * sm[JJ])/2
                HCoupling+= g[II]*g[JJ]/(w_q[II]-w_c)/(w_q[JJ]-w_c)*(w_q[II]-w_c+w_q[JJ]-w_c) * (sm[II]+sm[II].dag())*(sm[JJ]+sm[JJ].dag())/2/2
        H=HCoupling
    
    if APRType==0 or APRType==1:    
        for II in range(0,Num_Q):
            H+=-w_qa[II]*smm[II].dag()*smm[II]
    
    for II in range(0,Num_Q):
        H+=(Sz[II])*(w_q[II])
    
    w_f = [w_q[k] for k in range(Num_Q)]
    w_f.append(w_c)
    
    En= (H.eigenenergies())

    E_n = En.tolist()

    E_index=sorted(range(Num_Q+1), key=lambda k: w_f[k])
    cloc = np.where(np.array(E_index) == Num_Q)[0][0]

    E_n.pop(cloc+1)
    E_index.pop(cloc)
    
    E_q = np.zeros(Num_Q)
    for idx , II in enumerate(E_index):
#        JJ = np.where(E_index == II)[0][0]
#        E_q.append(En[JJ+1]-En[0])
#        print(E_q[II]/2/np.pi)
        E_q[II] = E_n[idx+1]-E_n[0]
    print(E_q/2/np.pi)
    H=[H]
    
    t0=0
    args={}
    for JJ in range(0,len(Operators)):
        Operator=Operators[JJ]           
        lenc=Num_Q
        
        LenSingleOperator=GetLen(Operator)
        t1=t0+LenSingleOperator
        
        w_t=[]
        D_t=[]
        
        for inxc in range(0,lenc):
            if Operator[inxc][0]=='X':
                if len(Operator[inxc])==1:
                    w_t_i,D_t_i,args_i=Gate_rx(inxc,t0,t1,JJ)
                else:
                    w_t_i,D_t_i,args_i=Gate_rx(inxc,t0,t1,JJ,phi=np.pi*float(Operator[inxc][1:]))
            elif Operator[inxc][0]=='Y':
                if len(Operator[inxc])==1:
                    w_t_i,D_t_i,args_i=Gate_ry(inxc,t0,t1,JJ)
                else:
                    w_t_i,D_t_i,args_i=Gate_ry(inxc,t0,t1,JJ,phi=np.pi*float(Operator[inxc][1:]))
            elif Operator[inxc][0]=='Z':
                if len(Operator[inxc])==1:
                    w_t_i,D_t_i,args_i=Gate_rz(inxc,t0,t1,JJ)
                else:
                    w_t_i,D_t_i,args_i=Gate_rz(inxc,t0,t1,JJ,phi=np.pi*float(Operator[inxc][1:]))
            elif Operator[inxc][0]=='H':
                w_t_i,D_t_i,args_i=Gate_H(inxc,t0,t1,JJ)
            elif Operator[inxc][0]=='i':
                w_t_i,D_t_i,args_i=Gate_iSWAP(inxc,int(Operator[inxc][5:]),t0,t1,JJ)
            elif Operator[inxc][0]=='C':
                if Operator[inxc][1]=='Z':
                    w_t_i,D_t_i,args_i=Gate_CZ(inxc,int(Operator[inxc][2:]),t0,t1,JJ)
            elif Operator[inxc][0]=='I':
                w_t_i,D_t_i,args_i=Gate_i(inxc,t0,t1,JJ)
            w_t.append(w_t_i)
            D_t.append(D_t_i)
            args=dict(args,**args_i)
        
        W=[]
        D=[]
        for II in range(0,Num_Q):
            W.append([Sz[II],w_t[II]])   
            D.append([Sx[II],D_t[II]])
            H.append(W[II])
            H.append(D[II])
            
        t0=t1
        
#        print(H[1].eigenenergies()/2/np.pi)
#    print(H)  
#    print(args)  
    return H,args,t1,E_q
    

def plotstate(states,tlist,E_q,PlotPts=1000):
    
    if np.size(tlist)>PlotPts*3:
        states=states[0:-1:int(np.size(tlist)/PlotPts)]
        tlist=tlist[0:-1:int(np.size(tlist)/PlotPts)]
    
    expect_z=[]
    fig, axes = plt.subplots(max(Num_Q,2), 1, figsize=(10,8))
    for ii in range(0,Num_Q):
        n_z = np.real(expect(1-Sz[ii]*2,states))
#        n_z = np.real(expect(Sz[ii],states))
        n_z=n_z.tolist()
        axes[ii].plot(tlist, n_z, label='Q'+str(ii)+'Z')
        axes[ii].set_ylim([-1.05,1.05])
        axes[ii].legend(loc=0)
        axes[ii].set_xlabel('Time')
        axes[ii].set_ylabel('P')
        plt.show()
        
        expect_z.append(n_z)
    
    U=[]
    for ii in range(0,Num_Q):
        for t in range(0,len(tlist)):
            U.append((np.exp(1j*E_q[ii]*tlist[t])*E_e[ii]+E_g[ii]).dag())
    
    expect_x=[]
    fig, axes = plt.subplots(max(Num_Q,2), 1, figsize=(10,8))
    for ii in range(0,Num_Q):
        n_x=[]
        U0 = U[(ii)*len(tlist)]
        for t in range(0,len(tlist)):
            
            Ui=U[(ii)*len(tlist)+t]
            op=Ui*Sx[ii]*Ui.dag()
#            op = Sx[ii]
            n_x.append(expect(op,states[t]))
        axes[ii].plot(tlist, n_x, label='Q'+str(ii)+'X')
        axes[ii].set_ylim([-1.05,1.05])
        axes[ii].legend(loc=0)
        axes[ii].set_xlabel('Time')
        axes[ii].set_ylabel('P')
        plt.show()
        
        expect_x.append(n_x)
    
    expect_y=[]
    fig, axes = plt.subplots(max(Num_Q,2), 1, figsize=(10,8))
    for ii in range(0,Num_Q):
        n_y=[]
        U0 = U[(ii)*len(tlist)]
        for t in range(0,len(tlist)):
            
            Ui=U[(ii)*len(tlist)+t]
            op=Ui*Sy[ii]*Ui.dag()
#            op = Sy[ii]
            n_y.append(expect(op,states[t]))
        axes[ii].plot(tlist, n_y, label='Q'+str(ii)+'Y')
        axes[ii].set_ylim([-1.05,1.05])
        axes[ii].legend(loc=0)
        axes[ii].set_xlabel('Time')
        axes[ii].set_ylabel('P')
        plt.show()
        
        expect_y.append(n_y)
    
    for ii in range(0,Num_Q):
        b=Bloch()
        b.add_points([expect_x[ii],expect_y[ii],expect_z[ii]])
        b.add_vectors([expect_x[ii][-1],expect_y[ii][-1],expect_z[ii][-1]])
        b.show()
        if SD:
            b.save('./saved results/Bloch_Q'+str(ii),'png')
    pltlist = tlist
    return expect_x,expect_y,expect_z,pltlist

def savedata(psi_total, tlist_total, expect_x,expect_y,expect_z):
    for ii in range(0,Num_Q):
        data_xyz=np.array([tlist_total.tolist(),expect_x[ii],expect_y[ii],expect_z[ii]])
        file_data_store('./saved results/data_xyz_Q'+str(ii)+'.csv',data_xyz,sep=',',numformat='exp')
    file_data_store('./saved results/data_dm.csv',psi_total[-1],sep=',',numformat='exp')

def MCaverage(inputstates):   
    states = inputstates
    outputstate = states-states
    lenc = len(states)
    lenr = len(states[0])
    for j in range(0,lenc):
        states1=states[j]
        for ii in range(0,lenr):
            states1[ii]=ket2dm(states[j][ii])
        outputstate[j] = np.mean(states1)
    
    return outputstate


def Measure(state,target):
    exp1 = expect(E_e[target],state)    #probability of excited state
    exp0 = expect(E_g[target],state)    #probability of ground state
    measurestate1 = E_e[target]*state*state.dag()*E_e[target].dag()/exp1 #state with 1
    measurestate0 = E_g[target]*state*state.dag()*E_g[target].dag()/exp0 #state with 0
                       
    return([[exp0,measurestate0],[exp1,measurestate1]])
    
def EvolveCircuits(Operators,psi):
    
    #List of Gates: 'I','X','Y', 'Z', 'H', 'CX', 'CY', 'CZ', 'CNOT', 'SWAP', 'iSWAP', 'CSWAP', 'CCNOT', 'sSWAP', 'siSWAP', 'sNOT'
    if APRType==0 or APRType==2:
        if np.size(psi0.dims[0])>Num_Q+1:
            lst=np.arange(0,Num_Q+1).tolist()
        elif np.size(psi0.dims[0])==Num_Q+1:
            lst=[]
        elif np.size(psi0.dims[0])<Num_Q+1:
            raise NameError('Initial state size too short!')
    else:
        if np.size(psi0.dims[0])>Num_Q:
            lst=np.arange(0,Num_Q).tolist()
        elif np.size(psi0.dims[0])==Num_Q:
            lst=[]
        elif np.size(psi0.dims[0])<Num_Q:
            raise NameError('Initial state size too short!')
    if lst!=[]:
        psi=psi.ptrace(lst)

    H,args,t1,E_q=GenerateH(Operators,APRType)
    
#    print(psi)
#    print(H)
#    print(args)
#    print(c_op_list)    
#    print(args)
    tlist=np.arange(0,t1+1)
#    tlist = np.linspace(0,t1,t1+1)
    if MC:
        output=mcsolve(H,psi,tlist,c_op_list,[],args=args,options=options)
        states=MCaverage2(output.states)
    else:
        output=mesolve(H,psi,tlist,c_op_list,[],args=args,options=options)
        states=output.states
##    print(ptrace(states[-1],1))
#    print(ptrace(states[250],1))
#    psi_target=tensor(basis(N,0) ,basis(3,1),basis(3,0))
#    psi_target=tensor(basis(N,0), (basis(3,0)-basis(3,1)).unit(),basis(3,1))
#    tar = (tensor(basis(3,0),basis(3,0))+tensor(basis(3,0),basis(3,1))+tensor(basis(3,1),basis(3,0))-tensor(basis(3,1),basis(3,1))).unit()
#    psi_target=tensor(basis(N,0),tar)
#
#    U0 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E_q[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()
#    U1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E_q[1])*tlist[-1])*basis(3,1)*basis(3,1).dag()
#    UT = tensor(qeye(3),U0,U1)
#
#    fid=fidelity(UT*states[-1]*states[-1].dag()*UT.dag(), psi_target)
#    print(fid)
    
#    print(ptrace(UT*states[-1]*UT.dag(),1))
#    
#    U0 = np.exp(1j*E_q[0]*tlist[-1])*E_e[0]+E_g[0]
#    fid=fidelity(ptrace(U0*states[-1],1), ptrace(psi_target[-1],1))
#    print(fid)
#    print('evolve end')
    expect_x=[]
    expect_y=[]
    expect_z=[]
    pltlist = []
    if PLT:
        expect_x,expect_y,expect_z,pltlist=plotstate(states,tlist,E_q)
    return states, tlist, expect_x,expect_y,expect_z,pltlist
    
    
def TestGate(w):
    #List of Gates: 'I','X','Y', 'Z', 'H', 'CX', 'CY', 'CZ', 'CNOT', 'SWAP', 'iSWAP', 'CSWAP', 'CCNOT', 'sSWAP', 'siSWAP', 'sNOT'
        
#    psi0=tensor(basis(N,0), (basis(3,0)+basis(3,1)).unit(),basis(3,1))
#    psi0=tensor(basis(N,0),basis(3,1), (basis(3,0)+basis(3,1)).unit())
    
    Operators=[
        ['CZ1','I','I','I','I'],
        ] 
    
#    psi_target=tensor(basis(N,0), (basis(3,0)-basis(3,1)).unit(),(basis(3,0)+basis(3,1)).unit())
#    psi_target=tensor(basis(N,0),basis(3,1), (basis(3,0)-basis(3,1)).unit())
    tar = (tensor(basis(3,0),basis(3,0))+tensor(basis(3,0),basis(3,1))+tensor(basis(3,1),basis(3,0))-tensor(basis(3,1),basis(3,1))).unit()
    psi_target=tensor(basis(N,0),tar).unit()
    
#    psi_target=tensor((basis(3,0)-basis(3,1)).unit(),basis(3,1))
#    psi_target=tensor(basis(3,1), (basis(3,0)-basis(3,1)).unit())
#    psi_target=(tensor(basis(3,0),basis(3,0))+tensor(basis(3,0),basis(3,1))+tensor(basis(3,1),basis(3,0))-tensor(basis(3,1),basis(3,1))).unit()
    
    
    def getfid(m):
        global fx
        fx=m
        psi_total, tlist_total, expect_x,expect_y,expect_z,pltlist=EvolveCircuits(Operators,psi0)
        
        U0 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E_q[0])*tlist_total[-1])*basis(3,1)*basis(3,1).dag()
        U1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E_q[1])*tlist_total[-1])*basis(3,1)*basis(3,1).dag()
#        U2 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(E_q[2])*tlist_total[-1])*basis(3,1)*basis(3,1).dag()
        UT = tensor(qeye(3),U0,U1)
        fid=fidelity(UT*psi_total[-1]*psi_total[-1].dag()*UT.dag(), psi_target)
#        print(fid)
        return 1-fid
    
    w_q[1] = w
    print(w_q/2/np.pi)
    # xopt,val=fminbound(getfid,10,400, xtol=1e-07)
    xopt = 0
    val = 0
    print(xopt,1-val)
    return(xopt,1-val)
    
    
        
        
    
#    fid=fminbound(getfid,10,400, xtol=1e-07,disp=3)
#    print(fid)

if __name__=='__main__':#fminbound
    starttime=time.time()
    
        
    ## Approximation Type{ 0: Consider Cavity but no direct coupling between qubits; 1: No cavity but consider direct coupling between qubits; 2: two level system; 3: two level system without cavity}
    APRType=0
    
    ## Qubit number
    Num_Q=2
    
    ## Cavity Ladder Number
    N=3
    
    ## Savedata 
    SD=False
    
    ## PlotData
    PLT=False
    
    ## Test Gate mode
    TG = 0
    
    ## add dissipation
    Dis=False
    
    ## Drag
    DRAG=True
    
    ## Use monte-carlo solver
    MC=False
    
    if Dis==False:
        MC=False
    
    ## Cavity Frequency
    w_c= 4.0 * 2 * np.pi
    
    ## Qubits frequency
    w_q = np.array([ 5.0 , 5.2 , 5.4 , 5.6 , 5.8 , 4.0 , 4.2 , 4.4 , 4.6 , 4.8 ]) * 2 * np.pi
    
    ## Coupling Strength
    g = np.array([0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03]) * 2 * np.pi
    
    ## Qubits Anharmonicity
    w_qa=  np.array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]) * 2 * np.pi
    
    ## Base Temperature(K)
    n_th=0.01
    
    ## 1/T_1
    gamma = np.array([1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 ]) * 1e-3    
        
    ## 1/T_phi                
    gamma_phi = np.array([1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 , 1./10 ]) * 1e-3  
                
    ## Define Circuit
    Operators=[
        ['CZ1','I','I','I','I'],
#        ['I','I','CZ3','I','I'],
#        ['I','CZ3','I','I','I'],
#        ['CZ3','I','I','I','I'],
#         ['H','H'],
#         ['CZ1','I'],
#         ['I','H'],
#         ['I','Y-0.333333'],
#         ['I','H'],
#         ['CZ1','I'],
#         ['H','H'],
#         ['I','Y0.333333'],
        ]                    
    
    if APRType==0:
#        psi0=tensor(basis(N,0), basis(3,1), basis(3,0), (basis(3,0)+basis(3,1)).unit(), basis(3,0) ,basis(3,1))
        psi0=tensor(basis(N,0), (basis(3,0)+basis(3,1)).unit(), (basis(3,0)+basis(3,1)).unit(), (basis(3,0)+basis(3,1)).unit(), basis(3,0) ,basis(3,1))
#        psi0=tensor(basis(N,0), (basis(3,1)+basis(3,0)).unit(), (basis(3,1)-basis(3,0)).unit(), (basis(3,0)+basis(3,1)).unit(), basis(3,0) ,basis(3,1))
#        psi0=tensor(basis(N,0), basis(N,0), basis(N,1), basis(N,0), basis(3,0) )
    elif APRType==1:
#        psi0=tensor((basis(3,0)+basis(3,1)).unit(), (basis(3,0)+basis(3,1)).unit()) 
        psi0=tensor(basis(3,1), (basis(3,0)+basis(3,1)).unit()) 
#        psi0=tensor((basis(3,0)+basis(3,1)).unit(), (basis(3,1)).unit()) 
    elif APRType==2:
        psi0=tensor(basis(N,0), (basis(2,0)+basis(2,1)).unit(), (basis(2,0)+basis(2,1)).unit())#,basis(2,0) ,basis(2,0) ,(basis(2,0)+basis(2,1)).unit())
    elif APRType==3:
        psi0=tensor((basis(2,0)+basis(2,1)).unit(), (basis(2,0)+basis(2,1)).unit())    
            
    ### Effective Coulnp.ping Strength                    
    #g_eff=g[0]*g[1]/(w_q[0]-w_c)/(w_q[1]-w_c)*(w_q[0]-w_c+w_q[1]-w_c)
    #print(g_eff/(2*np.pi)*1000)
    
    
    options=Options()
#    options.atol=1e-11
#    options.rtol=1e-9
    options.first_step=0.01
    options.num_cpus=4
    options.nsteps=1e6
    options.gui='True'
    options.ntraj=1000
    options.rhs_reuse=False
    
    
    ## Create Basic Operator
    a, sm, smm, E_e, E_g, Sz, Sx, Sy=CreateBasicOperator(APRType)             
    
    ## Create Dissipation Term
    c_op_list=GenerateDissipation(Dis)
    
    ## Do Evolution
    if TG==False:
        psi_total, tlist_total, expect_x,expect_y,expect_z,pltlist=EvolveCircuits(Operators,psi0)

    ## Test Gate
    if TG==True:
        PLT=False
#        TestGate()
        optdelta = []
        optfid = []
        x = np.array([4.0+0.1*i for i in range(3)])*2*np.pi 
        for i in x:
            w_q[0] = i
            p = Pool(3)
            A = p.map(TestGate,x)
            p.close()
            p.join()
            delta = [k[0] for k in A]
            fid = [k[1] for k in A]
            optdelta.append(delta)
            optfid.append(fid)
            
        optdelta = np.array(optdelta)
        optfid = np.array(optfid)
        np.savetxt('optdelta.txt', optdelta)
        np.savetxt('optfid.txt', optfid)
        np.savez('data.npz', optdelta=optdelta, optfid=optfid)
    
    ## save data
    if SD:
        savedata(psi_total, pltlist, expect_x,expect_y,expect_z)

#    measurement0 = Measure(psi_total[-1],1)
#    measurement1 = Measure(psi_total[-1],1)

    
    finishtime=time.time()
    print( 'Time used: ', (finishtime-starttime), 's')

