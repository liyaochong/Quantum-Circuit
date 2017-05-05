# -*- coding: utf-8 -*-

from qutip import *
from initialsetting import *

def Measure(state,target,setting = qusetting()):
    '''
    对某个state进行测量，测量比特为target
    返回值为测得为0,1的概率，以及测得0,1后体系分别的状态
    '''
    sm,E_uc,E_e,E_g,sn,sx,sxm,sy,sym,sz,En = initial(quset)[-11:]
    exp1 = expect(E_e[target],state)    #probability of excited state
    exp0 = expect(E_g[target],state)    #probability of ground state
    measurestate1 = E_e[target]*state*state.dag()*E_e[target].dag()/exp1 #state with 1
    measurestate0 = E_g[target]*state*state.dag()*E_g[target].dag()/exp0 #state with 0
                       
    return([[exp0,measurestate0],[exp1,measurestate1]])