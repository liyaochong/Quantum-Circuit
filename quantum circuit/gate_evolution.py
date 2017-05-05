# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 16:20:57 2017

@author: Chen
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from gate import *
from initialsetting import *
from GenerateH import *
from evolutionplot import *
from dissipation import *

#==============================================================================
def gate_evolution(psi0 , Operator , setting = qusetting()):#体系的演化过程
    options=Options()
#    options.atol=1e-11
#    options.rtol=1e-9
    options.first_step=0.01
    options.num_cpus= 4
    options.nsteps=1e6
    options.gui='True'
    options.ntraj=1000
    options.rhs_reuse=False
#==============================================================================

#==============================================================================
    H , args , tlist = GenerateH(Operator , setting)
    c_op_list = dissipation(setting)
#    print(c_op_list)
#    print(H)
#    print(psi0)
    result = mesolve(H,psi0,tlist,c_op_list,[],args = args,options = options)
#    print(ptrace(result.states[-1],1))
    return result , tlist
#==============================================================================

