#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 21:36:26 2017

@author: chen
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from gate_evolution import *


from time import clock
starttime=clock()

gateset = gatesetting()
#Operator = ['CZ1','I']
Operator = ['I600'  ,  'I600' ]
#Operator = ['iswap1'  ,  'I' ]
#psi0 = tensor(basis(N,n) , basis(3,1) ,  (basis(level,0)+basis(level,1)).unit())
psi0 = tensor(basis(N,n) , basis(3,1) ,  (basis(level,0)+basis(level,1)).unit())
#psi0 = tensor(basis(N,n) , basis(level,1) ,  basis(level,0))
#gateset.omega = 0.033348628532133037
#gateset.iswap_deltat = 20
result , tlist = gate_evolution(psi0 , Operator , gateset = gateset)
#evolutionplot(0 , result , tlist )
evolutionplot(1 , result , tlist )


finishtime=clock()
print ('Time used: ', (finishtime-starttime), 's')