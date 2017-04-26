# -*- coding: utf-8 -*-

from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *

w_c = 5.585  * 2 * np.pi  # cavity frequency
w_q = 5.849 * 2 * np.pi
g = 0.0198 * 2 * np.pi
eta_q = -0.244 * 2 * np.pi
N = 10              # number of cavity fock states
n= 0

delta = 0.004 * 2 * np.pi
omega = 0.001 * 2 * np.pi*np.sqrt(7)
wd = w_c+delta

a = tensor(destroy(N),qeye(3))
sm = tensor(qeye(N),destroy(3))

E_uc = tensor(qeye(N),basis(3,2)*basis(3,2).dag())
#用以表征非简谐性的对角线最后一项(非计算能级)
#E_uc1 = tensor(qeye(N),qeye(3), Qobj([[0,0],[0,1]]))

E_e = tensor(qeye(N),basis(3,1)*basis(3,1).dag())
#激发态

E_g = tensor(qeye(N),basis(3,0)*basis(3,0).dag())
#基态

sn = sm.dag()*sm

sx = sm.dag()+sm

sy = 1j*(sm.dag()-sm)

sz = E_g - E_e



HCoupling = g*(a+a.dag())*(sm+sm.dag())
Hc = w_c * a.dag() * a 
H_eta = eta_q * E_uc
Hq = w_q*sn
H = Hq + H_eta + Hc + HCoupling
Ee = H.eigenenergies()
print(Ee/2/np.pi)
Hd = [2*omega*(a+a.dag()),'np.cos(wd*t)']
H = [H,Hd]
args = {'wd':wd}

options=Options()
options.atol=1e-11
options.rtol=1e-9
options.first_step=0.01
options.num_cpus= 4
options.nsteps=1e6
options.gui='True'
options.ntraj=1000
options.rhs_reuse=False

psi0 = tensor(basis(N,n) , (basis(3,0)+basis(3,1)).unit())
#psi0 = tensor(basis(N,n) , basis(3,0))

tlist = np.linspace(0,250,251)

result = mesolve(H,psi0,tlist,[],[],args = args,options = options)
#result = mesolve(H,psi0,tlist,[],[],options = options)

n_x = [] ; n_y = [] ; n_z = [];n_a = []

for t in range(0,len(tlist)):
    U = 'basis(N,0)*basis(N,0).dag()'
    for i in range(1,n):
        U = U+'+np.exp(1j*'+str(i)+'*wd*tlist[t])*basis(N,'+str(i)+')*basis(N,'+str(i)+').dag()'      
    U = eval(U)
    RF = basis(3,0)*basis(3,0).dag()+np.exp(1j*(Ee[2])*tlist[t])*basis(3,1)*basis(3,1).dag()
    U = tensor(U,RF)
    
    opx = U.dag()*sx*U
    opy = U.dag()*sy*U
    opz = sz
    n_x.append(expect(opx,result.states[t]))
    n_y.append(expect(opy,result.states[t]))
    n_z.append(expect(opz,result.states[t]))
    n_a.append(expect(a.dag()*a,result.states[t]))
    
    
fig, axes = plt.subplots(4, 1, figsize=(10,6))
        
axes[0].plot(tlist, n_x, label='X');axes[0].set_ylim([-1.05,1.05])
axes[1].plot(tlist, n_y, label='Y');axes[1].set_ylim([-1.05,1.05])
axes[2].plot(tlist, n_z, label='Z');axes[2].set_ylim([-1.05,1.05])
axes[3].plot(tlist, n_a, label='N');
    
sphere = Bloch()
sphere.add_points([n_x , n_y , n_z])
sphere.add_vectors([n_x[-1],n_y[-1],n_z[-1]])
sphere.make_sphere() 
plt.show()













































