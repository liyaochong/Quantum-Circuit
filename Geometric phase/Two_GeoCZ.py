from qutip import *
import numpy as np
import matplotlib.pyplot as plt 
from pylab import *
from time import clock
from scipy.special import *
from scipy.optimize import *

def Two_GeoCZ(ome = 8.77841):
    
    #==============================================================================
    w_c = 5.585  * 2 * np.pi  # cavity frequency
    w_q = np.array([ 5.849 , 5.870]) * 2 * np.pi
    g = np.array([0.0209 , 0.0198]) * 2 * np.pi
    eta_q = np.array([-0.245 , -0.244]) * 2 * np.pi
    N = 10              # number of cavity fock states
    n= 0
    #==============================================================================
    
    delta = 0.004 * 2 * np.pi
    omega = 0.001 * 2 * np.pi*np.sqrt(ome)
    
    
    #==============================================================================
    a = tensor(destroy(N),qeye(3),qeye(3))
    sm = np.array([tensor(qeye(N),destroy(3),qeye(3)) , tensor(qeye(N),qeye(3),destroy(3))])
    
    E_uc = np.array([tensor(qeye(N),basis(3,2)*basis(3,2).dag(),qeye(3)) , tensor(qeye(N),qeye(3), basis(3,2)*basis(3,2).dag())])
    #用以表征非简谐性的对角线最后一项(非计算能级)
    #E_uc1 = tensor(qeye(N),qeye(3), Qobj([[0,0],[0,1]]))
    
    E_e = np.array([tensor(qeye(N),basis(3,1)*basis(3,1).dag(),qeye(3)),tensor(qeye(N),qeye(3),basis(3,1)*basis(3,1).dag())])
    #激发态
    
    E_g = np.array([tensor(qeye(N),basis(3,0)*basis(3,0).dag(),qeye(3)) , tensor(qeye(N),qeye(3),basis(3,0)*basis(3,0).dag())])
    #基态
    
    sn = np.array([sm[0].dag()*sm[0] , sm[1].dag()*sm[1]])
    
    sx = np.array([sm[0].dag()+sm[0],sm[1].dag()+sm[1]]);
    sxm = np.array([tensor(qeye(N),Qobj([[0,1,0],[1,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(N),qeye(3),Qobj([[0,1,0],[1,0,0],[0,0,0]]))])
    
    
    sy = np.array([1j*(sm[0].dag()-sm[0]) , 1j*(sm[1].dag()-sm[1])]);
    sym = np.array([tensor(qeye(N),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]),qeye(3)) , tensor(qeye(N),qeye(3),Qobj([[0,-1j,0],[1j,0,0],[0,0,0]]))])
    
    sz = np.array([E_g[0] - E_e[0] , E_g[1] - E_e[1]])
    #==============================================================================
    '''Hamilton'''
    
    HCoupling = g[0]*(a+a.dag())*(sm[0]+sm[0].dag()) + g[1]*(a+a.dag())*(sm[1]+sm[1].dag())
    Hc = w_c * a.dag() * a 
    H_eta = eta_q[0] * E_uc[0] + eta_q[1] * E_uc[1]
    Hq = w_q[0]*sn[0] + w_q[1]*sn[1]
    H = Hq + H_eta + Hc + HCoupling
    
    Ee = H.eigenenergies()
    wd = Ee[1]+delta
    w = 'omega*np.sin(wd*t)' 
#    w = '(omega*np.exp(-(t-125)**2/2.0/10**2)*np.cos(wd*t))*(0<t<=250)'
#    w = '(omega*(np.exp(-(t-125)**2/2.0/10**2)*np.cos(wd*t)+(t-125)/2/10**2*np.exp(-(t-125)**2/2.0/10**2)*np.cos(wd*t+np.pi/2)))*(0<t<=250)'
    Hd = [2*(a+a.dag()),w]
    H = [H,Hd]
    args = {'wd':wd,'omega':omega}
    
    #==============================================================================
    '''dissipation'''
    Q = 35000
    kappa = w_c/Q
    kappa_phi = w_c/Q
    gamma = np.array([1.0/10 , 1.0/10]) *1e-3
    gamma_phi = np.array([1.0/10 , 1.0/10]) *1e-3
    n_th = 0.01
    cops = []
#    for ii in range(2):
#        
#        cops.append(np.sqrt(gamma[ii] * (1+n_th)) * sm[ii])
#        cops.append(np.sqrt(gamma[ii] * n_th) * sm[ii].dag())
#        cops.append(np.sqrt(gamma_phi[ii]) * sm[ii].dag()*sm[ii])
#    cops.append(np.sqrt(kappa * (1+n_th)) * a)
#    cops.append(np.sqrt(kappa * n_th) * a.dag())
#    cops.append(np.sqrt(kappa_phi) * a.dag()*a)
    #==============================================================================
    '''evolution'''
    
    options=Options()
    options.atol=1e-11
    options.rtol=1e-9
    options.first_step=0.01
    options.num_cpus= 4
    options.nsteps=1e6
    options.gui='True'
    options.ntraj=1000
    options.rhs_reuse=False
    
    
    psi0 = tensor(basis(N,n) , (basis(3,0)+basis(3,1)).unit() , (basis(3,0)+basis(3,1)).unit())
    #psi0 = tensor(basis(N,n) , (basis(3,0)).unit() , (basis(3,0)+basis(3,1)).unit())
    #psi0 = tensor(basis(N,n) , (basis(3,0)+basis(3,1)).unit() , (basis(3,0)).unit())
    tlist = np.linspace(0,250,251)
    result = mesolve(H,psi0,tlist,cops,[],args = args,options = options)
    
    #==============================================================================
    '''
    PLot
    '''
#    n_x0 = [] ; n_y0 = [] ; n_z0 = [];n_a = [];nax = [];nay = [];
#    n_x1 = [] ; n_y1 = [] ; n_z1 = [];
#    
#    for t in range(0,len(tlist)):
#        U = 'basis(N,0)*basis(N,0).dag()'
#        for i in range(1,n):
#            U = U+'+np.exp(1j*'+str(i)+'*wd*tlist[t])*basis(N,'+str(i)+')*basis(N,'+str(i)+').dag()'      
#        U = eval(U)
#        RF1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(Ee[2]-Ee[0])*tlist[t])*basis(3,1)*basis(3,1).dag()+np.exp(1j*(Ee[7]-Ee[0])*tlist[t])*basis(3,2)*basis(3,2).dag()
#        RF2 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(Ee[3]-Ee[0])*tlist[t])*basis(3,1)*basis(3,1).dag()+np.exp(1j*(Ee[9]-Ee[0])*tlist[t])*basis(3,2)*basis(3,2).dag()
#        U = tensor(U,RF1,RF2)
#        
#        opx0 = U.dag()*sx[0]*U
#        opy0 = U.dag()*sy[0]*U
#        opz0 = sz[0]
#        opx1 = U.dag()*sx[1]*U
#        opy1 = U.dag()*sy[1]*U
#        opz1 = sz[1]
#        opnx = U.dag()*(a+a.dag())*U/2
#        opny = U.dag()*(a-a.dag())/2/1J*U
#        n_x0.append(expect(opx0,result.states[t]))
#        n_y0.append(expect(opy0,result.states[t]))
#        n_z0.append(expect(opz0,result.states[t]))
#        n_x1.append(expect(opx1,result.states[t]))
#        n_y1.append(expect(opy1,result.states[t]))
#        n_z1.append(expect(opz1,result.states[t]))
#        n_a.append(expect(a.dag()*a,result.states[t]))
#        nax.append(expect(opnx,result.states[t]))
#        nay.append(expect(opny,result.states[t]))
#        
#        
#    fig, axes = plt.subplots(3, 1, figsize=(10,6))
#            
#    axes[0].plot(tlist, n_x0, label='X');axes[0].set_ylim([-1.05,1.05])
#    axes[1].plot(tlist, n_y0, label='Y');axes[1].set_ylim([-1.05,1.05])
#    axes[2].plot(tlist, n_z0, label='Z');axes[2].set_ylim([-1.05,1.05])
#    
#    fig, axes = plt.subplots(3, 1, figsize=(10,6))
#            
#    axes[0].plot(tlist, n_x1, label='X');axes[0].set_ylim([-1.05,1.05])
#    axes[1].plot(tlist, n_y1, label='Y');axes[1].set_ylim([-1.05,1.05])
#    axes[2].plot(tlist, n_z1, label='Z');axes[2].set_ylim([-1.05,1.05])
#    
#    
#    fig, axes = plt.subplots(2, 1, figsize=(10,6))
#    axes[0].plot(nax, nay);
#    axes[1].plot(tlist, n_a, label='N');
#    
#    
#    sphere = Bloch()
#    sphere.add_points([n_x0 , n_y0 , n_z0])
#    sphere.add_vectors([n_x0[-1],n_y0[-1],n_z0[-1]])
#    sphere.make_sphere() 
#    plt.show()
#    sphere = Bloch()
#    sphere.add_points([n_x1 , n_y1 , n_z1])
#    sphere.add_vectors([n_x1[-1],n_y1[-1],n_z1[-1]])
#    sphere.make_sphere() 
#    plt.show()
    #==============================================================================
    #==============================================================================
    '''
    fidelity
    '''
    U = 'basis(N,0)*basis(N,0).dag()'
    for i in range(1,n):
        U = U+'+np.exp(1j*'+str(i)+'*wd*tlist[-1])*basis(N,'+str(i)+')*basis(N,'+str(i)+').dag()'      
    U = eval(U)
    RF1 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(Ee[2]-Ee[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()+np.exp(1j*(Ee[7]-Ee[0])*tlist[-1])*basis(3,2)*basis(3,2).dag()
    RF2 = basis(3,0)*basis(3,0).dag()+np.exp(1j*(Ee[3]-Ee[0])*tlist[-1])*basis(3,1)*basis(3,1).dag()+np.exp(1j*(Ee[9]-Ee[0])*tlist[-1])*basis(3,2)*basis(3,2).dag()
    U = tensor(U,RF1,RF2)
    
    
    tar = (-tensor(basis(3,0) , basis(3,0))+tensor(basis(3,1) , basis(3,0))+tensor(basis(3,0) , basis(3,1))+tensor(basis(3,1) , basis(3,1))).unit()
    target = tensor(basis(N,n) , tar)
    #target = tensor(basis(N,n) , (basis(3,0)).unit() , (-basis(3,0)+basis(3,1)).unit())
    #target = tensor(basis(N,n) , (-basis(3,0)+basis(3,1)).unit() , (basis(3,0)).unit())
    fid = fidelity(U*result.states[-1]*result.states[-1].dag()*U.dag(),target)
    print('fidelity = ',fid)
    #==============================================================================
    return(1-fid)

if __name__ == '__main__':
    

    starttime = clock()
#    Two_GeoCZ(ome = 8.77841)
    fid=fminbound(Two_GeoCZ,8.0,100.0, xtol=1e-07,disp=3)
    
    endtime = clock()
    print('Time used :',endtime-starttime,'s' )

























