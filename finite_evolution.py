# Imports
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
from scipy.fft import fft, fftfreq, fftshift
import csv
from IPython.display import display, Latex
import scipy as sc
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num_spins', nargs='+', default=[50,200], type=int, help='list of number of spins')
parser.add_argument('-k', '--kappa', default=0.4, type=float, help='kappa, see the master equation')
parser.add_argument('-theta', '--theta', default=0.7854, type=float, help='theta')
parser.add_argument('-phi', '--phi', default=0.0, type=float, help='phi')
parser.add_argument('-fs', '--fs', default=101, type=int, help='samples per second')
parser.add_argument('-q', action='store_true', help='quiet mode')
args = parser.parse_args()

def state_evolution(N=36,Ω=1,κ=0.2,theta=0,phi=0,t=None,q=False):
    j = N/2
    Jx = jmat(N/2,"x")
    Jy = jmat(N/2,"y")
    Jz = jmat(N/2,"z")
    Jm = jmat(N/2,"-")
    Jp = jmat(N/2,"+")
    state = (1/np.sqrt(2))*(spin_coherent(j, np.pi-theta, phi)+spin_coherent(j, theta, phi))
    ρ0=state*state.dag()
    ρ0=ρ0/ρ0.tr()
    if t is None:
        t=np.linspace(0,1000,80000)
    H = Ω*Jx
    c_ops=[np.sqrt(κ/j)*Jm]
    L=liouvillian(H,c_ops)
    if q == True:
        sol=mesolve(H, ρ0, t, c_ops, [Jx/j, Jy/j, (Jz/j), (Jx*Jx)/(j**2), (Jy*Jy)/(j**2) , (Jz*Jz)/(j**2), (Jx*Jy)/(j**2), (Jy*Jz)/(j**2), (Jz*Jx)/(j**2), (Jy*Jx)/(j**2), (Jz*Jy)/(j**2), (Jx*Jz)/(j**2)])
    else:
        sol=mesolve(H, ρ0, t, c_ops, [Jx/j, Jy/j, (Jz/j), (Jx*Jx)/(j**2), (Jy*Jy)/(j**2) , (Jz*Jz)/(j**2), (Jx*Jy)/(j**2), (Jy*Jz)/(j**2), (Jz*Jx)/(j**2), (Jy*Jx)/(j**2), (Jz*Jy)/(j**2), (Jx*Jz)/(j**2)],progress_bar=True)
    return L,sol

def evolve(N_list=[10],k=0.4,fs=101,theta=np.pi/4,phi=0,q=False):
    Ω = 1;
    κ = k;
    t_end = round(100/κ);
    tSpan = np.linspace(0, t_end, fs*t_end);
    
    for N in N_list:
        filename = 'data/N00N_theta_{0:.2f}_k_{1:.2f}_N_{2:03d}.csv'.format(theta,κ,N)
        if q==False:
            print(f'N = {N}',end='\r')
        if os.path.exists(filename):
                continue
        L,sol = state_evolution(N=N,Ω=Ω,κ=κ,theta=theta,t=tSpan,phi=phi,q=q)
    
        ## Expectations
        expectations = dict()
    
        expectations['m_x'] = sol.expect[0]
        expectations['m_y'] = sol.expect[1]
        expectations['m_z'] = sol.expect[2]
    
        expectations['χ_xx'] = np.real(0.5*(sol.expect[3] + sol.expect[3])) - sol.expect[0]*sol.expect[0]
        expectations['χ_xy'] = np.real(0.5*(sol.expect[6] + sol.expect[9])) - sol.expect[0]*sol.expect[1]
        expectations['χ_xz'] = np.real(0.5*(sol.expect[8] + sol.expect[11])) - sol.expect[0]*sol.expect[2]
        expectations['χ_yy'] = np.real(0.5*(sol.expect[4] + sol.expect[4])) - sol.expect[1]*sol.expect[1]
        expectations['χ_yz'] = np.real(0.5*(sol.expect[7] + sol.expect[10])) - sol.expect[1]*sol.expect[2]
        expectations['χ_zz'] = np.real(0.5*(sol.expect[5] + sol.expect[5])) - sol.expect[2]*sol.expect[2]
        
        
        ## Saving data
        data=np.vstack((κ*tSpan,expectations['m_x'],expectations['m_y'],expectations['m_z'],
               expectations['χ_xx'],expectations['χ_xy'],expectations['χ_xz'],
              expectations['χ_yy'],expectations['χ_yz'],expectations['χ_zz'])).T
        
        with open(filename, 'w') as csvfile: 
            writer = csv.writer(csvfile)
            writer.writerows(data)

# Running the code
evolve(N_list=args.num_spins,k=args.kappa,fs=args.fs,theta=args.theta,phi=args.phi,q=args.q)
if args.q==False:
    print('Process completed successfully')
