import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.signal import find_peaks
import csv
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-theta', '--theta', default=0.0, type=float, help='theta')
parser.add_argument('-fs', '--fs', default=1001, type=int, help='samples per second')
parser.add_argument('-q', action='store_true', help='quiet mode')
parser.add_argument('--format', type=str, default='eps', help='file format (eps/jpg/png etc.)')
args = parser.parse_args()


plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['font.family'] = 'serif'

# Make Directory for Plots
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
directory = "data/plots"
create_directory(directory)

theta = args.theta
fs = args.fs

# κ='0.40'
κ_list=np.concatenate((np.arange(0.5,1.0,0.05),np.arange(1.0,5.01,0.05)))


Y_N=dict()
Y_D=dict()
Y_D_mf=dict()
for κ in κ_list:
    if not args.q:
        print(f'k = {κ:.2f}', end='\r')
    κ_string='{:.2f}'.format(κ)
    theta = '0.00'
    data = np.genfromtxt(f"data/N00N_theta_{theta}_k_{κ_string}_N_inf_matlab.csv", delimiter=",")
    num_points = data.shape[0]
    Y_N[κ]=np.average(data[num_points//3:,3])
    theta = '0.79'
    data = np.genfromtxt(f"data/N00N_theta_{theta}_k_{κ_string}_N_inf_matlab.csv", delimiter=",")
    num_points = data.shape[0]
    Y_D[κ]=np.average(data[num_points//3:,3])
    data = np.genfromtxt(f"data/N00N_theta_{theta}_k_{κ_string}_N_inf_matlab_mf.csv", delimiter=",")
    num_points = data.shape[0]
    Y_D_mf[κ]=np.average(data[num_points//3:,3])
    


plt.rcParams.update({'text.color': 'k',
                    'axes.labelcolor': 'k',
                    'axes.edgecolor': 'k',
                    'xtick.color':'k',
                    'ytick.color':'k',
                    'legend.facecolor':'w'})



c=0

fig, ax = plt.subplots(1, 1, figsize=(6,5))
ax.plot(1.0/κ_list, np.zeros(len(κ_list)), label=r'Mean Field $\theta=0$',linewidth=1.132,c='#cb3f10', linestyle=(0,(4,2.5)))
ax.plot(1.0/κ_list, np.abs(np.array(list(Y_D_mf.values()))), label=r'Mean Field $\theta=\pi/4$',linewidth=1.132,c='#D98378', linestyle=(0,(4,2.5)))
ax.plot(1.0/κ_list, np.abs(np.array(list(Y_N.values()))), label=r'Second Order $\theta=0$',linewidth=1.132,c='#1b485e', linestyle='-')
ax.plot(1.0/κ_list, np.abs(np.array(list(Y_D.values()))), label=r'Second Order $\theta=\pi/4$',linewidth=1.132,c='#568b87', linestyle='-')

ax.plot([1,1],[0.01,1.01],linewidth=1.132,c='#444444',linestyle=':')
ax.set_xlabel(r'$\Omega/\kappa$', fontsize=24)
ax.set_ylabel(r'$\left|\langle m_z\rangle_t\right|$', fontsize=24)
ax.legend(loc='upper right', fontsize=12)
ax.tick_params(which='both', direction='in', labelsize=16)
plt.savefig(f'data/plots/figure_2.{args.format}', transparent=True, bbox_inches='tight')
    