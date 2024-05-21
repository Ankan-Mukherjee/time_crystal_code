import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.signal import find_peaks
import csv
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-fs', '--fs', default=1001, type=int, help='samples per second')
parser.add_argument('-state', '--state', default=[1,2], nargs='+', type=int, help='state number for initial state')
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

state_numbers = args.state
fs = args.fs

# κ='0.40'
κ_list=np.concatenate((np.arange(0.5,1.0,0.05),np.arange(1.0,5.01,0.05)))


Y=dict()
for state_number in state_numbers:
    Y[state_number]=dict()
    Y[state_number]['so']=dict()
    Y[state_number]['mf']=dict()
    for κ in κ_list:
        if not args.q:
            print(f'k = {κ:.2f}', end='\r')
        κ_string='{:.2f}'.format(κ)
        data = np.genfromtxt(f"data/evolution/state_{state_number:03d}_k_{κ_string}_N_inf_matlab.csv", delimiter=",")
        num_points = data.shape[0]
        Y[state_number]['so'][κ]=np.average(data[num_points//3:,3])
        data = np.genfromtxt(f"data/evolution/state_{state_number:03d}_k_{κ_string}_N_inf_matlab_mf.csv", delimiter=",")
        num_points = data.shape[0]
        Y[state_number]['mf'][κ]=np.average(data[num_points//3:,3])
    


plt.rcParams.update({'text.color': 'k',
                    'axes.labelcolor': 'k',
                    'axes.edgecolor': 'k',
                    'xtick.color':'k',
                    'ytick.color':'k',
                    'legend.facecolor':'w'})

color_list=[
    '#cb3f10',
    '#D98378',
    '#1b485e',
    '#568b87'
]
c=0

fig, ax = plt.subplots(1, 1, figsize=(6,5))
for state_number in state_numbers:
    ax.plot(1.0/κ_list, np.abs(np.array(list(Y[state_number]['mf'].values()))), label=f'Mean Field, state={state_number}',linewidth=1.132,c=color_list[c], linestyle=(0,(4,2.5)))
    c+=1
    ax.plot(1.0/κ_list, np.abs(np.array(list(Y[state_number]['so'].values()))), label=f'Second Order, state={state_number}',linewidth=1.132,c=color_list[c], linestyle='-')
    c+=1

ax.plot([1,1],[0.01,1.01],linewidth=1.132,c='#444444',linestyle=':')
ax.set_xlabel(r'$\Omega/\kappa$', fontsize=24)
ax.set_ylabel(r'$\left|\langle m_z\rangle_t\right|$', fontsize=24)
ax.legend(loc='upper right', fontsize=12)
ax.tick_params(which='both', direction='in', labelsize=16)
plt.savefig(f'data/plots/figure_2.{args.format}', transparent=True, bbox_inches='tight')
    