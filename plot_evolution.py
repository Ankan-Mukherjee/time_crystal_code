import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.signal import find_peaks
import csv
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num_spins', nargs='+', default=[50,200], type=int, help='list of number of spins')
parser.add_argument('-k', '--kappa', default=0.4, type=float, help='kappa, see the master equation')
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

color_list=[
    '#06592A',
    '#46AEA0',
    '#9A95B8',
    '#000000',
    '#56B4E9'
]

κ = args.kappa
Nb = args.num_spins
theta = args.theta
fs = args.fs

Y = dict()

for N in Nb:
    Y[N]=np.genfromtxt(f"data/N00N_theta_{theta:.2f}_k_{κ:.2f}_N_{N:03d}.csv", delimiter=",")
Y_inf_matlab = np.genfromtxt(f"data/N00N_theta_{theta:.2f}_k_{κ:.2f}_N_inf_matlab.csv", delimiter=",")


y_labels=[r'$\kappa t$',
    r'$m_x$',
    r'$m_y$',
    r'$m_z$',
    r'$\chi_{xx}$',
    r'$\chi_{xy}$',
    r'$\chi_{xz}$',
    r'$\chi_{yy}$',
    r'$\chi_{yz}$',
    r'$\chi_{zz}$',
]

plt.rcParams.update({'text.color': 'k',
                    'axes.labelcolor': 'k',
                    'axes.edgecolor': 'k',
                    'xtick.color':'k',
                    'ytick.color':'k',
                    'legend.facecolor':'w'})

for index in [3,9]:
    fig, ax = plt.subplots(1, 1, figsize=(6,5))
    c=0
    for N in Nb:
        Y_N = Y[N]
        ax.plot(Y_N[:,0], Y_N[:,index], label=fr'$N={N}$',linewidth=1.5,c=color_list[c])
        c+=1
    # ax.plot(Y_200[:,0], Y_200[:,index], label=r'$N=200$',linewidth=1.5,c=color_list[1])
    ax.plot(Y_inf_matlab[:,0], Y_inf_matlab[:,index], label=r'Second Order',linewidth=1.5,c=color_list[c])
    c=c+1
    ax.plot(Y_inf_matlab[:,0], np.zeros(len(Y_inf_matlab[:,index])), label=r'Mean Field',linewidth=1.5,linestyle='--',c=color_list[c])

    ax.set_xlabel(r'$\kappa t$', fontsize=24)
    ax.set_ylabel(y_labels[index], fontsize=24)
    ax.legend(loc='upper right', fontsize=12)
    ax.tick_params(which='both', direction='in', labelsize=16)
    ax.set_xlim([15,30])
    plt.savefig(f'data/plots/figure_1_{index}.{args.format}', transparent=True, bbox_inches='tight')