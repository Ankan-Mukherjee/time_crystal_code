import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.signal import find_peaks
import csv
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num_spins', default=60, type=int, help='list of number of spins')
parser.add_argument('-k', '--kappa', default=0.4, type=float, help='kappa, see the master equation')
parser.add_argument('-theta', '--theta', default=0.7854, type=float, help='theta')
parser.add_argument('-fs', '--fs', default=1001, type=int, help='samples per second')
parser.add_argument('-neig', '--num_eigs', default=80, type=int, help='number of eigenvalues to truncate to')
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

color_list = [
    '#06592A',
    '#46AEA0',
    '#56B4E9',
    '#9A95B8']

κ = args.kappa
N = args.num_spins
theta = args.theta
fs = args.fs
numeigs = args.num_eigs

Y_inf_matlab = np.genfromtxt(f"data/N00N_theta_{theta:.2f}_k_{κ:.2f}_N_inf_matlab.csv", delimiter=",")
Y_inf_matlab_mf = np.genfromtxt(f"data/N00N_theta_{theta:.2f}_k_{κ:.2f}_N_inf_matlab_mf.csv", delimiter=",")

overlaps1 = np.genfromtxt(f'data/liouvillian/overlap_jz_N00N_theta_{theta:.2f}_numeigs_{numeigs:02d}_N_{N:03d}.csv', delimiter=",")
overlaps2 = np.genfromtxt(f'data/liouvillian/overlap_jzjz_N00N_theta_{theta:.2f}_numeigs_{numeigs:02d}_N_{N:03d}.csv', delimiter=",")

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

def fft(Y, index=3):    
    label = r'$N\rightarrow\infty$'
    peak_ind = find_peaks(Y[:,index])[0]
    first_ind = peak_ind[2]
    last_ind = peak_ind[-2 if len(peak_ind)%2==0 else -1]-1
    Y=Y[first_ind:last_ind,:]
    N = Y.shape[0]
    amplitude = np.fft.fft(Y[:,index])
    freqs = np.fft.fftfreq(len(Y[:,index]), 1/fs)
    fft_x = np.fft.fftshift(2*np.pi*freqs)
    fft_y = np.fft.fftshift(np.abs(amplitude)/np.amax(np.abs(amplitude)))
    return (fft_x, fft_y)

fft1_inf = fft(Y_inf_matlab,3)
fft1_inf_mf = fft(Y_inf_matlab_mf,3)
fft2_inf = fft(Y_inf_matlab,9)

from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)



# First Order
fig, ax = plt.subplots(1, 1, figsize=(10,7))
ax.plot([0],[0],c='w',label=fr'$\Omega/\kappa={1/0.4}$')
ax.scatter(overlaps1[:,0],overlaps1[:,1]*2.4, s=80, edgecolors=color_list[2], marker='d',facecolors='none',label=r'$\hat{\mathcal{L}}$ Eigenvalues')
ax.plot(fft1_inf_mf[0], fft1_inf_mf[1], label=r'Mean Field',linewidth=1,alpha=0.9,c=color_list[3], linestyle='-')
ax.plot(fft1_inf[0], fft1_inf[1], label=r'Second Order',linewidth=1,alpha=0.9,c=color_list[0], linestyle='-')
ax.set_xlim([0,4])
ax.set_xlabel(r'Frequency', fontsize=22)
ax.set_ylabel(r'$\mathcal{F}\ [m_z]$', fontsize=32)
ax.legend(fontsize=16, bbox_to_anchor=(0.70,0.04,0.3,0.4), bbox_transform=ax.transAxes)
ax.tick_params(which='both', direction='in', labelsize=16)
ax.tick_params(which='both', direction='in')
axins1 = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=(0.53,0.53,0.45,0.45), bbox_transform=ax.transAxes)
axins1.plot(Y_inf_matlab_mf[:,0], Y_inf_matlab_mf[:,3], label=r'Mean Field',linewidth=1.0,alpha=0.9,c=color_list[3])
axins1.plot(Y_inf_matlab[:,0], Y_inf_matlab[:,3], label=r'Second Order',linewidth=1.0,alpha=0.9,c=color_list[0])
axins1.text(0.5,-0.08,fr'time', ha='center', va='top', transform=plt.gca().transAxes, fontsize=16)
axins1.text(-0.08,0.5,fr'$m_z$', rotation='vertical', ha='right', va='center', transform=plt.gca().transAxes, fontsize=20)
axins1.set_xlim([-0.1,0.4*np.max(Y_inf_matlab[:,0])])
axins1.set_ylim([0.8*np.min(Y_inf_matlab[:,3]),1.1*np.max(Y_inf_matlab[:,3])])
axins1.legend(loc='upper center', ncol=2)
axins1.set_yticks([-0.5,0,0.5])
axins1.tick_params(direction='in', which='both')

axins2 = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=(0.05,0.65,0.15,0.3), bbox_transform=ax.transAxes)
for spine in axins2.spines.values():
    spine.set_edgecolor('0.5')
mark_inset(ax, axins2, loc1=2, loc2=4, ec='0.5')
axins2.scatter(overlaps1[:,0],overlaps1[:,1]*2.4, s=100, edgecolors=color_list[2], marker='d',facecolors='none',label='Eigenvalues')
axins2.plot(fft1_inf_mf[0], fft1_inf_mf[1], label=r'Mean Field',linewidth=2,alpha=0.9,c=color_list[3], linestyle='-')
axins2.plot(fft1_inf[0], fft1_inf[1], label=r'Second Order',linewidth=2,alpha=0.9,c=color_list[0], linestyle='-')
axins2.set_xlim([0.87,1.02])
axins2.set_ylim([0.9,1.05])
axins2.set_xticks([])
axins2.set_yticks([])

axins3 = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=(0.28,0.35,0.15,0.35), bbox_transform=ax.transAxes)
for spine in axins3.spines.values():
    spine.set_edgecolor('0.5')
mark_inset(ax, axins3, loc1=3, loc2=1, ec='0.5')
axins3.scatter(overlaps1[:,0],overlaps1[:,1]*2.4, s=100, edgecolors=color_list[2], marker='d',facecolors='none',label='Eigenvalues')
axins3.plot(fft1_inf_mf[0], fft1_inf_mf[1], label=r'Mean Field',linewidth=2,alpha=0.9,c=color_list[3], linestyle='-')
axins3.plot(fft1_inf[0], fft1_inf[1], label=r'Second Order',linewidth=2,alpha=0.9,c=color_list[0], linestyle='-')
# axins1.set_xlabel('y')
axins3.set_xlim([1.75,2.05])
axins3.set_ylim([-0.02,0.28])
axins3.set_xticks([])
axins3.set_yticks([])

plt.savefig(f'data/plots/figure_3_1.{args.format}', transparent=True, bbox_inches='tight')





## Second Order
fig, ax = plt.subplots(1, 1, figsize=(10,7))
ax.plot([0],[0],c='w',label=fr'$\Omega/\kappa={1/0.4}$')
ax.scatter(overlaps2[:,0],overlaps2[:,1]*2.5, s=80, edgecolors=color_list[2], marker='d',facecolors='none',label=r'$\hat{\mathcal{L}}$ Eigenvalues')
ax.plot(fft2_inf[0], fft2_inf[1], label=r'Second Order',linewidth=1,alpha=0.9,c=color_list[0], linestyle='-')
ax.set_xlim([0,4])
ax.set_xlabel(r'Frequency', fontsize=22)
ax.set_ylabel(r'$\mathcal{F}\ [\chi_{zz}]$', fontsize=32)
ax.legend(fontsize=16, bbox_to_anchor=(0.05,0.6,0.3,0.4), bbox_transform=ax.transAxes)
ax.tick_params(which='both', direction='in', labelsize=16)
ax.tick_params(which='both', direction='in')

axins1 = inset_axes(ax, width='100%', height='100%', bbox_to_anchor=(0.53,0.53,0.45,0.45), bbox_transform=ax.transAxes)
axins1.plot(Y_inf_matlab[:,0], Y_inf_matlab[:,9], label=r'Second Order',linewidth=1.0,alpha=0.9,c=color_list[0])
axins1.text(0.5,-0.08,fr'time', ha='center', va='top', transform=plt.gca().transAxes, fontsize=16)
axins1.text(-0.05,0.5,fr'$\chi_{{zz}}$', rotation='vertical', ha='right', va='center', transform=plt.gca().transAxes, fontsize=20)
axins1.set_xlim([-0.1,0.4*np.max(Y_inf_matlab[:,0])])
axins1.set_ylim([-0.01,1.3*np.max(Y_inf_matlab[:,9])])
axins1.legend(loc='upper center', ncol=2)
axins1.set_yticks([0,0.5])
axins1.tick_params(direction='in', which='both')

plt.savefig(f'data/plots/figure_3_2.{args.format}', transparent=True, bbox_inches='tight')

