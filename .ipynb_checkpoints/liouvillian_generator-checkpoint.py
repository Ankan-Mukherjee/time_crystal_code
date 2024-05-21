#imports
import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from scipy.integrate import odeint
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.signal import find_peaks
import csv
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num_spins', nargs='+', default=[10], type=int)
parser.add_argument('-k', '--kappa', nargs='+', default=[0.4], type=float)
parser.add_argument('-ops', action='store_true', help='save operator matrices')
parser.add_argument('-q', action='store_true', help='quiet mode')
args = parser.parse_args()


# Make Directory for Liouvillian
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
directory = "data/liouvillian"
create_directory(directory)



# Actual Code
def liouvillian_generate(Nb=[10],k_list=[0.4],q=False,ops=False):
    ω=1
    for index,N in enumerate(Nb):
            for κ in k_list:
                if q==False:
                    print(f'N = {N}, k = {κ}', end='\r')
                lv_filename = 'data/liouvillian/lv_k_{0:.2f}_N_{1:03d}.csv'.format(κ,N)
                
                
                Jx = jmat(N/2, 'x')
                Jy = jmat(N/2, 'y')
                Jz = jmat(N/2, 'z')
                Jp = jmat(N/2, '+')
                Jm = jmat(N/2, '-')

                ops_dict = {'jx':Jx,'jy':Jy,'jz':Jz,'jp':Jp,'jm':Jm}

                if ops==True:
                    for op_name in ops_dict.keys():
                        op_filename = 'data/liouvillian/{0}_N_{1:03d}.csv'.format(op_name,N)
                        with open(op_filename, 'w') as csvfile: 
                            writer = csv.writer(csvfile)
                            writer.writerows(ops_dict[op_name].full())
                    
                H_strong = ω * Jx
                cops = [np.sqrt(κ/(N/2)) * Jm]
                liouv = liouvillian(H_strong, cops)
                data=liouv.full()
                with open(lv_filename, 'w') as csvfile: 
                    writer = csv.writer(csvfile)
                    writer.writerows(data)
                eigs_list=np.unique(np.round(np.abs(np.imag(liouv.eigenenergies())),4))[:20]
                eigs_filename = 'data/liouvillian/eigenvalues_k_{0:.2f}_N_{1:03d}.csv'.format(κ,N)
                data = eigs_list
                data = data.reshape(len(data),1)
                with open(eigs_filename, 'w') as csvfile: 
                    writer = csv.writer(csvfile)
                    writer.writerows(data)

# Running it
liouvillian_generate(Nb=args.num_spins,k_list=args.kappa,q=args.q,ops=args.ops)
if args.q==False:
    print('Process completed successfully!')