This is the code used to generate the plots in the paper

A. Mukherjee et. al. "Symmetries and Correlations in Continuous Time Crystals"

*************************************************************************************************************

Requirements:

The code requires MATLAB and python to generate the datasets and plot the same. We have tested the code on python=3.12.2 and MATLAB R2023b.

For python, the requirements are listed in requirements.txt and can be installed using

pip install -r requirements.txt

For MATLAB, install the following toolboxes:
Signal Processing Toolbox
Symbolic Math Toolbox

Make sure there is a folder data/ in the current directory.

*************************************************************************************************************

Make sure there is a folder called init in data/. This folder contains the initial states. Note that there are three files corresponding to each state number.

1. <state_number>_params.csv: Contains the parameters of the state. The state is described by a1*|theta=a2,phi=a3> + a4*|theta=a5,phi=a6> + a7*|theta=a8,phi=a9> + ......, where a1, a2, ... are the elements of the list.
2. <state_number>_mf.csv: Contains the mean field initial values in the order m_x, m_y, m_z.
3. <state_number>_so.csv: Contains the second order initial values in the order m_x, m_y, m_z, chi_xx, chi_xy, chi_xz, chi_yy, chi_yz, chi_zz.

Three initial states are provided as a reference (and for the plots in the paper). These are
0: All spins pointing up (|theta=0, phi=0>)
1: Cat state (equal superposition of all spins pointing up and all spins pointing down)
2: Equal superposition of |theta=pi/4,phi=0> and |theta=3pi/4,phi=0>

One can create files corresponding to any state in this fashion.

*************************************************************************************************************

Run the code from the folder with suitable arguments

python3 <code name>.py -arguments

Example:
python3 finite_evolution.py -n 50 200 -k 0.4 -state 1 -fs 1001

*************************************************************************************************************

The codes run as follows:

finite_evolution.py
options:
  -h, --help            show this help message and exit
  -n NUM_SPINS [NUM_SPINS ...], --num_spins NUM_SPINS [NUM_SPINS ...]
                        list of number of spins
  -k KAPPA, --kappa KAPPA
                        kappa, see the master equation
  -state STATE, --state STATE
                        state number for initial state
  -fs FS, --fs FS       samples per second
  -q                    quiet mode

Example:
python3 finite_evolution.py -n 50 200 -k 0.4 -state 1 -fs 1001 

*********************

liouvillian_generator.py
options:
  -h, --help            show this help message and exit
  -n NUM_SPINS [NUM_SPINS ...], --num_spins NUM_SPINS [NUM_SPINS ...]
  -k KAPPA [KAPPA ...], --kappa KAPPA [KAPPA ...]
  -ops                  save operator matrices
  -q                    quiet mode

Example:
python3 liouvillian_generator.py -n 60 -k 0.4 -ops


*********************

plot_evolution.py
options:
  -h, --help            show this help message and exit
  -n NUM_SPINS [NUM_SPINS ...], --num_spins NUM_SPINS [NUM_SPINS ...]
                        list of number of spins
  -k KAPPA, --kappa KAPPA
                        kappa, see the master equation
  -state STATE, --state STATE
                        state number for initial state
  -fs FS, --fs FS       samples per second
  -q                    quiet mode
  --format FORMAT       file format (eps/jpg/png etc.)

Example:
python3 plot_evolution.py --format eps


***********************

plot_fft.py
options:
  -h, --help            show this help message and exit
  -n NUM_SPINS, --num_spins NUM_SPINS
                        list of number of spins
  -k KAPPA, --kappa KAPPA
                        kappa, see the master equation
  -state STATE, --state STATE
                        state number for initial state
  -fs FS, --fs FS       samples per second
  -neig NUM_EIGS, --num_eigs NUM_EIGS
                        number of eigenvalues to truncate to
  -q                    quiet mode
  --format FORMAT       file format (eps/jpg/png etc.)

Example:
python3 plot_fft.py --format jpg

***********************

plot_phase_transition.py
options:
  -h, --help            show this help message and exit
  -fs FS, --fs FS       samples per second
  -state STATE [STATE ...], --state STATE [STATE ...]
                        state number for initial state
  -q                    quiet mode
  --format FORMAT       file format (eps/jpg/png etc.)

Example:
python3 plot_phase_transition.py --format jpg

************************

order1_ode_num.m
k_num = kappa as per master equation
state = state number
fs = sampling rate (set to 1001)

Set the initial values and run

************************

order2_ode_num.m
k_num = kappa as per master equation
state = state number
fs = sampling rate (set to 1001)

Set the initial values and run

************************

liouville_switch.m
k_num = kappa as per master equation
state = state number
N = number of spins (default at 60)
num_eigs = number of eigenvalues to use (sorted by abs(imaginary parts), provided real part is 0))
operator = operator (set as a list of operators multiplied from left to right)

************************

order_parameter.m
k_list = list of kappa to plot
state = state number
fs = sampling rate (set to 1001)

*************************************************************************************************************

Plots 1(a) and 1(b)
Run finite_evolution.py for N=[50,200], and kappa=0.4
Run order1_ode_num.m with k_num=0.4, state_number = 1
Run order2_ode_num.m with k_num=0.4, state_number = 1
Run plot_evolution.py

Code:
Execute python3 finite_evolution.py -n 50 200 -k 0.4 -state 1 -fs 1001
Run order1_ode_num.m with said parameters
Run order2_ode_num.m with said parameters
Execute python3 plot_evolution.py

*******************************************************

Plot 1(c)
Run order_parameter.m with the parameters provided
Run plot_phase_transition.py
 
Code:
Run order_parameter.m with state_number = 1
Run order_parameter.m with state_number = 2
Run order_parameter_mf.m with state_number = 1
Run order_parameter_mf.m with state_number = 2
Execute python3 plot_phase_transition.py 

*******************************************************

Plots 2(a) and 2(b)
Run liouvillian_generator.py for N=60, and kappa=0.4
Run liouville_switch.m with N=60, num_eigs=80, kappa=0.4, operator=["jz"], state_number=2, fs=10001
Run liouville_switch.m with N=60, num_eigs=80, kappa=0.4, operator=["jz","jz"], state_number=2, fs=10001
Run order2_ode_num.m with k_num=0.4, state_number = 2
Run order1_ode_num.m with k_num=0.4, state_number = 2
Run plot_fft.py

Code:
Execute python3 liouvillian_generator.py -n 60 -k 0.4 -ops 
Run liouville_switch.m with the said parameters
Run order1_ode_num.m with the said parameters
Run order2_ode_num.m with the said parameters
Execute python3 plot_fft.py -n 60 -k 0.4 -state 2 -fs 10001

*************************************************************************************************************

For all purposes, w_num = Omega is hard coded to 1 (the ratio of Omega to kappa determines the dynamics, we set kappa to be the independent parameter).
The Jupyter notebooks provide the code for second and third order cumulant equation derivation.

*************************************************************************************************************
