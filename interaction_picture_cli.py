import qutip
import argparse
import pandas as pd
import numpy as np
import os 
import time
from interaction import simulation
#Argument Parser Section
#Argument Parser Section
#Argument Parser Section
#define our argument parser
parser = argparse.ArgumentParser(
	description='Atom-Cavity System Simulation'
)

#initial parameters
GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
ETA = 4
G_SI = 2*np.pi*np.sqrt(4*ETA*GAMMA_SI*KAPPA_SI)
# OMEGA_SI = 2*np.pi*5e6
OMEGA_SI=10*KAPPA_SI
LAMBDA_SI = 0.1*KAPPA_SI
SCALE=KAPPA_SI

def gg(x,gamma,kappa):
	return 2*np.pi*np.sqrt(4*x*gamma*kappa)
parser.add_argument('-N', type=int, nargs='?', default=3, help='number of atoms in the system (default: 1)')
parser.add_argument('-d', type=float, required=True, help='detuning parameter, in units of the cavity linewidth')
parser.add_argument('-o', type=str, required=True, help='output file name (.csv)')
parser.add_argument('-T', type=float, nargs='?', default=60, help='simulation time in the scale units.')
parser.add_argument('--gamma', type=float, nargs='?', default=GAMMA_SI, help='atomic spontaneous emission rate')
parser.add_argument('--kappa', type=float, nargs='?', default=KAPPA_SI, help='cavity decay rate')
parser.add_argument('--eta', type=float, nargs='?', default=ETA, help='pump parameter')
parser.add_argument('-g', type=float, nargs='?', default=G_SI, help='atom-cavity coupling strength, the code uses the larger of either the eta calculated value or g.')
parser.add_argument('--omega', type=float, nargs='?', default=OMEGA_SI, help='cavity resonance frequency')
parser.add_argument('--driving_strength', type=float, nargs='?', default=LAMBDA_SI, help='SHO driving strength')
parser.add_argument('--scale', type=float, nargs='?', default=SCALE, help='simulation time scaling factor (default: 1)')
parser.add_argument('--spin_state_command', type=str, default=None, help='qutip command for defining the initial spin state, example command: "dd = qutip.tensor([(1/np.sqrt(2))*(u-1j*d),]*Natoms); spin_state = \'y\'"')
args = parser.parse_args()

N           	= args.N
detuning    	= args.d
filename    	= args.o
T           	= args.T
gamma       	= args.gamma
kappa       	= args.kappa
eta         	= args.eta
g           	= np.max([gg(eta, gamma, kappa), args.g])
omega       	= args.omega
Lambda      	= args.driving_strength
spin_command	= args.spin_state_command
scale       	= args.scale

if not filename.endswith('.csv'):
	filename = filename + '.csv'
print("filename: ", filename)
#Simulation Section
#Simulation Section
#Simulation Section
print(f"Running simulation with {N} atoms and detuning {detuning}")
data = simulation(
	Natoms=N, 
	detuning=detuning,
	Gamma_si=gamma,
	kappa_si=kappa,
	eta=eta,
	g_si=g,
	omega_si=omega,
	Lambda=Lambda,
	scale=scale,
	spin_command=spin_command,
	T=T)
# data = {"fake data": np.random.random(), "intracavity_photon_number": np.nan, "detuning": detuning}
intracavity_photon_number = data['intracavity_photon_number']
print(f"Simulation returned {intracavity_photon_number}")




#Save/Append Data Section
#Save/Append Data Section
#Save/Append Data Section
print(f"Opening the pandas dataframe \"{filename}\" file so far...")
try:
	df = pd.read_csv(filename)
	append = True
except FileNotFoundError:
	append = False

print("Making new data dataframe...")
df2 = pd.DataFrame(data, index = [0])

if append:
	print("Concatenating dataframes...")
	df = pd.concat([df, df2], ignore_index=True)
else:
	df = df2
print(f"Saving the result to a pandas dataframe...")
df.to_csv(filename, index=False)

print(f"Result saved to {filename} as a pandas dataframe in CSV format:")
print(df)