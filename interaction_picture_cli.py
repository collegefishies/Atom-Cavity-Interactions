import qutip
import argparse
import pandas as pd
import numpy as np
import os 
from interaction import simulation


#Argument Parser Section
#Argument Parser Section
#Argument Parser Section
#define our argument parser
parser = argparse.ArgumentParser(
	description='Atom-Cavity System Simulation'
)
parser.add_argument('-N', type=int, nargs='?', default=3, help='number of atoms in the system (default: 1)')
parser.add_argument('-d', type=float, required=True, help='detuning parameter, in units of the cavity linewidth')
parser.add_argument('-o', type=str, required=True, help='output file name (.csv)')
parser.add_argument('-T', type=float, nargs='?', default=60, help='simulation time in the scale units.')
args = parser.parse_args()

N       	= args.N
detuning	= args.d
filename	= args.o
T       	= args.T

if not filename.endswith('.csv'):
	filename = filename + '.csv'

#Simulation Section
#Simulation Section
#Simulation Section
print(f"Running simulation with {N} atoms and detuning {detuning}")
data = simulation(Natoms=N, detuning=detuning, T=T)
intracavity_photon_number = data['intracavity_photon_number']
print(f"Simulation returned {intracavity_photon_number}")




#Save/Append Data Section
#Save/Append Data Section
#Save/Append Data Section
print(f"Opening the pandas dataframe \"{filename}\" file so far...")
try:
	df = pd.read_csv(filename)
	print("Array opened: ", np.array(df['array'][0]))
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