import subprocess
import numpy as np
from time import sleep, time
from multiprocessing import Semaphore
import os
import pandas as pd
import argparse

#Argument Parser
parser = argparse.ArgumentParser(
	description='Atom-Cavity System Simulation'
)
parser.add_argument('-p', type=int, nargs='?', default=10, help='number of atoms in the system (default: 1)')
parser.add_argument('--onlycompile', action='store_true', help='a flag to only compile old csvs in cased simulation crashed earlier')

args = parser.parse_args()

max_processes = args.p



N = 3
Tmax = 60
folder = 'batch'
fname_root = f'{folder}/data'
# GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
# ETA = 4
# G_SI = 2*np.pi*ETA*np.sqrt(GAMMA_SI*KAPPA_SI)/4
OMEGA_SI = 10*KAPPA_SI
# LAMBDA_SI = 0.1*KAPPA_SI

#parameters to test
detuning_list = np.arange(-2.1*OMEGA_SI, +OMEGA_SI, KAPPA_SI/5)

number_of_params = len(detuning_list)
if not args.onlycompile:
	answer = input(f"This code has {number_of_params} parameters. Do you still want to run it? (y/n)")
else:
	answer = 'n'

if answer.lower() == 'y' and not args.onlycompile:
	# create a semaphore to limit the number of concurrent processes
	semaphore = Semaphore(max_processes)
	# list to store the processes
	processes = []
	start_times = []

	process_id = 0
	remaining_processes = len(detuning_list)
	for detuning in detuning_list:
		semaphore.acquire()  # acquire a permit from the semaphore
		command = f'python interaction_picture_cli.py -N {N} -d {detuning} -T {Tmax} -o {fname_root}{process_id}.csv'
		process_id += 1
		print(f"Running command '{command}'...")
		# start the process and add it to the list of processes
		process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		processes.append(process)
		start_times.append(time())

		while len(processes) == max_processes:
			for i in range(len(processes)):
				p = processes[i]
				if p.poll() is not None:
					return_code = p.poll()
					stdout, stderr = process.communicate()
					processes.remove(processes[i])
					semaphore.release()  # release the permit to the semaphore
					print(f"Process finished. {remaining_processes} processes remaining.")
					remaining_processes -= 1
					print(f"Process took {time()-start_times[i]:.2f}s to complete.")
					del start_times[i]
					if p.poll() !=0:
						print(f"Std Err: {stderr.decode('utf-8')}")
						# raise Exception("Subprocess returned error.")
					break; #restart loop so it checks through all the processes again.

	# check if any of the processes have finished and remove them from the list
	while len(processes) > 0:
		for p in processes:
			if p.poll() is not None:
				stdout, stderr = process.communicate()
				processes.remove(p)
				semaphore.release()  # release the permit to the semaphore
				print(f"Process finished. {len(processes)} processes remaining.")
				print(f"Std Err: {str(stderr)}")
if answer.lower() == 'y' or args.onlycompile:
	#compile all the csv's into a single dataframe.
	dataframes = []
	try:
		combined_df = pd.read_csv(fname_root + '.csv', index_col=0, header=0)
		combined_df.drop(combined_df.columns[0])
		dataframes.append(combined_df)
	except:
		pass
	filenames = os.listdir(folder)
	for fname in filenames:
		if fname.endswith('.csv'):
			filepath = os.path.join(folder, fname)
			df = pd.read_csv(filepath)
			dataframes.append(df)
			os.remove(filepath)

	combined_df = pd.concat(dataframes, ignore_index=True)
	combined_df.to_csv(fname_root + '.csv', index=False)

else:
	print("Code execution cancelled by the user.")
