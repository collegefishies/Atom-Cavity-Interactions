import subprocess
import numpy as np
from time import sleep
from multiprocessing import Semaphore

N = 1
Tmax = 0.1
fname_root = 'batch/data'
GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
ETA = 4
G_SI = 2*np.pi*ETA*np.sqrt(GAMMA_SI*KAPPA_SI)/4
OMEGA_SI = 10*KAPPA_SI
LAMBDA_SI = 0.1*KAPPA_SI

#parameters to test
# detuning_list = np.arange(-3*OMEGA_SI, 3*OMEGA_SI, KAPPA_SI/2)
detuning_list = np.arange(21)

number_of_params = len(detuning_list)
answer = input(f"This code has {number_of_params} parameters. Do you still want to run it? (y/n)")

if answer.lower() == 'y':
	# create a semaphore to limit the number of concurrent processes
	max_processes = 10
	semaphore = Semaphore(max_processes)
	# list to store the processes
	processes = []

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

		while len(processes) == max_processes:
			for p in processes:
				if p.poll() is not None:
					processes.remove(p)
					semaphore.release()  # release the permit to the semaphore
					remaining_processes -= 1
					print(f"Process finished. {remaining_processes} processes remaining.")

	# check if any of the processes have finished and remove them from the list
	while len(processes) > 0:
		for p in processes:
			if p.poll() is not None:
				processes.remove(p)
				semaphore.release()  # release the permit to the semaphore
				print(f"Process finished. {len(processes)} processes remaining.")


else:
	print("Code execution cancelled by the user.")
