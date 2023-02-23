import subprocess
import numpy as np
from time import sleep

GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
ETA = 4
G_SI = 2*np.pi*ETA*np.sqrt(GAMMA_SI*KAPPA_SI)/4
OMEGA_SI = 10*KAPPA_SI
LAMBDA_SI = 0.1*KAPPA_SI

#parameters to test
# detuning_list = np.arange(-3*OMEGA_SI, 3*OMEGA_SI, KAPPA_SI/2)
detuning_list = [1,2,3]

number_of_params = len(detuning_list)
answer = input(f"This code has {number_of_params} parameters. Do you still want to run it? (y/n)")


if answer.lower() == 'y':
	processes = []
	for detuning in detuning_list:
		command = f'python interaction_picture_cli.py -N 1 -d {detuning} -T 1 -o batch/test{detuning}.pd'
		
		print(f"Running command '{command}'...")
		#call the command and capture the output
		process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		processes.append(process)
		# sleep(10)

	num_completed = 0
	while processes:
		for process in processes:
			if process.poll() is not None:
				processes.remove(process)
				num_completed += 1
				print(f"{num_completed} processes completed out of {number_of_params}.")
				break # process list has been modified, so restart loop
else:
	print("Code execution cancelled by the user.")