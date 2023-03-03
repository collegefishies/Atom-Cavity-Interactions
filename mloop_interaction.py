import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv
import numpy as np
import time
from interaction import simulation
import logging

l = logging

#initial parameters
GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
ETA = 4
G_SI = 2*np.pi*np.sqrt(4*ETA*GAMMA_SI*KAPPA_SI)
# OMEGA_SI = 2*np.pi*5e6
OMEGA_SI=10*KAPPA_SI
LAMBDA_SI = 0.1*KAPPA_SI
SCALE=KAPPA_SI


class ExampleInterface(mli.Interface):
	#initialization of the interface, including this method is optional
	def __init__(self):
		#you must include the super command to call the parent class, Interface, constructor
		super(ExampleInterface, self).__init__()

		#attributes can be added here
		#if you want to precalculate any variables this is the place to do it.

	#you must include the get_next_cost_dict method in your class
	#this method is called whenever M-LOOP wants to run an experiment
	def get_next_cost_dict(self, params_dict):

		#get parameters from the provided dictionary
		params = params_dict['params']


		#get the cost from your experiment
		x = params[0]
		y = params[1]

		cost = -np.sin(np.pi*x)*np.sin(np.pi*y) - 10*np.exp((-(x-0.7)**2 -(y-0.7)**2)/0.0001)
		#there is no uncertainty in our result
		uncer=0
		#the evaluation will always be a success
		bad = False
		#add a small time delay to mimic a real experiment
		# time.sleep(0.1)

		#the cost, uncertainty and bad boolean must all be returned as a dictionary
		#you can include other variables you want to record as well if you want
		cost_dict = {'cost': cost, 'uncer': uncer, 'bad': bad}
		return cost_dict

class SimulationInterface(mli.Interface):
	def __init__(self):
		super(SimulationInterface, self).__init__()

	def get_next_cost_dict(self, params_dict):
		params = params_dict['params']

		#scale em up, cus the function will descale em later.
		omega = params[0]*SCALE
		detuning = params[1]*SCALE
		ETA = params[2]
		# N = int(params[3])
		G_SI = 2*np.pi*np.sqrt(4*ETA*GAMMA_SI*KAPPA_SI)
		try:
			results = simulation(
				Natoms=3,
				detuning=detuning,
				Gamma_si=GAMMA_SI, 
				kappa_si=KAPPA_SI,
				eta=ETA, 
				g_si=G_SI,
				omega_si=omega,
				Lambda=LAMBDA_SI, #this parameter is the SHO driving strength
				scale=SCALE,
				spin_command='',
				T=60,
				SAVE_ARRAYS=True)

			# logging.info("Max Wineland Parameter: ", results['max_wineland_parameter'])
			# logging.info("Max Wineland Time: ", results['max_wineland_time'])
			cost = -results['max_wineland_parameter']
			uncer = 0
			bad = False
		except:
			cost = np.nan
			uncer = 0
			bad = False

		cost_dict = {'cost': cost, 'uncer': uncer, 'bad': bad}
		return cost_dict



def sample_main():
	#M-LOOP can be run with three commands

	#first create your interface
	interface = ExampleInterface()
	#next create the controller, provide it with your interface
	#and any options you want to set

	controller = mlc.create_controller(
			interface,
			controller_type='neural_net',
			max_num_runs=1000,
			num_params = 2,
			min_boundary = [0,0],
			max_boundary = [1,1],
		)

	#to run m-loop and find the optimal parameters just use the controller method optimize
	controller.optimize()

	#the results of the optimization will be saved to files and can also be accessed as attributes of the controller
	print("Best parameters found: ")
	print(controller.best_params)

def main():
	#M-LOOP	can be run with three commands
	logging.basicConfig(filename='mloop_sim.log', level=logging.DEBUG)
	#first create your interface
	interface = SimulationInterface()
	#next create the controller, provide it with your interface
	#and any options you want to set

	params = {
		'omega':   	[0,     	0,     	0.00001],
		'detuning':	[0,     	-20,   	0],
		'eta':     	[1.2276,	1.2276,	1.2277],
		# 'N' :    	[3,     	3,     	3.1],
	}
	param_names = list(params.keys())
	zipped = list(zip(*params.values()))
	first_params = list(zipped[0])
	min_boundary = list(zipped[1])
	max_boundary = list(zipped[2])
	print(param_names)
	print(first_params)
	print(min_boundary)
	print(max_boundary)
	l.info("Param range: ", params)
	controller = mlc.create_controller(
			interface,
			controller_type='neural_net',
			max_num_runs=1000,
			max_num_runs_without_better_params = 100,
			max_duration=36000,
			num_params = len(first_params),
			min_boundary = min_boundary,
			max_boundary = max_boundary,
			# trust_region = 1, #freedom to randomly sample the full range.
			first_params = first_params,
			training_type='random',
			param_names = param_names
		)

	#to run m-loop and find the optimal parameters just use the controller method optimize
	controller.optimize()

	#the results of the optimization will be saved to files and can also be accessed as attributes of the controller
	print("Best parameters found: ")
	print(controller.best_params)
if __name__ == "__main__":
	# sample_main()
	main()