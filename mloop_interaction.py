import mloop.interfaces as mli
import mloop.controllers as mlc
import mloop.visualizations as mlv
import numpy as np
import time
from interaction import simulation
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


def sample_main():
	#M-LOOP can be run with three commands

	#first create your interface
	interface = ExampleInterface()
	#next create the controller, provide it with your interface
	#and any options you want to set

	params = {
		'omega': []
	}
	controller = mlc.create_controller(
			interface,
			controller_type='neural_net',
			max_num_runs=1000,
			num_params = 2,
			min_boundary = [0,0],
			max_boundary = [1,1]
		)

	#to run m-loop and find the optimal parameters just use the controller method optimize
	controller.optimize()

	#the results of the optimization will be saved to files and can also be accessed as attributes of the controller
	print("Best parameters found: ")
	print(controller.best_params)

if __name__ == "__main__":
	sample_main()