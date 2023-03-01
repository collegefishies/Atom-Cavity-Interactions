import qutip
from qutip import mesolve
from qutip.ui.progressbar import TextProgressBar, BaseProgressBar
import numpy as np

SAVE_ARRAYS = False

GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
ETA = 20
G_SI = 2*np.pi*np.sqrt(4*ETA*GAMMA_SI*KAPPA_SI)
# OMEGA_SI = 2*np.pi*5e6
OMEGA_SI=10*KAPPA_SI
LAMBDA_SI = 0.1*KAPPA_SI
SCALE=KAPPA_SI
#print a coupling constant in engineering notation
print(f"Atom Coupling Constant: {G_SI:.2e} rad/s")

def make_natoms_operator(sx, id_atom, Natoms):
	'''
		Makes the natoms version of the single particle number operator, using
		the identity, and corresponding operator.
	'''
	array_of_identities = [id_atom,]*Natoms
	Sx = 0*qutip.tensor(*array_of_identities)
	for i in range(Natoms):
		temp = [id_atom,]*Natoms
		temp[i] = sx

		Sx = Sx + qutip.tensor(*temp)
	return Sx

def pretty_print_dict(dictionary):
	for key, value in dictionary.items():
		print(f"{key}: {value}")

def simulation(Natoms=3,detuning=0,
	Gamma_si=GAMMA_SI, 
	kappa_si=KAPPA_SI,
	eta=ETA, 
	g_si=G_SI,
	omega_si=OMEGA_SI,
	Lambda=LAMBDA_SI, #this parameter is the SHO driving strength
	scale=SCALE,
	spin_command='',
	T=60):
	#create the atom states
	u = qutip.basis(3,0)
	d = qutip.basis(3,1)
	e = qutip.basis(3,2)
	id_atom = qutip.qeye(3)

	#create the cavity fock states
	zero = qutip.basis(3,0)
	one = qutip.basis(3,1)
	two = qutip.basis(3,2)
	id_cav = qutip.qeye(3)

	#define the spin operators using the spin states
	sz = (u*u.dag() - d*d.dag())/2
	sx = (u*d.dag() + d*u.dag())/2
	sy = (u*d.dag() - d*u.dag())/(2*1j)

	#define the many atom spin operators
	sx = make_natoms_operator(sx, id_atom, Natoms)
	sy = make_natoms_operator(sy, id_atom, Natoms)
	sz = make_natoms_operator(sz, id_atom, Natoms)
	sz = qutip.tensor(sz, id_cav)
	sx = qutip.tensor(sx, id_cav)
	sy = qutip.tensor(sy, id_cav)
	id_atoms = qutip.tensor(*([id_atom,]*Natoms))

	#make the spin up to e and spin down to e excitation operators.
	nu = u * e.dag()
	mu = d * e.dag()
	nu = make_natoms_operator(nu, id_atom, Natoms)
	mu = make_natoms_operator(mu, id_atom, Natoms)
	mu = qutip.tensor(mu, id_cav)
	nu = qutip.tensor(nu, id_cav)

	#define the annihilation operator
	a = qutip.destroy(3)
	a = qutip.tensor(id_atoms, a)

	#rescaling the parameters by the scale parameter
	print(f"Rescaling parameters by the scale parameter {scale}")
	g = g_si/scale
	g = g_si/scale
	kappa = kappa_si/scale
	Gamma = Gamma_si/scale
	omega = omega_si/scale
	detuning = detuning/scale
	Lambda = Lambda/scale
	parameters_to_save = {
		'N': Natoms,
		'detuning': detuning,
		'g': g,
		'kappa': kappa,
		'Gamma': Gamma,
		'omega': omega,
		'Lambda': Lambda,
		'scale' : scale,
		'Tmax' : T,
	}
	print("These are the parameters to use...")
	pretty_print_dict(parameters_to_save)

	#define the Hamiltonian in the interaction picture
	def H(t,args):
		try:
			w = args['w']
		except:
			w = 0
		Hpart = g*(a*nu.dag()*np.cos(omega*t)- 1j*a*mu.dag()*np.sin(omega*t)) + Lambda*(np.exp(1j*w*t)*a)

		Ht = Hpart + Hpart.dag()
		return Ht

	#define the initial state
	# dd = qutip.tensor([(1/np.sqrt(2))*(u-d),]*Natoms); spin_state = 'x' #x state
	dd = qutip.tensor([(1/np.sqrt(2))*(u-1j*d),]*Natoms); spin_state = 'y' #y state
	# dd = qutip.tensor([(d),]*Natoms); spin_state='z' #z state
	if len(spin_command) > 3:
		print("Evaluating the spin command...")
		dd = eval(eval(spin_command))
		print("Spin Command evaluated!")
		# raise Exception(str(dd))
	psi0 = qutip.tensor(dd, zero)
	rho0 = psi0*psi0.dag()
	print(f"Spin State is initially {spin_state}")

	#define measurement operators
	id_spin = u*u.dag() + d*d.dag()
	id_spin = qutip.tensor(*[id_spin,]*Natoms, id_cav)

	#define the lindblad operators
	cavity_decay = a*np.sqrt(kappa)
	atom_decay = nu*np.sqrt(Gamma)

	#set up the simulation
	print("Setting up the simulation...")
	# opts=None
	opts = qutip.solver.Options(nsteps=120000)
	print("Executing simulation...")
	solution = mesolve(
		H=H,
		rho0=rho0,
		tlist=np.linspace(0, T*2*np.pi/kappa/10, 100),
		c_ops=[cavity_decay,atom_decay],
		e_ops=None,
		args={'w': detuning},
		options=opts,
		progress_bar=TextProgressBar()
	)
	print("Simulation Done!")

	#calculating expectation values
	states = solution.states
	t = solution.times
	e_ops = [sx, sy, sz, sx*sx, sy*sy, sz*sz, id_spin, a.dag()*a, sy*sz+sz*sy, sz*sx + sx*sz, sx*sy + sy*sz]
	e_vals = [np.array([qutip.expect(e_op, state) for state in states]) for e_op in e_ops]
	sx, sy, sz, sxx, syy, szz, norm, photons, syz, szx, sxy = e_vals
	entropy = [qutip.entropy_vn(states[i]) for i in range(len(states))]

	#save the calculated single scalar parameters
	#calculate the average photon number after 3*kappa
	steady_state_photons = photons[t > 3*kappa]
	intracavity_photon_number = np.mean(steady_state_photons)

	#save the calculated parameters
	parameters_to_save['intracavity_photon_number'] = intracavity_photon_number

	#adding the parameters to save
	if SAVE_ARRAYS:
		parameters_to_save['t'] = str(t)
		parameters_to_save['sx'] = str(sx)
		parameters_to_save['sy'] = str(sy)
		parameters_to_save['sz'] = str(sz)
		parameters_to_save['sx^2'] = str(sxx)
		parameters_to_save['sy^2'] = str(syy)
		parameters_to_save['sz^2'] = str(szz)
		parameters_to_save['norm'] = str(norm)
		parameters_to_save['photons'] = str(photons)
		parameters_to_save['sy*sz'] = str(syz)
		parameters_to_save['sz*sx'] = str(szx)
		parameters_to_save['sx*sy'] = str(sxy)
		parameters_to_save['entropy'] = str(entropy)
		parameters_to_save['spin_state'] = str(spin_state)
		parameters_to_save['spin_command'] = str(spin_command)

	return parameters_to_save

