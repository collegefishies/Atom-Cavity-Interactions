import qutip
from qutip import mesolve
from qutip.ui.progressbar import TextProgressBar, BaseProgressBar
import numpy as np

SAVE_ARRAYS = False

GAMMA_SI = 180e3/2/np.pi
KAPPA_SI = 500E3/2/np.pi
ETA = 20
G_SI = 2*np.pi*np.sqrt(4*ETA*GAMMA_SI*KAPPA_SI)
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

def rotated_variance(sx, sy, sz, sxy, syz, szx, sxx, syy, szz, axis, angle):
	"""Calculates the rotated moments over time of a spin-1/2 system from its expectation values of spin operators.
	
	Args:
		sx (ndarray): Array of expectation values of sigma_x.
		sy (ndarray): Array of expectation values of sigma_y.
		sz (ndarray): Array of expectation values of sigma_z.
		sxx (float): Expectation value of sigma_x^2.
		syy (float): Expectation value of sigma_y^2.
		szz (float): Expectation value of sigma_z^2.
		axis (ndarray): String defining the rotation axis 'x', 'y', 'z'.
		angle (float): Angle of rotation in radians.
		
	Returns:
		float: The rotated variance.
	"""
	
	if axis == 'x' or axis == 'X':
		s0  = sy
		s1  = sz
		s01 = sxy
		s00 = syy
		s11 = szz
	elif axis == 'y' or axis == 'Y':
		s0  = sz
		s1  = sx
		s01 = szx
		s00 = szz
		s11 = sxx
	elif axis == 'z' or axis == 'Z':
		s0  = sx
		s1  = sy
		s01 = sxy
		s00 = sxx
		s11 = syy
	else:
		raise ValueError("Invalid rotation axis. Must be 'x', 'y', or 'z'.")
		
	
	#calculate the rotation moments sa, saa, sb, sbb. For 'x', sa = s0, sb = s1 at no rotation
	#sa = s0*cos(angle) + s1*sin(angle)
	#saa = s00*cos(angle)^2 + s11*sin(angle)^2 + (s01+s10)*cos(angle)*sin(angle).
	
	c = np.cos(angle)
	s = np.sin(angle)
	
	sa = s0 * c + s1 * s
	sb = -s0 * s + s1 * c
	
	saa = s00 * c ** 2 + s11 * s ** 2 +  s01 * c * s
	sbb = s00 * s ** 2 + s11 * c ** 2 -  s01 * c * s
	
	#return the variances over time
	return (saa - sa**2, sbb - sb**2)
	 
def max_min_variance(sx, sy, sz, sxy, syz, szx, sxx, syy, szz, axis):
	variances = []
	for angle in np.linspace(0,2*np.pi, 180):
		x,y = rotated_variance(sx, sy, sz, sxy, syz, szx, sxx, syy, szz, axis, angle)
		
		variances.append(x)
		
	#convert this to a numpy array
	variances = np.array(variances)


	#and take max over the new dimension
	max_variance_over_time = np.max(variances, axis=0)
	min_variance_over_time = np.min(variances, axis=0)
		
		
	return max_variance_over_time, min_variance_over_time

def simulation(Natoms=3,detuning=0,
	Gamma_si=GAMMA_SI, 
	kappa_si=KAPPA_SI,
	eta=ETA, 
	g_si=G_SI,
	omega_si=OMEGA_SI,
	Lambda=LAMBDA_SI, #this parameter is the SHO driving strength
	scale=SCALE,
	spin_command='',
	T=60,
	SAVE_ARRAYS=SAVE_ARRAYS):
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

	#calculate the min max var
	axis = 'Y'
	s_maxmax, s_minmin  = max_min_variance(sx, sy, sz, sxy, syz, szx, sxx, syy, szz, axis)
	stdVar = Natoms/4
	Qsq = np.sqrt(s_maxmax/s_minmin)
	contrast = np.sqrt(sx**2+sy**2 + sz**2)/(Natoms/2) #normalized to maximum spin vector length
	area = np.sqrt(s_maxmax*s_minmin)/stdVar #normalized to the CSS std dev.
	wineland_parameter = Qsq*contrast**2/(area)
	max_wineland_time = t[np.argmax(wineland_parameter)]

	#save the calculated single scalar parameters
	#calculate the average photon number after 3*kappa
	steady_state_photons = photons[t > 3*kappa]
	intracavity_photon_number = np.mean(steady_state_photons)

	#save the calculated parameters
	parameters_to_save['intracavity_photon_number'] = intracavity_photon_number
	parameters_to_save['spin_state'] = str(spin_state)
	parameters_to_save['spin_command'] = str(spin_command)
	parameters_to_save['max_wineland_parameter'] = np.max(wineland_parameter)
	parameters_to_save['max_wineland_time'] = max_wineland_time

	print("Max Wineland Parameter: ", parameters_to_save['max_wineland_parameter'])
	print("Max Wineland Time: ", parameters_to_save['max_wineland_time'])
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

	return parameters_to_save

class JuliaSimulator():
	def __init__(self):
		self._import_packages()
		self._build_quantum_system()
		self._define_hamiltonian()
		self._calculate_equations()
		self._get_ode()

	def _import_packages(self):
		#import the necessary packages
		pkgs = {
			"QuantumCumulants": "0.2.13",
			"ModelingToolkit": "8.21.0",
			"OrdinaryDiffEq": "6.11.2",
			"PyPlot": "2.11.0"
		}
		jl.eval("using Pkg"); commands = []
		for key,val in pkgs.items():
			commands.append(f'Pkg.add(PackageSpec(name="{key}", version="{val}"))')
		for cmd in commands:
			print(cmd)
		for cmd in commands:
			jl.eval(cmd)
	def _build_quantum_system(self):
		jl.eval("hc = FockSpace(:cavity)")
		jl.eval("ha = NLevelSpace(:atom, 3)")
		jl.eval("h = hc ⊗ ha")
		jl.eval("Transition(h,:σ,1,1)")

		jl.eval("@qnumbers a::Destroy(h)")
		jl.eval("σ(α,β,i) = IndexedOperator(Transition(h,:σ,α,β),i)")
		jl.eval("@cnumbers Γ κ Ω g λ N Δc Δ2 Δ3")
		# jl.eval("@syms t::Real")
		# jl.eval("@register_symbolic Ωt(t)")
		# jl.eval("@register_symbolic λt(t)")

		i = jl.eval("Index(h,:i,N,ha)")
		j = jl.eval("Index(h,:j,N,ha)")

	def _define_hamiltonian(self):
		println("Adding Hamiltonian")
		jl.eval("H0 = -Δc*a'a - Δ2*∑(σ(2,2,i),i) - Δ3*∑(σ(3,3,i),i)")
		jl.eval("Hda = Ω*∑( σ(2,1,i) + σ(1,2,i), i)")
		jl.eval("Hdc = λ*(a' + a)")
		jl.eval("Hint = g*∑( (a'*σ(2,3,i) + a*(σ(3,2,i))) ,i)")
		jl.eval("H = H0 + Hda + Hdc + Hint")

		jl.eval("J = [a, σ(2,3,i)]")
		jl.eval("rates = [κ, Γ]")

	def _calculate_equations(self):
		print("Calculating equations")
		jl.eval("ops = [a'a, σ(2,2,j)]")
		jl.eval("eqs = meanfield(ops, H, J;rates, order=2)")
		jl.eval("eqs_c = complete(eqs)")
		jl.eval("eqs_sc = scale(eqs_c)")
		jl.eval("println("Length of eqs_sc: ", length(eqs_sc))")

	def _get_ode(self):
		jl.eval("@named sys = ODESystem(eqs_sc)")