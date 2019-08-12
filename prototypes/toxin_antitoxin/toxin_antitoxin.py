import numpy as np
class TypeIIToxinAntitoxin(object):
	def __init__(self, params):
		pass
		# TODO figure out params!

	def evolve_state(self, time_step, state):
		A = state['A']
		T = state['T']
		kt1 = state['kt1']
		kt2 = state['kt2']
		kp1 = state['kp1'] 
		kp2 = state['kp2'] 
		kh = state['kh']
		umax = state['umax']
		lambda_A = state['lambda_A'] 
		lambda_T = state['lambda_T']
		alpha = state['alpha']
		n = state['n'] 
		sigma = state['sigma'] 
		p = state['p']
		d = state['d']

		X1 = 1 + T**n/kt1**n #Hill equation
		X2 = 1 + T**n/kt2**n #Hill equation
		Y = 1 + A**2/kp1**2 + ((2 * A**2 * T)/(kp2**2 * kh))**p + (A**2 * T**2)/(kp1**2 * kh**2)

		dAdt = sigma * alpha /(Y * X2) - umax * A/X1 - lambda_A*A + d * np.random.normal(0,np.sqrt(time_step)) 
		dTdt = alpha/(Y * X2) - umax * T/X1 - lambda_T*T + d * np.random.normal(0,np.sqrt(time_step)) 

		state['A'] = A + dAdt * time_step 
		state['T'] = T + dTdt * time_step
		return state


class Integrator(object):
	def __init__(self, system, params):
		self.system = system
		self.time_step = params.get('time_step', 1.0)

	def integrate(self, time_step, duration, state):
		time = 0

		while time < duration:
			state = self.system.evolve_state(time_step, state)
			time += time_step

		return state

if __name__ == '__main__':
	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import copy

	system = TypeIIToxinAntitoxin({})
	integrator = Integrator(system, {})
	state = {
		'A' : 7,
		'T' : 7,
		'kt1' : 10,
		'kt2' : 100,
		'kp1' : 1000,  # 1uM = 1000 nM
		'kp2' : 10,
		'kh' : 100,
		'umax' : np.log(2)/30,  # tu = 30 min
		'lambda_A' : np.log(2)/60,  # ta = 60 min
		'lambda_T' : np.log(2)/2880, # tt = 48hrs = 2880 min
		'alpha' : 1, # 1 nM/min
		'n' : 2.5,  # hill coeff 
		'sigma' : 10,
		'p' : 2,
		'd' : 0.2}  # noise magnitude

	# T = []
	# A = []
	# time_vec = range(0,100)
	# evolved_state = state
	# for t in time_vec:	
	# 	evolved_state = integrator.integrate(1.0, 60, evolved_state)
	# 	T.append(evolved_state['T'])
	# 	A.append(evolved_state['A'])

	# plt.plot(time_vec, T, label='T', color='r')
	# plt.plot(time_vec, A, label='A', color='b')
	# plt.xlabel('Time (hrs)')
	# plt.legend()
	# plt.savefig('AT')
	# plt.clf()



	time_vec = range(0,200)
	num_sims = range(0,1000)
	for n in num_sims:
		evolved_state = copy.deepcopy(state)
		T = [evolved_state['T']]
		for t in time_vec[1:]:	
			evolved_state = integrator.integrate(1, 60, evolved_state)
			T.append(evolved_state['T'])

		plt.plot(time_vec,T, alpha=0.2, color=[0.5,0.5,0.5])
	plt.savefig('stochastic_sims')

