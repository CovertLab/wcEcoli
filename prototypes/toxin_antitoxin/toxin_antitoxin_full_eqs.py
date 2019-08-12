import numpy as np
class TypeIIToxinAntitoxin(object):
	def __init__(self, params):
		pass
		# TODO figure out params!

	def evolve_state(self, time_step, state):
		A = state['A']
		T = state['T']
		T2 = state['T2']
		A2 = state['A2']
		A2T = state['A2T']
		A2T2 = state['A2T2']
		kt1 = state['kt1']
		kt2 = state['kt2']
		kp1 = state['kp1']
		kp2 = state['kp2'] 
		kh = state['kh']
		tu = state['tu']
		ta = state['ta'] 
		tt = state['tt']
		alpha = state['alpha']
		rho = state['rho']
		n = state['n'] 
		sigma = state['sigma'] 
		p = state['p']

		ka2f = state['ka2f']
		ka2r = state['ka2r']
		ka2tf = state['ka2tf']
		ka2tr = state['ka2tr']
		ka2t2f = state['ka2t2f']
		ka2t2r = state['ka2t2r']

		# convert half lives to rates
		umax = np.log(2)/tu
		lambda_A = np.log(2)/ta
		lambda_T = np.log(2)/tt

		# Eqs from Fasani and Savageau
		X1 = 1 + T**n/kt1**n #Hill equation
		X2 = 1 + T**n/kt2**n #Hill equation
		# Y = 1 + A2/kp1 + ((2 * A2T)/(kp2 * kh))**p + (A2T2)/(kp1 * kh)
		Y = 1 + (A**2)/(kp1**2) + ((2 * A**2 * T)/(kp2**2 * kh))**p + (A**2 * T**2)/(kp1**2 * kh**2)


		dAdt = sigma * alpha /(Y * X2) - umax * A/X1 - lambda_A*A + np.random.normal(0,1) * 10**(-1)
		dTdt = alpha/(Y * X2) - umax * T/X1 - lambda_T*T + np.random.normal(0,1) * 10**(-1)

		# complexes
		# A2 = A2 + ka2f * A - ka2r*A2
		# A2T = A2T + ka2tf * A2 * T - ka2tr*A2T
		# A2T2 = A2T2 + ka2t2f * A2T * T - ka2t2r * A2T2



		state['A'] = state['A'] + dAdt
		state['T'] = state['T'] + dTdt
		state['dAdt'] = dAdt
		state['dTdt'] = dTdt

		state['A2'] = A2
		state['A2T'] = A2T
		state['A2T2'] = A2T2
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
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	system = TypeIIToxinAntitoxin({})
	integrator = Integrator(system, {})
	state = {
		'A' : 100,
		'T' : 100,
		'T2' : 10,
		'A2' : 10,
		'A2T' : 10,
		'A2T2' : 10,
		'kt1' : 10,
		'kt2' : 100,
		'kp1' : 1,
		'kp2' : 10,
		'kh' : 100,
		'tu' : 30,
		'ta' : 60,  # min
		'tt' : 2880,  # 48hrs = 2880 min
		'alpha' : 1,
		'rho' : 2,
		'n' : 2,
		'sigma' : 10,
		'p' : 5,
		'dAdt' : 0,
		'dTdt' : 0,
		# totally made up stuff
		'ka2f' : 100,   # A -> A2
		'ka2r' : 10,    # A <- A2
		'ka2tf' : 100,  # A2 + T -> A2T
		'ka2tr' : 10,   # A2 + T <- A2T
		'ka2t2f' : 100, # A2T + T -> A2T2
		'ka2t2r' : 10}  # A2T + T <- A2T2

	T = []
	A = []
	dAdt = []
	dTdt = []	
	time_vec = range(0,1000)
	for t in time_vec:	
		state = integrator.integrate(1.0, 2, state)
		T.append(state['T'])
		A.append(state['A'])
		dAdt.append(state['dAdt'])
		dTdt.append(state['dTdt'])

	plt.plot(time_vec, T, label='T')
	plt.plot(time_vec, A, label='A')
	plt.legend()
	plt.savefig('AT')
	plt.clf()

	plt.plot(time_vec, dAdt, label='dAdt')
	plt.plot(time_vec, dTdt, label='dTdt')
	plt.legend()
	plt.savefig('AT_derivatives')

	
