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

		umax = np.log(2)/tu
		lambda_A = np.log(2)/ta
		lambda_T = np.log(2)/tt

		X1 = 1 + T**2/kt1**2 #Hill equation, n = 2
		X2 = 1 + T**2/kt2**2 #Hill equation, n = 2
		Y = 1 + A2/kp1 + ((2 * A2 * T)/(kp2 * kh))**p + (A2 * T2)/(kp1 * kh)

		state['dAdt'] = sigma * alpha /(Y * X2) - umax * A/X1 - lambda_A*A 
		state['dTdt'] = alpha/(Y * X2) - umax * T/X1 - lambda_T*T

		state['A'] = state['A'] + state['dAdt']
		state['T'] = state['T'] + state['dTdt']
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
		'ta' : 60,
		'tt' : 48,
		'alpha' : 1,
		'rho' : 2,
		'n' : 2,
		'sigma' : 10,
		'p' : 2,
		'dAdt' : 0,
		'dTdt' : 0} 
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

	
