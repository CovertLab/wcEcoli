import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import copy

class TypeIIToxinAntitoxin(object):
	def __init__(self, params):
		pass
		# TODO figure out params!

	def evolve_state(self, time_step, state):
		A = state['A']
		T = state['T']
		kt = state['kt']
		n = state['n']
		kp1 = state['kp1']
		kp2 = state['kp2']
		kh = state['kh']
		p = state['p']
		sigma = state['sigma']
		alpha = state['alpha']
		umax = state['umax']
		lambda_A = state['lambda_A']
		lambda_T = state['lambda_T']
		d = state['d']


		X = 1 + (T**n)/(kt**n)
		Y = 1 + (A**2)/(kp1 ** 2) + ((2 * A**2 * T) / (kp2**2 * kh))**p + (A**2 * T**2)/(kp1**2 * kh**2)

		xi_A = d * np.random.normal(0,1)
		xi_T = d * np.random.normal(0,1)

		dAdt = (sigma * alpha) / (X * Y) - umax/X * A - lambda_A * A + xi_A
		dTdt = alpha / (X * Y) - umax/X * T - lambda_T * T + xi_T

		A = A + dAdt * time_step
		T = T + dTdt * time_step

		state['A'] = A
		state['T'] = T
		state['dTdt'] = dTdt
		state['dAdt'] = dAdt

		return state


class Integrator(object):
	def __init__(self, system, params):
		self.system = system
		self.time_step = params.get('time_step', 1.0)

	def integrate(self, time_step, dTdt_cutoff, state, samples=None):
		time = 0
		t_vec = [0]
		dTdt = 1
		dAdt = 1

		series = {}
		if samples is not None:
			series = {key : [state[key]] for key in samples}


		# while (dTdt_cutoff < dTdt) and (dTdt_cutoff < dAdt):
		while dTdt_cutoff < (dTdt + dAdt):
			state = self.system.evolve_state(time_step, state)
			for key in series:
				series[key].append(state[key])

			dTdt = np.abs(state['dTdt']) 
			dAdt = np.abs(state['dAdt'])

			time += time_step
			t_vec.append(time)
			

		return state, series, t_vec

def hysteresis_graph(param_name, param_range, AT_inital, color='r'):
	system = TypeIIToxinAntitoxin({})
	integrator = Integrator(system, {})


	state = {
			'A' : AT_inital,
			'T' : AT_inital,
			'lambda_A' : (np.log(2)/60),
			'lambda_T' : np.log(2)/(48*60)*0,
			'sigma' : 10.0,
			'alpha' : 1.0,
			'kh' : 100.0,
			'kp1' : 1000.0,
			'kp2' : 10.0,
			'kt' : 10.0,
			'd' : 0, # noise scaling parameter
			'p' : 2.0,
			'n' : 2.5, 
			'umax' : (np.log(2)/30),  # umax/umax0 = 0.75
			'dTdt' : 1,
			'dAdt' : 1} 
	# ---------------------- reproduce figure 3 from 2013 paper
	time_step = 0.01
	dTdt_cutoff = 1e-3
	evolved_state = copy.deepcopy(state)

	param_ratio = []
	T = []
	for m in np.linspace(param_range[0], param_range[1] , 20):
		evolved_state = copy.deepcopy(state)
		evolved_state[param_name] = m * state[param_name]
		evolved_state, series, t_vec = integrator.integrate(time_step, dTdt_cutoff, evolved_state, ['T', 'dTdt'])  
		param_ratio.append(m)
		T.append(series['T'][-1])

	plt.plot(param_ratio, T, color=color)


# ------------------------------------------------------------
if __name__ == '__main__':
	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

	plt.figure()

	hysteresis_graph('umax', [0.2, 1.2], 200)
	hysteresis_graph('umax', [1.2, 0.2], 1, 'b')

	plt.xlim((1.2,0.2))
	plt.ylabel('[T]')
	plt.xlabel('umax/umax0')
	plt.savefig('plots/fig3B_2013.png')
	plt.clf()

# # lambda A
# 	plt.figure()
# 	p_range = [2**-8, 2**8]
# 	hysteresis_graph('lambda_A', [p_range[1], p_range[0]], 200)
# 	hysteresis_graph('lambda_A', [p_range[0], p_range[1]], 1, 'b')

# 	plt.xlim((p_range[0], p_range[1]))
# 	plt.xscale('log')
# 	plt.ylabel('[T]')
# 	plt.xlabel('lambda A/lambda A0')
# 	plt.savefig('plots/fig3A_2013.png')
# 	plt.clf()

# # kh
# 	plt.figure()
# 	p_range = [2**-8, 2**8]
# 	hysteresis_graph('kh', [p_range[1], p_range[0]], 200)
# 	hysteresis_graph('kh', [p_range[0], p_range[1]], 1, 'b')

# 	plt.xlim((p_range[0], p_range[1]))
# 	plt.ylabel('[T]')
# 	plt.xlabel('kh/kh0')
# 	plt.savefig('plots/kh_vs_T.png')
# 	plt.clf()