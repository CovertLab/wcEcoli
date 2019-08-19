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
		# print('X1: ' + str(X1))
		# print('X2: ' + str(X2))
		# print('Y: ' + str(Y))

		xi_A = d * np.random.normal(0,1)
		xi_T = d * np.random.normal(0,1)

		# print('xi_A: ' + str(xi_A))
		
		dAdt = sigma * alpha /(Y * X2) - umax/X1 * A - lambda_A*A + xi_A
		dTdt = alpha/(Y * X2) - umax/X1 * T - lambda_T*T + xi_T


		state['A'] = A + dAdt * time_step 
		state['T'] = T + dTdt * time_step
		return state


class Integrator(object):
	def __init__(self, system, params):
		self.system = system
		self.time_step = params.get('time_step', 1.0)

	def integrate(self, time_step, duration, state, samples=None):
		time = 0
		t_vec = [0]

		series = {}
		if samples is not None:
			series = {key : [state[key]] for key in samples}

		while time < duration:
			state = self.system.evolve_state(time_step, state)
			for key in series:
				series[key].append(state[key])

			time += time_step
			t_vec.append(time)
			

		return state, series, t_vec

if __name__ == '__main__':
	import numpy as np
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import copy

	system = TypeIIToxinAntitoxin({})
	integrator = Integrator(system, {})
	state = {
		'A' : 100.0,
		'T' : 155.0,
		'kt1' : 10.0,
		'kt2' : 10.0,
		'kp1' : 1000.0,  # 1uM = 1000 nM
		'kp2' : 10.0,
		'kh' : 100.0,
		'umax' : np.log(2)/30,  # tu = 30 min
		'lambda_A' : np.log(2)/60,  # ta = 60 min
		'lambda_T' : np.log(2)/2880, # tt = 48hrs = 2880 min
		'alpha' : 1.0, # 1 nM/min
		'n' : 2,  # hill coeff 
		'sigma' : 10.0,
		'p' : 2,
		'd' : 0.5}  # noise magnitude

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



	max_time = 2000
	t_step = 0.1
	num_sims = range(0,10)
	for n in num_sims:
		evolved_state = copy.deepcopy(state)

		evolved_state, series, t_vec = integrator.integrate(t_step, max_time, evolved_state,['T'])
		plt.plot(t_vec,series['T'], alpha=0.2, color=[0.5,0.5,0.5])
	
	plt.savefig('stochastic_sims')


# #--------------  PARAMETER SHIFTING  ------------------------
# # * UMAX
# state = {
# 	'A' : 10,
# 	'T' : 10,
# 	'kt1' : 10,
# 	'kt2' : 100,
# 	'kp1' : 1000,  # 1uM = 1000 nM
# 	'kp2' : 10,
# 	'kh' : 100,
# 	'umax' : np.log(2)/30,  # tu = 30 min
# 	'lambda_A' : np.log(2)/60,  # ta = 60 min
# 	'lambda_T' : np.log(2)/2880, # tt = 48hrs = 2880 min
# 	'alpha' : 1, # 1 nM/min
# 	'n' : 2.5,  # hill coeff 
# 	'sigma' : 10,
# 	'p' : 2,
# 	'd' : 0}  # noise magnitude

# umax = state['umax']
# umax_multipliers = np.linspace(0.2,3,400)
# time_step = 0.1
# evolved_state = copy.deepcopy(state)
# T = []
# umax_ratio = []
# for s in umax_multipliers:
# 	evolved_state['umax'] = state['umax'] * s
# 	# print(evolved_state['umax']/state['umax'])
# 	evolved_state, series = integrator.integrate(time_step, time_step, evolved_state)
# 	T.append(evolved_state['T'])
# 	# print(evolved_state['T'])
# 	umax_ratio.append(s)

# plt.plot(umax_ratio, T)
# plt.savefig('umax_vs_T')
# plt.close()

# # * lambda_A

# state = {
# 	'A' : 10,
# 	'T' : 10,
# 	'kt1' : 10,
# 	'kt2' : 100,
# 	'kp1' : 1000,  # 1uM = 1000 nM
# 	'kp2' : 10,
# 	'kh' : 100,
# 	'umax' : np.log(2)/30,  # tu = 30 min
# 	'lambda_A' : np.log(2)/60,  # ta = 60 min
# 	'lambda_T' : np.log(2)/2880, # tt = 48hrs = 2880 min
# 	'alpha' : 1, # 1 nM/min
# 	'n' : 2.5,  # hill coeff 
# 	'sigma' : 10,
# 	'p' : 2,
# 	'd' : 0}  # noise magnitude

# umax = state['umax']
# # slowly increase
# lambda_A_multipliers = np.linspace(1,2.5,400)
# time_step = 0.1
# evolved_state = copy.deepcopy(state)
# T = []
# lambdaA_ratio = []
# for s in lambda_A_multipliers:
# 	evolved_state['umax'] = state['umax'] * s
# 	# print(evolved_state['umax']/state['umax'])
# 	evolved_state, series = integrator.integrate(time_step, time_step, evolved_state)
# 	T.append(evolved_state['T'])
# 	# print(evolved_state['T'])
# 	lambdaA_ratio.append(s)

# plt.plot(lambdaA_ratio, T, color='b')

# # slowly decrease
# lambda_A_multipliers = np.linspace(2.5,1,400)
# evolved_state = copy.deepcopy(state)
# T = []
# lambdaA_ratio = []
# for s in lambda_A_multipliers:
# 	evolved_state['umax'] = state['umax'] * s
# 	# print(evolved_state['umax']/state['umax'])
# 	evolved_state, series = integrator.integrate(time_step, time_step, evolved_state)
# 	T.append(evolved_state['T'])
# 	# print(evolved_state['T'])
# 	lambdaA_ratio.append(s)

# plt.plot(lambdaA_ratio, T, color='r')
# plt.savefig('lambdaA_ratio_vs_T')
