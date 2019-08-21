import numpy as np
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
	from matplotlib import cm
	import copy

	system = TypeIIToxinAntitoxin({})
	integrator = Integrator(system, {})

# ---------------------- ORIGINAL PARAMS
	state = {
			'A' : 100.0,
			'T' : 100.0,
			'lambda_A' : (np.log(2)/60),
			'lambda_T' : np.log(2)/(48*60),
			'sigma' : 10.0,
			'alpha' : 1.0,
			'kh' : 100.0,
			'kp1' : 1000.0,
			'kp2' : 10.0,
			'kt' : 10.0,
			'd' : 0, # noise scaling parameter
			'p' : 2.0,
			'n' : 2.5,
			'umax' : 0.75*(np.log(2)/30)} # umax/umax0 = 0.75

	time_step = 0.1
	duration = 2000
	n_sims = 1

	for n in range(0, n_sims):
		evolved_state = copy.deepcopy(state)
		evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
		plt.plot(t_vec, series['T'])
		# print(series['T'][-1])

	plt.ylabel('[T]')
	plt.savefig('plots/T_original.png')
	plt.clf()

# ## -------------------------- alternate params from table 1: S2
# 	state = {
# 			'A' : 55,
# 			'T' : 55,
# 			'lambda_A' : np.log(2)/120,
# 			'lambda_T' : np.log(2)/(96*60),
# 			'sigma' : 10,
# 			'alpha' : 1,
# 			'kh' : 100,
# 			'kp1' : 1000,
# 			'kp2' : 10,
# 			'kt' : 10,
# 			'd' : 0.2,
# 			'p' : 2,
# 			'n' : 2,
# 			'umax' : np.log(2)/30}

# 	time_step = 0.1
# 	duration = 2000
# 	n_sims = 20

# 	for n in range(0, n_sims):
# 		evolved_state = copy.deepcopy(state)
# 		evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 		plt.plot(t_vec, series['T'])

# 	plt.savefig('plots/T_S2.png')
# 	plt.clf()

# ## --------------------------------- alternate params from table 1: S3
# 	state = {
# 			'A' : 40,
# 			'T' : 40,
# 			'lambda_A' : np.log(2)/240,
# 			'lambda_T' : np.log(2)/(192*60),
# 			'sigma' : 2.5,
# 			'alpha' : 1,
# 			'kh' : 100,
# 			'kp1' : 1000,
# 			'kp2' : 2.5,
# 			'kt' : 10,
# 			'd' : 0.2,
# 			'p' : 2,
# 			'n' : 2,
# 			'umax' : np.log(2)/30}

# 	time_step = 0.1
# 	duration = 2000
# 	n_sims = 20

# 	for n in range(0, n_sims):
# 		evolved_state = copy.deepcopy(state)
# 		evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 		plt.plot(t_vec, series['T'])

# 	plt.savefig('plots/T_S3.png')
# 	plt.clf()

# # --------------------------------- vary lambda_A
# N = 15
# cm_vals = [cm.rainbow(x) for x in np.linspace(0.0, 1.0, N)]
# k = 0
# for m in np.linspace(.00001, 60, N):
# 	state = {
# 			'A' : 55,
# 			'T' : 55,
# 			'lambda_A' : np.log(2)/m,
# 			'lambda_T' : np.log(2)/(48*60),
# 			'sigma' : 10,
# 			'alpha' : 1,
# 			'kh' : 100,
# 			'kp1' : 1000,
# 			'kp2' : 10,
# 			'kt' : 10,
# 			'd' : 0,
# 			'p' : 2,
# 			'n' : 2,
# 			'umax' : np.log(2)/30}

# 	A_deg = np.log(2)/state['lambda_A']
# 	line_label = '{0:.4g}'.format(A_deg) + ' min'
			
# 	time_step = 0.1
# 	duration = 5000
# 	n_sims = 1

# 	evolved_state = copy.deepcopy(state)
# 	evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 	plt.plot(t_vec, series['T'], color=cm_vals[k], label=line_label)
# 	k += 1
# plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# plt.ylabel('[T]')
# plt.title('ODE with varied lambda A')
# plt.savefig('plots/vary_lambdaA.png', bbox_inches = "tight")
# plt.clf()

# # --------------------------------- vary lambda_T
# N = 15
# cm_vals = [cm.rainbow(x) for x in np.linspace(0.0, 1.0, N)]
# k = 0
# for m in np.linspace(.5, 480, N):
# 	state = {
# 			'A' : 55.0,
# 			'T' : 55.0,
# 			'lambda_A' : np.log(2)/60,
# 			'lambda_T' : np.log(2)/(m*60),
# 			'sigma' : 10.0,
# 			'alpha' : 1.0,
# 			'kh' : 100.0,
# 			'kp1' : 1000.0,
# 			'kp2' : 10.0,
# 			'kt' : 10.0,
# 			'd' : 0,
# 			'p' : 2,
# 			'n' : 2,
# 			'umax' : np.log(2)/30}


# 	line_label = '{0:.4g}'.format(m) + ' hrs'
			
# 	time_step = 0.1
# 	duration = 5000
# 	n_sims = 1

# 	evolved_state = copy.deepcopy(state)
# 	evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 	plt.plot(t_vec, series['T'], color=cm_vals[k], label=line_label)
# 	k += 1
# plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# plt.ylabel('[T]')
# plt.title('ODE with varied lambda T')
# plt.savefig('plots/vary_lambdaT.png', bbox_inches = "tight")
# plt.clf()

# # --------------------------------- vary sigma
# N = 15
# cm_vals = [cm.rainbow(x) for x in np.linspace(0.0, 1.0, N)]
# k = 0
# for m in np.linspace(.000001, 10, N):
# 	state = {
# 			'A' : 55,
# 			'T' : 55,
# 			'lambda_A' : np.log(2)/60,
# 			'lambda_T' : np.log(2)/(48*60),
# 			'sigma' : m,
# 			'alpha' : 1,
# 			'kh' : 100,
# 			'kp1' : 1000,
# 			'kp2' : 10,
# 			'kt' : 10,
# 			'd' : 0,
# 			'p' : 2,
# 			'n' : 2,
# 			'umax' : np.log(2)/30}

# 	line_label = '{0:.4g}'.format(m) 
			
# 	time_step = 0.1
# 	duration = 5000
# 	n_sims = 1

# 	evolved_state = copy.deepcopy(state)
# 	evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 	plt.plot(t_vec, series['T'], color=cm_vals[k], label=line_label)
# 	k += 1
# plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# plt.ylabel('[T]')
# plt.title('ODE with varied sigma')
# plt.savefig('plots/vary_sigma.png', bbox_inches = "tight")
# plt.clf()

# # --------------------------------- k params

# params_to_vary = {'kh' : 100, 'kp1' : 1000,'kp2' : 10, 'kt' : 10}

# N = 15
# cm_vals = [cm.rainbow(x) for x in np.linspace(0.0, 1.0, N)]
# for param in params_to_vary.keys():
# 	state = {
# 		'A' : 55,
# 		'T' : 55,
# 		'lambda_A' : np.log(2)/60,
# 		'lambda_T' : np.log(2)/(48*60),
# 		'sigma' : 10,
# 		'alpha' : 1,
# 		'kh' : 100,
# 		'kp1' : 1000,
# 		'kp2' : 10,
# 		'kt' : 10,
# 		'd' : 0,
# 		'p' : 2,
# 		'n' : 2,
# 		'umax' : np.log(2)/30}

# 	time_step = 0.1
# 	duration = 5000
# 	evolved_state = copy.deepcopy(state)
# 	evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 	plt.plot(t_vec, series['T'], color='k', label='default', linewidth=2)
	
# 	k = 0
# 	for m in np.logspace( 10**-5, 5, N):

# 		state[param] = m

# 		line_label = '{0:.4g}'.format(m) 

# 		evolved_state = copy.deepcopy(state)
# 		evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 		plt.plot(t_vec, series['T'], color=cm_vals[k], label=line_label)
# 		k += 1
# 	plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# 	plt.ylabel('[T]')
# 	plt.title('ODE with varied ' + param + ', default = ' + str(params_to_vary[param]))
# 	plt.savefig('plots/vary_' + param + '.png', bbox_inches = "tight")
# 	plt.clf()

# # --------------------------------- vary A and T initial
# params_to_vary = {'A' : 55, 'T' : 55}

# N = 10
# cm_vals = [cm.rainbow(x) for x in np.linspace(0.0, 1.0, N)]
# for param in params_to_vary.keys():
# 	state = {
# 		'A' : 55,
# 		'T' : 55,
# 		'lambda_A' : np.log(2)/60,
# 		'lambda_T' : np.log(2)/(48*60),
# 		'sigma' : 10,
# 		'alpha' : 1,
# 		'kh' : 100,
# 		'kp1' : 1000,
# 		'kp2' : 10,
# 		'kt' : 10,
# 		'd' : 0,
# 		'p' : 2,
# 		'n' : 2,
# 		'umax' : np.log(2)/30}

# 	time_step = 0.01
# 	duration = 2000
# 	evolved_state = copy.deepcopy(state)
# 	evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 	plt.plot(t_vec, series['T'], color='k', label='default', linewidth=2)
	
# 	k = 0
# 	for m in np.linspace( 1,200, N):

# 		state[param] = m

# 		line_label = '{0:.4g}'.format(m) 

# 		evolved_state = copy.deepcopy(state)
# 		evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 		plt.plot(t_vec, series['T'], color=cm_vals[k], label=line_label)
# 		k += 1
# 	plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
# 	plt.ylabel('[T]')
# 	plt.title('ODE with varied ' + param + ', default = ' + str(params_to_vary[param]))
# 	plt.savefig('plots/vary_' + param + '.png', bbox_inches = "tight")
# 	plt.clf()

# # --------------------------------- T vs lambda_A 
# N = 20
# cm_vals = [cm.rainbow(x) for x in np.linspace(0.0, 1.0, N)]
# k = 0
# A_ratio = []
# T_ss = []
# for m in np.logspace(np.log10(1/256.0), np.log10(256), N):
# 	state = {
# 			'A' : 55,
# 			'T' : 55,
# 			'lambda_A' : np.log(2)/(60*m),
# 			'lambda_T' : np.log(2)/(48*60),
# 			'sigma' : 10.0,
# 			'alpha' : 1.0,
# 			'kh' : 100.0,
# 			'kp1' : 1000.0,
# 			'kp2' : 10.0,
# 			'kt' : 10.0,
# 			'd' : 0,
# 			'p' : 2.0,
# 			'n' : 2.0,
# 			'umax' : np.log(2)/30}

# 	A_ratio.append(np.log2(state['lambda_A']/(np.log(2)/60)))
			
# 	time_step = 1
# 	duration = 2000
# 	n_sims = 1

# 	evolved_state = copy.deepcopy(state)
# 	evolved_state, series, t_vec = integrator.integrate(time_step, duration, evolved_state, ['T'])  
# 	T_ss.append(np.log10(series['T'][-1]))
# 	k += 1

# # print(A_ratio)
# # print(T_ss)
# plt.plot(A_ratio, T_ss, marker='.')
# plt.xlabel('log2 lambdaA / lambdaA0')
# plt.ylabel('[T] steady state')
# plt.ylabel('log10 [T]')
# plt.xlim((-8,8))
# plt.title('T vs lambda A')
# plt.savefig('plots/lambdaA_vs_Tss.png', bbox_inches = "tight")
# plt.clf()
