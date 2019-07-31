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


		X1 = 1 + T2/kt1

		state['X1'] = X1
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
		''} 
	post = integrator.integrate(1.0, 20, state)
	print(post)
