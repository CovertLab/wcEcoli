class TypeIIToxinAntitoxin(object):
	def __init__(self, params):
		pass
		# TODO figure out params!

	def evolve_state(self, time_step, state):
		T = state['T']
		KT1 = state['KT1']

		X1 = 1 + T/KT1

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
	state = {'T' : 23, 'KT1' : 4, 'X1' : 2}
	post = integrator.integrate(1.0, 20, state)
	print(post)
