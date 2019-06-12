import os

class WcmWorkflow(object):
	def __init__(self, root, options):
		self.root = root
		self.options = options
		self.paths = self.rootify({
			'raw-data': 'raw-data',
			'raw-validation-data': 'raw-validation-data'
			'validation-data': 'validation-data'
			'sim-data': 'sim-data'})

	def rootify(self, paths):
		return {
			key: os.path.join(self.root, path)
			for key, path in paths}

	def process(self, key, command, inputs, outputs, var={}):
		outputs = self.rootify(outputs)
		spec = {
			key: os.path.join(self.root, key),
			command: command,
			outputs: outputs}

		if len(inputs) > 0:
			inputs = self.rootify(inputs)
			spec['inputs'] = inputs
		if len(var) > 0:
			spec['vars'] = var

		return spec

	def self.paths_for(self, keys):
		return {
			key: self.paths[key]
			for key in keys}

	def build_workflow(variants, simulations, generations, options):
		init_raw_data = self.process(
			'init-raw-data',
			'init-raw-data',
			{},
			self.paths_for(['raw-data']))

		init_raw_validation_data = self.process(
			os.path.join(self.root, 'init-raw-validation-data'),
			'init-raw-validation-data',
			{},
			self.paths_for(['raw-validation-data']))

		init_validation_data = self.process(
			os.path.join(self.root, 'init-validation-data'),
			'init-validation-data',
			self.paths_for(['raw-data', 'raw-validation-data']),
			self.paths_for(['validation-data']))

		fit_sim_data = self.process(
			'fit-sim-data',
			'fit-sim-data',
			self.paths_for(['raw-data']),
			self.paths_for(['sim-data']))

		sim_datas = [fit_sim_data]

		for index, variant in enumerate(sim_datas):
			key = 'sim-data-{}'.format(index)
			self.paths[key] = os.path.join(self.root, key)
			variant_sim_data.append(self.process(
				key,
				'variant-sim-data',
				{'input-sim-data': self.paths['sim-data']},
				{'output-sim-data': self.paths[key]}
				{'variant-function': variant,
				 'variant-index': index}))

		simulations = []

		for simulation in simulations:
			for variant, sim_data in enumerate(sim_datas):
				for generation in generations:
					branch = 0
					key = 'simulation-{}-variant-{}-generation-{}-daughter-{}'.format(simulation, variant, generation, branch)
					endow_a = 'simulation-{}-variant-{}-generation-{}-daughter-{}-endow-{}'.format(simulation, variant, generation, branch, 0)
					endow_b = 'simulation-{}-variant-{}-generation-{}-daughter-{}-endow-{}'.format(simulation, variant, generation, branch, 1)

					self.paths[key] = os.path.join(self.root, key)
					self.paths[endow_a] = os.path.join(self.root, endow_a)
					self.paths[endow_b] = os.path.join(self.root, endow_b)

					inputs = [{'sim-data': self.paths['sim-data-{}'.format(variant)]}]
					base_var = {'sim-index': "{:06d}".format(simulation)
						   'variant-index': "{:06d}".format(variant),
						   'generation': "{:06d}".format(generation)}

					if generation > 0:
						inputs[0]['inherited-state'] = self.paths[endow_a]
						if options['branch']:
							inputs.append(
								{'sim-data': self.paths['sim-data-{}'.format(variant)],
								 'inherited-state': self.paths[endow_b]})

					for branch, input in enumerate(inputs):
						var = base_var
						if generation > 0:
							var = dict(var, 'branch' = "{:06d}".format(branch))

						simulation_process = self.process(
							key,
							'simulation' if generation == 0 else 'simulation-daughter'
							inputs,
							{'sim-out': self.paths[key],
							 'daughter-a': self.paths[endow_a]
							 'daughter-b': self.paths[endow_b]},
							var)
