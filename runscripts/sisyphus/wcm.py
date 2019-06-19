import os
import pprint

class WcmWorkflow(object):
	def __init__(self, root, options={}):
		self.root = root
		self.options = options
		self.paths = self.rootify({
			'raw-data': 'raw-data',
			'raw-validation-data': 'raw-validation-data',
			'validation-data': 'validation-data',
			'sim-data': 'sim-data'})

	def add_root(self, key):
		return os.path.join(self.root, key)

	def rootify(self, paths):
		return {
			key: self.add_root(path)
			for key, path in paths.iteritems()}

	def paths_for(self, keys):
		return {
			key: self.paths[key]
			for key in keys}

	def add_path(self, key):
		self.paths[key] = self.add_root(key)

	def build_workflow(self, variants, simulations, generations, options={}):
		'''
		Build the WCM workflow with the given options.

		Args:
		    variants (dict(str, list)): A dict where keys are variant functions and
		        values are lists of variant indexes.
		    simulations (int): Number of simulations to run.
		    generations (int): Number of generations to run.
		    options (dict(str, *)): Options, including `branch` which specifies whether
		        or not to make a new simulation for each daughter or just follow a single
		        lineage (which is the default).
		'''

		init_raw_data = {
			'key': self.add_root('init-raw-data'),
			'command': 'init-raw-data',
			'outputs': self.paths_for(['raw-data'])}

		init_raw_validation_data = {
			'key': self.add_root('init-raw-validation-data'),
			'command': 'init-raw-validation-data',
			'outputs': self.paths_for(['raw-validation-data'])}

		init_validation_data = {
			'key': self.add_root('init-validation-data'),
			'command': 'init-validation-data',
			'inputs': self.paths_for(['raw-data', 'raw-validation-data']),
			'outputs': self.paths_for(['validation-data'])}

		fit_sim_data = {
			'key': self.add_root('fit-sim-data'),
			'command': 'fit-sim-data',
			'inputs': self.paths_for(['raw-data']),
			'outputs': self.paths_for(['sim-data'])}

		variants_tasks = []

		for variant, indexes in variants.iteritems():
			for index in indexes:
				output_key = 'sim-data-{}-{}'.format(variant, index)
				self.add_path(output_key)
				key = 'variant-{}'.format(output_key)

				variants_tasks.append({
					'key': self.add_root(key),
					'command': 'variant-sim-data',
					'inputs': {'input-sim-data': self.paths['sim-data']},
					'outputs': {'output-sim-data': self.paths[output_key]},
					'vars': {
						'variant-function': variant,
						'variant-index': index}})

		simulations_tasks = []

		for simulation in range(simulations):
			for variant, indexes in variants.iteritems():
				for variant_index in indexes:
					for generation in range(generations):
						branches = [0]
						if options.get('branch') and generation > 0:
							branches = range(2 ** (generations - 1))

						for branch in branches:
							sim_out_key = 'simulation-{}-variant-{}-{}-generation-{}-daughter-{}'.format(simulation, variant, variant_index, generation, branch)
							endow_a = '{}-endow-{}'.format(sim_out_key, 0)
							endow_b = '{}-endow-{}'.format(sim_out_key, 1)
							daughter = branch // 2
							inherit_key = 'simulation-{}-variant-{}-{}-generation-{}-daughter-{}-endow-{}'.format(
								simulation,
								variant,
								variant_index,
								generation-1,
								daughter,
								daughter % 2)
						
							self.add_path(sim_out_key)
							self.add_path(endow_a)
							self.add_path(endow_b)
							self.add_path(inherit_key)

							key = 'run-{}'.format(sim_out_key)
							variant_key = self.paths['sim-data-{}-{}'.format(variant, variant_index)]
							inputs = {'sim-data': variant_key}
							base_var = {
								'sim-index': "{:06d}".format(simulation),
								'variant-function': variant,
								'variant-index': "{:06d}".format(variant_index),
								'generation': "{:06d}".format(generation)}

							if generation > 0:
								inputs['inherited-state'] = self.paths[inherit_key]

							var = base_var
							if generation > 0:
								var = dict(var, branch="{:06d}".format(branch))

							simulations_tasks.append({
								'key': self.add_root(key),
								'command': 'simulation' if generation == 0 else 'simulation-daughter',
								'inputs': inputs,
								'outputs': {
									'sim-out': self.paths[sim_out_key],
									'daughter-a': self.paths[endow_a],
									'daughter-b': self.paths[endow_b]},
								'vars': var})
		return [
			init_raw_data,
			init_raw_validation_data,
			init_validation_data,
			fit_sim_data] + variants_tasks + simulations_tasks

if __name__ == '__main__':
	workflow = WcmWorkflow('sisyphus:base')
	processes = workflow.build_workflow({'wildtype': [0], 'other': [11, 12]}, 2, 3)
	pp = pprint.PrettyPrinter()
	pp.pprint(processes)
