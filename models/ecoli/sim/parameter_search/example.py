import pickle

from models.ecoli.sim.parameter_search.base_parameter_search import BaseParameterSearch


class Example(BaseParameterSearch):
	_sim_params = ('constants.a',)
	sims_to_run = (
		{
			'jit': False,
			'length_sec': 2,
			},
		)
	_init_sim_params = {'constants.a': 1}

	def update_raw_data(self, objectives, paths):
		print(self.raw_params)

	def update_sim_data(self, objectives, paths):
		print(self.sim_params)

	def get_objective(self, sim_out_dirs, sim_data_files):
		objectives = []
		for sim_data_file in sim_data_files:
			with open(sim_data_file, 'rb') as f:
				sim_data = pickle.load(f)
			objectives.append(self.get_attr(sim_data, 'constants.a')**2)

		return objectives
