import pickle

from models.ecoli.sim.parameter_search.base_parameter_search import BaseParameterSearch


class Example(BaseParameterSearch):
	_sim_params = ('constants.a',)
	_init_sim_params = {'constants.a': 1}

	def get_objective(self, sim_out_dirs, sim_data_files):
		objectives = []
		for sim_data_file in sim_data_files:
			with open(sim_data_file, 'rb') as f:
				sim_data = pickle.load(f)
			objectives.append(self.get_attr(sim_data, 'constants.a')**2)

		return objectives
