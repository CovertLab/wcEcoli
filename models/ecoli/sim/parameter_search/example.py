import pickle

from models.ecoli.sim.parameter_search.base_parameter_search import BaseParameterSearch


class Example(BaseParameterSearch):
    raw_params = ()
    sim_params = ()
    sims_to_run = (
        {
            'jit': False,
            'length_sec': 2,
        },
    )

    def update_raw_data(self, objectives, paths):
        print(self.raw_params)

    def update_sim_data(self, objectives, paths):
        print(self.raw_params)

    def get_objective(self, sim_out_dirs, sim_data_files):
        objectives = []
        for sim_data_file in sim_data_files:
            with open(sim_data_file, 'rb') as f:
                sim_data = pickle.load(f)
            objectives.append(self.get_attrs(sim_data, 'constants.a'))

        return objectives

    def initial_sim_data(self, iteration, index):
        return {'constants.a': 1}
