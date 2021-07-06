"""

"""

import pickle

from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS


DEFAULT_CLI_KWARGS = {
    'init_sims': 1,
    'generations': 1,
    }

class BaseParameterSearch():
    parca_args = {'cpus': 8}
    _raw_params = ()
    _sim_params = ()
    sims_to_run = ()

    def __init__(self):
        self.variant_name = self.__class__.__name__
        self.raw_params = {p: None for p in self._raw_params}
        self.sim_params = {p: None for p in self._sim_params}
        self.initialized = False

    def update_raw_data(self, objectives, paths):
        raise NotImplementedError('Need to implement in a subclass.')

    def update_sim_data(self, objectives, paths):
        raise NotImplementedError('Need to implement in a subclass.')

    def get_objective(self, sim_out_dirs, sim_data_files):
        raise NotImplementedError('Need to implement in a subclass.')

    def initialize(self, raw_data_file, sim_data_file):
        with open(raw_data_file, 'rb') as f:
            raw_data = pickle.load(f)
        for param in self.raw_params:
            self.raw_params = self.get_attrs(raw_data, param)

        with open(sim_data_file, 'rb') as f:
            sim_data = pickle.load(f)
        for param in self.sim_params:
            self.sim_params = self.get_attrs(sim_data, param)

        self.initialized = True

    def initial_raw_data(self, iteration, index):
        return {}

    def initial_sim_data(self, iteration, index):
        return {}

    def get_sim_params(self, sim_dir, variants):
        all_params = []
        for variant in variants:
            for index, sim_params in enumerate(self.sims_to_run):
                params = DEFAULT_SIMULATION_KWARGS.copy()
                params.update(DEFAULT_CLI_KWARGS)
                params.update(sim_params)

                params['variant type'] = self.variant_name
                params['variant'] = variant
                params['sim dir'] = sim_dir
                params['index'] = index

                all_params.append(params)

        return all_params

    def get_attrs(self, obj, attr):
        attrs = attr.split('.')
        for a in attrs:
            obj = getattr(obj, a)
        return obj

    def print_update(self):
        # TODO: pretty print
        print(self.raw_params)
        print(self.sim_params)
