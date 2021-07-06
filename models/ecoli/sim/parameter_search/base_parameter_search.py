"""

"""

from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS


DEFAULT_CLI_KWARGS = {
    'init_sims': 1,
    'generations': 1,
    }

class BaseParameterSearch():
    parca_args = {'cpus': 8}
    raw_params = ()
    sim_params = ()
    sims_to_run = ()

    def get_sim_params(self, sim_dir, variants):
        all_params = []
        for variant in variants:
            for index, sim_params in enumerate(self.sims_to_run):
                params = DEFAULT_SIMULATION_KWARGS.copy()
                params.update(DEFAULT_CLI_KWARGS)
                params.update(sim_params)

                params['variant type'] = self.__class__.__name__
                params['variant'] = variant
                params['sim dir'] = sim_dir
                params['index'] = index

                all_params.append(params)

        return all_params

    def get_objective(self, sim_data_file, sim_out_dirs):
        raise NotImplementedError('Need to implement in a subclass.')
