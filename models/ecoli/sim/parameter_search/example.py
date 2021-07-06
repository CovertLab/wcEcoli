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

    def get_objective(self, sim_out_dirs):
        return 0
