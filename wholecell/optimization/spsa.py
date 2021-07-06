from wholecell.optimization.base_solver import BaseSolver


class SPSA(BaseSolver):
    def __init__(self, method, sim_dir, **options):
        super().__init__(method, sim_dir, **options)

        # TODO: handle solver specific options

    def parameter_update(self):
        raise NotImplementedError('Need to implement in a subclass.')

    def get_parameter_perturbations(self):
        return [{}]*2, [{}]*2

    def n_variants_per_iteration(self):
        return 2
