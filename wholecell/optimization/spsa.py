from wholecell.optimization.base_solver import BaseSolver


class SPSA(BaseSolver):
    def __init__(self, method, sim_dir, **options):
        super().__init__(method, sim_dir, **options)

        # TODO: handle solver specific options

    def parameter_update(self, param, objectives, paths, reference_path):
        objective_diff = objectives[1] - objectives[0]
        parameter_diff = self.get_param(param, paths[1]) - self.get_param(param, paths[0])
        original_val = self.get_param(param, reference_path)

        return objective_diff / parameter_diff * original_val

    def get_parameter_perturbations(self):
        return [{}]*2, [{}]*2

    def n_variants_per_iteration(self):
        return 2
