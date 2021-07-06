from wholecell.optimization.base_solver import BaseSolver


class SPSA(BaseSolver):
    def __init__(self, method, sim_dir, **options):
        super().__init__(method, sim_dir, **options)

        self.lr = 0.1
        # TODO: handle solver specific options

    def parameter_updates(self, original_values, objectives, paths):
        # TODO: use actual SPSA update algorithm
        for param, original_value in original_values.items():
            objective_diff = objectives[1] - objectives[0]
            parameter_diff = self.get_param(param, paths[1]) - self.get_param(param, paths[0])
            original_values[param] -= self.lr * objective_diff / parameter_diff * original_value

    def get_parameter_perturbations(self, iteration, index):
        raw_data_perturbations = {}
        sim_data_perturbations = {}

        # TODO: generalize
        if index == 0:
            sim_data_perturbations['constants.a'] = self._method.sim_params['constants.a'] * 1.1
        else:
            sim_data_perturbations['constants.a'] = self._method.sim_params['constants.a'] * 0.9

        return raw_data_perturbations, sim_data_perturbations

    def n_variants_per_iteration(self):
        return 2
