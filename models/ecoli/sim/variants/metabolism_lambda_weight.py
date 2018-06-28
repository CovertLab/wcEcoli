# lambda from 0 to 1
LAMBDA = [0] + [10**x for x in range(-8, 1)]

def metabolism_lambda_weight_indices(sim_data):
	return len(LAMBDA)

def metabolism_lambda_weight(sim_data, index):
	weight = LAMBDA[index]
	sim_data.constants.metabolismKineticObjectiveWeight = weight

	return dict(
		shortName = "lambda={:.0E}".format(weight),
		desc = "Simulation with lambda = {}.".format(weight)
		), sim_data
