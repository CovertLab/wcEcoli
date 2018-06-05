# lambda from 0 to 1
LAMBDA = [0] + [10**x for x in range(-8, 1)]

def lambdaWeightIndices(sim_data):
	return len(LAMBDA)

def lambdaWeight(sim_data, index):
	weight = LAMBDA[index]
	sim_data.constants.metabolismKineticObjectiveWeight = weight

	return dict(
		shortName = "lambda={:.0E}".format(weight),
		desc = "Simulation with lambda = {}.".format(weight)
		), sim_data
