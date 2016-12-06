DIVISION_MASS_VARIANCE = [0., 0.3, 0.5]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def divisionMassVarianceTotalIndices(sim_data):
	return len(DIVISION_MASS_VARIANCE) + 1


def divisionMassVariance(sim_data, index):

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	sim_data.divisionMassVariance = DIVISION_MASS_VARIANCE[index - 1]

	return dict(
		shortName = "{}".format(DIVISION_MASS_VARIANCE[index - 1]),
		desc = "Division mass variance is {}".format(DIVISION_MASS_VARIANCE[index - 1])
		), sim_data
