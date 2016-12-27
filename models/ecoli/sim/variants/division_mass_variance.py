DIVISION_MASS_VARIANCE = [1., 0.5, 1.5]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def divisionMassVarianceTotalIndices(sim_data):
	return len(DIVISION_MASS_VARIANCE) + 1


def divisionMassVariance(sim_data, index):

	conditionLabels = sorted(sim_data.conditionActiveTfs)
	conditionLabel = conditionLabels[2]
	sim_data.condition = conditionLabel
	# TODO: add new column to condition defs to replace this?
	if sim_data.conditions[conditionLabel]["nutrients"] == "minimal_plus_amino_acids":
		sim_data.nutrientsTimeSeriesLabel = "000003_aa"

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	sim_data.divisionMassVariance = DIVISION_MASS_VARIANCE[index - 1]

	return dict(
		shortName = "{}".format(DIVISION_MASS_VARIANCE[index - 1]),
		desc = "Division mass variance is {}".format(DIVISION_MASS_VARIANCE[index - 1])
		), sim_data
