import collections

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def tcsActivityTotalIndices(sim_data):
	nNutrientTimeSeries = len(sim_data.tcsToActiveInactiveConds)
	return nNutrientTimeSeries


def tcsActivity(sim_data, index):

	nTcsActivityTimeSeries = tcsActivityTotalIndices(sim_data)

	if index % nTcsActivityTimeSeries == 0:
		return CONTROL_OUTPUT, sim_data

	tcsList = ["basal (no TCS)"] + sorted(sim_data.tcsToActiveInactiveConds)
	tcs = tcsList[(index + 1) // 2]
	tcsStatus = None
	if index % 2 == 1:
		tcsStatus = "active"
	else:
		tcsStatus = "inactive"

	sim_data.condition = tcs + "__" + tcsStatus

	sim_data.nutrientsTimeSeriesLabel = tcs + "__" + tcsStatus
	sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel] = collections.deque()
	sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel].append((
		0.0,
		sim_data.tcsToActiveInactiveConds[tcs][tcsStatus + " nutrients"]
		))

	return dict(
		shortName = "{}_phenotype".format(tcs + "__" + tcsStatus),
		desc = "Simulation of phenotype {}.".format(tcs + "__" + tcsStatus)
		), sim_data

