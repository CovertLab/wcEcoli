from wholecell.utils import units
from reconstruction.ecoli.fitkb1 import fitKb_1
#DOUBLING_TIME = [25. * units.min, 27. * units.min, 30. * units.min, 33. * units.min, 35. * units.min]
DOUBLING_TIME = [80. * units.min, 90. * units.min]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def growthRateTotalIndices(kb):
	nTimeSteps = len(DOUBLING_TIME)
	nConditions = nTimeSteps + 1
	return nConditions


def growthRate(kb, index):
	# Vary time step used for simulation

	nConditions = growthRateTotalIndices(kb)

	if index == 0:
		return CONTROL_OUTPUT

	doubling_time = DOUBLING_TIME[index - 1]

	fitKb_1(kb, doubling_time)

	return dict(
		shortName = "{} sec".format(DOUBLING_TIME[index - 1]),
		desc = "Simulation uses growth rate of {} minutes.".format(DOUBLING_TIME[index - 1])
		)
