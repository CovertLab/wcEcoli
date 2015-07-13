PPGPP_POWER = [2., 2.5, 3., 3.5, 4.]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def ppGppPowerTotalIndices(kb):
	return len(PPGPP_POWER) + 1


def ppGppPower(kb, index):
	# Vary time step used for simulation

	nConditions = ppGppPowerTotalIndices(kb)

	if index == 0:
		return CONTROL_OUTPUT

	kb.ppGpp_power = PPGPP_POWER[index - 1]

	return dict(
		shortName = "{} power".format(PPGPP_POWER[index - 1]),
		desc = "Simulation uses feedback power of {} for ppGpp in transcription initiation.".format(PPGPP_POWER[index - 1])
		)
