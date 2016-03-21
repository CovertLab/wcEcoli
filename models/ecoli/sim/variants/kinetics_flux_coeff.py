N_COEFFICIENTS = 9
HIGHEST = 100000
LOWEST = 100

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def kineticsFluxCoeffTotalIndices(sim_data):
	nTimeSteps = N_COEFFICIENTS
	nConditions = nTimeSteps + 1
	return nConditions


def kineticsFluxCoeff(sim_data, index):
	# Sets the coefficient of the upper limit of flux used for enzymeKinetics estimates

	nConditions = kineticsFluxCoeffTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	fluxLimit = float(LOWEST) + ((index - 1) * float(HIGHEST - LOWEST) / N_COEFFICIENTS)

	sim_data.constants.kineticRateLimitFactorUpper = fluxLimit

	return dict(
		shortName = "upper_flux %.1f" % (fluxLimit),
		desc = "Rate limit upper coefficient set to %.2f." % (fluxLimit)
		), sim_data
