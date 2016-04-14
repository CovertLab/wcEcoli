RELAXATION_VALUE = 100000000

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def kineticsIndivCoeffTotalIndices(sim_data):
	# Only adjusts reactions which have a non-one kinetics coefficient
	nUniqueReactions = len(set([x["reactionID"] for x in sim_data.process.metabolism.reactionRateInfo.values() if x["constraintMultiple"] > 1]))
	nConditions = nUniqueReactions + 1
	return nConditions


def kineticsIndivCoeff(sim_data, index):
	# Sets the coefficient of the upper limit of flux used for enzymeKinetics estimates

	nConditions = kineticsIndivCoeffTotalIndices(sim_data)

	reactionID = sorted(list(set([x["reactionID"] for x in sim_data.process.metabolism.reactionRateInfo.values()])))[index]

	relaxAll = False
	shortName = "Removed constraint on %s" % (reactionID),
	desc = "Rate limit on constraint %s set to %f." % (reactionID, RELAXATION_VALUE)

	# Relax all constraints in the control simulation
	if index % nConditions == 0:
		relaxAll = True
		shortName = "Control simulation"
		desc = "Control simulation"

	for reactionInfo in sim_data.process.metabolism.reactionRateInfo.values():
		if reactionInfo["reactionID"] == reactionID or relaxAll:
			reactionInfo["constraintMultiple"] = RELAXATION_VALUE

	return dict(
		shortName = shortName,
		desc = desc
		), sim_data
