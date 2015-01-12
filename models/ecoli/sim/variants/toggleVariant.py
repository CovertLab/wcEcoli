
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def toggleVariantTotalIndices(kb):
	nVariants = 1
	nConditions = nVariants + 1
	return nConditions


def geneKnockout(kb, index):
	# Toggles variants you specify here

	nConditions = toggleVariantTotalIndices(kb)

	if index % nConditions == 0:
		return CONTROL_OUTPUT

	kb.turnOnSynthetaseContraints = False
	kb.turnOnGlucoseLimitation = False

	if index == 1:
		# Turn on synthetase constraints
		kb.turnOnSynthetaseContraints = True

		return dict(
			shortName = "machineryLimitNoGlucoseLimit",
			desc = "Elongation limited by machinery or AA. No glucose limitation."
			)
	elif index == 2:
		# Turn on synthetase constraints and glucose limit
		kb.turnOnSynthetaseContraints = True
		kb.turnOnGlucoseLimitation = True

		return dict(
			shortName = "machineryLimitGlucoseLimit",
			desc = "Elongation limited by machinery or AA. 1/2 Glucose limitation at 10 min."
			)

	elif index == 3:
		# Glucose limit
		kb.turnOnGlucoseLimitation = True

		return dict(
			shortName = "noMachineryLimitGlucoseLimit",
			desc = "Elongation limited by AA. 1/2 Glucose limitation at 10 min."
			)
