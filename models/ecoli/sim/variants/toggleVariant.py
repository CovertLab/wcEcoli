
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def toggleVariantTotalIndices(kb):
	nVariants = 3
	nConditions = nVariants + 1
	return nConditions


def toggleVariant(kb, index):
	# Toggles variants you specify here

	nConditions = toggleVariantTotalIndices(kb)

	kb.ppGppFeedback = False
	kb.turnOnGlucoseLimitation = False

	if index % nConditions == 0:
		return CONTROL_OUTPUT

	if index == 1:
		# Turn on transcription initiation feedback
		kb.ppGppFeedback = True

		return dict(
			shortName = "growthRateControl",
			desc = "Transcription initiation feedback by ppGpp. 1/2 glucose at 10 min."
			)

	elif index == 2:
		# Turn on glucose limit
		# Turn on transcription initiation feedback
		kb.turnOnGlucoseLimitation = True
		kb.ppGppFeedback = True

		return dict(
			shortName = "growthRateControl_glucoseLimit",
			desc = "Transcription initiation feedback by ppGpp. 1/2 glucose at 10 min."
			)

	elif index == 3:
		# Glucose limit
		kb.turnOnGlucoseLimitation = True

		return dict(
			shortName = "glucoseLimit",
			desc = "1/2 glucose at 10 min."
			)
