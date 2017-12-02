FACTORS_SYNTH_PROB = [0.1, 1., 10.]
FACTORS_FLUX = [0.1, 1., 10.]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def meneParamsTotalIndices(sim_data):
	return None


def meneParams(sim_data, index):
	if index == 0:
		return CONTROL_OUTPUT, sim_data

	factor_synthProb = None

	# Get index of menE gene in rnaData
	factor_synthProb = FACTORS_SYNTH_PROB[index - 1]
	sim_data.process.transcription.meneFactor = factor_synthProb

	return dict(
		shortName = "{}_meneParams".format(index),
		desc = "Simulation with menE synthesis probability increased by factor {}.".format(FACTORS_SYNTH_PROB[index - 1])
		), sim_data
