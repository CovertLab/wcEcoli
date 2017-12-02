FACTORS_SYNTH_PROB = [0.1, 1., 10.]
FACTORS_TRANSL_EFF = [0.1, 1., 10.]

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def meneParamsTotalIndices(sim_data):
	return None


def meneParams(sim_data, index):
	if index == 4:
		return CONTROL_OUTPUT, sim_data

	# Store factor for changing synthesis prob
	factor_synthProb = FACTORS_SYNTH_PROB[index / 3]
	sim_data.process.transcription.meneFactor = factor_synthProb

	# Store factor for changing translation efficiency
	factor_translEff = FACTORS_TRANSL_EFF[index % 3]
	sim_data.process.translation.meneFactor = factor_translEff

	return dict(
		shortName = "{}_meneParams".format(index),
		desc = "Simulation with menE synthesis probability increased by factor {}, translation efficiency increased by factor {}.".format(FACTORS_SYNTH_PROB[index / 3], FACTORS_TRANSL_EFF[index % 3])
		), sim_data
