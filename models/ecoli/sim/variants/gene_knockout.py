"""
Knockout expression of a gene

Modifies:
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob

Expected variant indices (depends on length of sim_data.process.transcription.rna_data):
	0: control
	1-4692: gene index to knockout
"""

import numpy as np


CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def gene_knockout(sim_data, index):
	transcription = sim_data.process.transcription
	transcription_regulation = sim_data.process.transcription_regulation

	nGenes = len(transcription.rna_data)
	nConditions = nGenes + 1

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	geneIndex = (index - 1) % nConditions

	factor = 0  # Knockout expression
	recruitment_mask = np.array([i == geneIndex
		for i in transcription_regulation.delta_prob['deltaI']])
	for synth_prob in transcription.rna_synth_prob.values():
		synth_prob[geneIndex] *= factor
	for exp in transcription.rna_expression.values():
		exp[geneIndex] *= factor
	transcription.exp_free[geneIndex] *= factor
	transcription.exp_ppgpp[geneIndex] *= factor
	transcription_regulation.basal_prob[geneIndex] *= factor
	transcription_regulation.delta_prob['deltaV'][recruitment_mask] *= factor

	# Renormalize parameters
	for synth_prob in transcription.rna_synth_prob.values():
		synth_prob /= synth_prob.sum()
	for exp in transcription.rna_expression.values():
		exp /= exp.sum()
	transcription.exp_free /= transcription.exp_free.sum()
	transcription.exp_ppgpp /= transcription.exp_ppgpp.sum()

	geneID = transcription.rna_data["id"][geneIndex]

	return dict(
		shortName = "{}_KO".format(geneID),
		desc = "Complete knockout of {}.".format(geneID)
		), sim_data
