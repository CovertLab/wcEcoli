"""
Compare the impacts of increasing expression level of new genes.

Modifies:
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob

Expected variant indices (int):
	0: control (knockout new gene expression)
	x > 0: multiply new gene expression by a factor of 10^(x-1)
"""

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def new_gene_expression(sim_data, index):
	if index == 0:
		factor = 0
	else:
		factor = 10**(index - 1)

	rna_data = sim_data.process.transcription.rna_data

	rna_ids = rna_data['id']
	new_gene_indices = [i for i in range(len(rna_ids)) if rna_ids[i].startswith('NG')]

	for geneIndex in new_gene_indices:
		sim_data.adjust_final_expression([geneIndex], [factor])
		geneID = rna_data["id"][geneIndex]

	if index == 0:
		return CONTROL_OUTPUT, sim_data

	return dict(
		shortName = "{}_NGEXP",
		desc = "Changed expression of new genes by a factor of {}.".format(factor)
		), sim_data
