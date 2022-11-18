"""
Variant to compare the impacts of knocking out a subset of the rRNA operons.

Modifies:
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob

Expected variant indices:
	0: control
	1-7: knockout of one rRNA operon
	8-28: knockout of two rRNA operons
	29-63: knockout of three rRNA operons
	64-98: knockout of four rRNA operons
	99-119: knockout of five rRNA operons
	120-126: knockout of six rRNA operons
	127: knockout of all seven rRNA operons
"""

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rrna_operon_knockout(sim_data, index):
	if index > 0:
		def int_to_bool_array(i):
			return np.array([bool(i & (1 << n)) for n in range(7)])

		# Build arrays for every combination of rRNAs to knock out
		rrna_ko_array = np.array([
			int_to_bool_array(i) for i in range(1, 128)
			])

		# Sort by number of rRNAs being knocked out
		rrna_ko_array = rrna_ko_array[np.argsort(rrna_ko_array.sum(axis=1)), :]

		# Get genes to knock out for this variant
		rRNA_operons_to_ko = [
			sim_data.molecule_groups.rRNA_operons[i] for i in
			np.where(rrna_ko_array[index - 1, :])[0]]
		genes_to_ko = []
		for rRNA_operon in rRNA_operons_to_ko:
			genes_to_ko.extend(sim_data.molecule_groups.__dict__[rRNA_operon])

		# Set expression levels of all RNAs encoding for the genes to zero
		rna_indexes_to_ko = np.where(
			sim_data.process.transcription.cistron_tu_mapping_matrix.T.dot(
				np.isin(
					sim_data.process.transcription.cistron_data['gene_id'],
					genes_to_ko)
				))[0]
		sim_data.adjust_final_expression(
			rna_indexes_to_ko, np.zeros_like(rna_indexes_to_ko))

		return dict(
			shortName = f"{'_'.join(rRNA_operons_to_ko)}_rrna_knockout",
			desc = f"Simulation with the {', '.join(rRNA_operons_to_ko)} rRNA operons knocked out"
			), sim_data

	else:
		return CONTROL_OUTPUT, sim_data
