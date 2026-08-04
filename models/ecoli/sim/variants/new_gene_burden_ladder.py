"""
Sweep new-gene (GFP) translation efficiency at a single fixed expression
factor to produce one clean burden axis.

This variant is deliberately minimal: it is reused unchanged across every
insertion-position and copy-number batch so that results from different
genomes are directly comparable. The only thing that varies between those
batches is NEW_GENES (which flat-file directory supplies the construct), not
the variant definition.

Layout (7 indices total):

  0:     Control (GFP knockout, expression factor = 0)
  1-6:   exp = EXPRESSION_FACTOR, trl_eff = TRL_EFF_VALUES[index - 1]

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Modifies (after shift):
	sim_data.process.transcription.rna_synth_prob
	sim_data.process.transcription.rna_expression
	sim_data.process.transcription.exp_free
	sim_data.process.transcription.exp_ppgpp
	sim_data.process.transcription.attenuation_basal_prob_adjustments
	sim_data.process.transcription_regulation.basal_prob
	sim_data.process.transcription_regulation.delta_prob
	sim_data.process.translation.translation_efficiencies_by_monomer
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	determine_new_gene_ids_and_indices,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
	CONTROL_OUTPUT,
)

# Ascending burden. Index 0 of this list is variant index 1.
TRL_EFF_VALUES = [0.1, 0.5, 1.0, 2.0, 3.5, 5.0]

# Per the convention in new_gene_internal_shift, an expression factor of x
# multiplies new gene expression by 10^(x - 1). 8.5 matches the ladder the
# PLAN's burden-range and copy-number predictions are calibrated against.
EXPRESSION_FACTOR = 8.5

N_TRL_EFF = len(TRL_EFF_VALUES)  # 6
N_VARIANTS = N_TRL_EFF + 1  # 7 (1 control + 6 burden levels)


def is_control(index):
	"""Return True if this variant is the control (GFP knockout)."""
	return index == 0


def _induce(sim_data, index):
	"""
	Induce new genes at this variant's expression factor and translation
	efficiency.

	For index == 0 (control), expression is set to 0 (knockout). For
	index >= 1, expression is set to 10^(EXPRESSION_FACTOR - 1) and
	translation efficiency is set to TRL_EFF_VALUES[index - 1].

	Every new gene index returned by determine_new_gene_ids_and_indices is
	given the same settings, so tandem copies at one locus are all induced
	together with identical per-copy expression.
	"""
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	if index == 0:
		# Control: knock out new gene expression
		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			sim_data.adjust_new_gene_final_expression([gene_index], [0])
	else:
		expression_factor = 10 ** (EXPRESSION_FACTOR - 1)
		trl_eff_value = TRL_EFF_VALUES[index - 1]

		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			monomer_index = new_gene_monomer_indices[i]

			sim_data.adjust_new_gene_final_expression(
				[gene_index], [expression_factor])
			sim_data.process.translation.translation_efficiencies_by_monomer[
				monomer_index] = trl_eff_value


def new_gene_burden_ladder(sim_data, index):
	"""
	Apply variant. Sets minimal media, then schedules new gene induction at
	this variant's expression factor and translation efficiency.
	"""
	assert 0 <= index < N_VARIANTS, (
		f"new_gene_burden_ladder index must be in [0, {N_VARIANTS - 1}],"
		f" got {index}")

	# Always minimal media
	condition(sim_data, 0)

	# Initialize internal shift dictionary
	setattr(sim_data, 'internal_shift_dict', {})

	# Add the new gene induction to the internal_shift instructions
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, index)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(_induce, 0)]  # knockout = set expression to 0

	# Build description
	if is_control(index):
		return dict(CONTROL_OUTPUT), sim_data

	trl_eff_value = TRL_EFF_VALUES[index - 1]
	shortName = f'trl_eff_{trl_eff_value}_exp_{EXPRESSION_FACTOR}'
	desc = (
		f'Expression factor 10^{EXPRESSION_FACTOR - 1}, '
		f'translation efficiency {trl_eff_value}.')
	return dict(shortName=shortName, desc=desc), sim_data
