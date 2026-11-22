"""
Max kcat constraints + new gene internal shift.

Combines 1.00 * max kcat upper bounds (all generations) with phased new
gene induction from new_gene_internal_shift.

Expected variant indices:
	0:     control (new gene knockout) + max kcat constraints
	1-N:   new gene expression/translation combos + max kcat constraints,
	       N = (number of expression factors - 1) x number of translation
	       efficiencies in new_gene_internal_shift
	N + 1: wildtype control (no kcat constraints, no new gene modifications)

The wildtype control's index is derived from the grid size, so changing the
grid in new_gene_internal_shift cannot make it collide with a GFP design.
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	induce_new_genes,
	knockout_induced_new_gene_expression,
	get_new_gene_expression_factor_and_translation_efficiency,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
	CONTROL_OUTPUT,
)

# GFP indices span 0 to (len(EXPRESSION_FACTORS)-1) * SEPARATOR, where
# index 0 is control and 1..max_gfp_index cover expression x translation.
# The wildtype control comes right after them. It was a fixed 21, which fails
# the moment the grid grows past 20 designs; this module is imported for every
# simulation (the variant registry loads all variants), so that would have
# stopped every run, not just this variant's.
from models.ecoli.sim.variants.new_gene_internal_shift import (
	NEW_GENE_EXPRESSION_FACTORS as _EXP,
	NEW_GENE_TRANSLATION_EFFICIENCY_VALUES as _TRL,
)
_MAX_GFP_INDEX = (len(_EXP) - 1) * len(_TRL)
WILDTYPE_CONTROL_INDEX = _MAX_GFP_INDEX + 1


def new_gene_internal_shift_with_kcat(sim_data, index):
	# Wildtype control: no kcat, no GFP
	if index == WILDTYPE_CONTROL_INDEX:
		return dict(
			shortName='wildtype_no_kcat',
			desc='Wildtype control with no kcat constraints and no new gene modifications.',
		), sim_data

	# --- Apply max kcat constraints (all gens) ---
	metabolism = sim_data.process.metabolism
	metabolism.kcat_estimate_quantile = 'max'
	metabolism.kcat_estimate_multiplier = 1.0
	metabolism.selected_kcat_estimates = dict(metabolism.kcat_estimates['max'])

	# --- New gene logic (identical to new_gene_internal_shift) ---
	condition_index = index // 1000
	condition(sim_data, condition_index)

	index_remainder = index - condition_index * 1000
	(expression_list_index, trl_eff_list_index, expression_factor,
	 trl_eff_value) = get_new_gene_expression_factor_and_translation_efficiency(
		sim_data, index_remainder)

	setattr(sim_data, 'internal_shift_dict', {})
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(induce_new_genes, index_remainder)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(knockout_induced_new_gene_expression, index_remainder)]

	if index == 0:
		return dict(
			shortName='control_max_kcat',
			desc='Control (new gene knockout) with 1.0x max kcat constraints.',
		), sim_data

	if index == WILDTYPE_CONTROL_INDEX:
		return CONTROL_OUTPUT, sim_data

	return dict(
		shortName="{}_NGEXP_TRLEFF_INTERNAL_SHIFT."
			+ "{}_env_max_kcat".format(sim_data.ordered_conditions[condition_index]),
		desc="Changed expression of new genes by variant index {}.".format(
			expression_factor) + ". Set translation efficiency of new genes "
			"to {}".format(trl_eff_value)
			+ "Simulation of condition {}.".format(sim_data.ordered_conditions[condition_index])
			+ " With 1.0x max kcat constraints."
		), sim_data
