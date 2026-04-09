"""
Sweep translation efficiencies at a fixed expression factor (8.5), with half
of the indices running without kcat constraints and half running with kcat
constraints. Enables direct comparison of doubling time impacts when kcat
constraints are active vs inactive at matching GFP expression levels.

Variant index layout (44 indices total):

No-kcat half (indices 0-21):
  0:     Control — exp=0 (GFP knockout), no kcat constraints
  1-21:  exp=8.5, trl_eff = [5.00, 4.75, ..., 0.25, 0], no kcat constraints

Kcat half (indices 22-43):
  22:    Control — exp=0 (GFP knockout), WITH max kcat constraints
  23-43: exp=8.5, trl_eff = [5.00, 4.75, ..., 0.25, 0], WITH max kcat constraints

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

Modifies (kcat half only):
	sim_data.process.metabolism.kcat_estimate_quantile
	sim_data.process.metabolism.kcat_estimate_multiplier
	sim_data.process.metabolism.selected_kcat_estimates
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	determine_new_gene_ids_and_indices,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
	CONTROL_OUTPUT,
)

TRL_EFF_VALUES = [
	5.00, 4.75, 4.50, 4.25, 4.00,
	3.75, 3.50, 3.25, 3.00, 2.75,
	2.50, 2.25, 2.00, 1.75, 1.50,
	1.25, 1.00, 0.75, 0.50, 0.25,
	0]
EXPRESSION_FACTOR = 8.5
N_TRL_EFF = len(TRL_EFF_VALUES)  # 21
KCAT_HALF_START = N_TRL_EFF + 1  # 22


def _induce(sim_data, local_index):
	"""
	Induce new genes using this variant's EXPRESSION_FACTOR and TRL_EFF_VALUES.
	For local_index == 0 (control), expression is set to 0 (knockout).
	For local_index >= 1, expression is set to 10^(EXPRESSION_FACTOR-1) and
	translation efficiency is set to TRL_EFF_VALUES[local_index - 1].
	"""
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	if local_index == 0:
		# Control: knockout new gene expression
		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			sim_data.adjust_new_gene_final_expression([gene_index], [0])
	else:
		expression_factor = 10 ** (EXPRESSION_FACTOR - 1)
		trl_eff_value = TRL_EFF_VALUES[local_index - 1]

		for i in range(len(new_gene_indices)):
			gene_index = new_gene_indices[i]
			monomer_index = new_gene_monomer_indices[i]

			sim_data.adjust_new_gene_final_expression(
				[gene_index], [expression_factor])
			sim_data.process.translation.translation_efficiencies_by_monomer[
				monomer_index] = trl_eff_value


def new_gene_trl_eff_sweep(sim_data, index):
	"""
	Apply variant. Determines whether this index belongs to the no-kcat half
	or the kcat half, applies kcat constraints if needed, then sets up
	internal shift for new gene induction at the appropriate expression
	factor and translation efficiency.
	"""
	is_kcat = (index >= KCAT_HALF_START)
	local_index = index - KCAT_HALF_START if is_kcat else index

	# Apply kcat constraints for kcat half
	if is_kcat:
		metabolism = sim_data.process.metabolism
		metabolism.kcat_estimate_quantile = 'max'
		metabolism.kcat_estimate_multiplier = 1.0
		metabolism.selected_kcat_estimates = dict(
			metabolism.kcat_estimates['max'])

	# Always minimal media
	condition(sim_data, 0)

	# Initialize internal shift dictionary
	setattr(sim_data, 'internal_shift_dict', {})

	# Add the new gene induction to the internal_shift instructions
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, local_index)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(_induce, 0)]  # knockout = set expression to 0

	# Build description
	kcat_str = "WITH" if is_kcat else "WITHOUT"
	if local_index == 0:
		return dict(
			shortName=f"control_{'kcat' if is_kcat else 'no_kcat'}",
			desc=f"Control (new gene knockout) {kcat_str} max kcat constraints.",
		), sim_data

	trl_eff_value = TRL_EFF_VALUES[local_index - 1]
	return dict(
		shortName=f"trl_eff_{trl_eff_value}_exp_{EXPRESSION_FACTOR}"
			f"_{'kcat' if is_kcat else 'no_kcat'}",
		desc=f"Expression factor 10^{EXPRESSION_FACTOR - 1}, "
			f"translation efficiency {trl_eff_value}, "
			f"{kcat_str} max kcat constraints.",
	), sim_data
