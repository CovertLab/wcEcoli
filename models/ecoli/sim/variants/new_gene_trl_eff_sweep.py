"""
Sweep new-gene (GFP) translation efficiencies at a fixed expression factor
(8.5) across four kcat-constraint regimes. Enables direct comparison of
doubling-time impacts as the kcat bound is progressively tightened at matched
GFP expression levels.

Layout (48 indices total, 4 categories x 12 variants each):

Category 0 — no kcat constraint (indices 0-11):
  0:       Control (GFP knockout, exp factor = 0)
  1-11:    exp=8.5, trl_eff = TRL_EFF_VALUES

Category 1 — kcat = 1.0 x max (indices 12-23):
  12:      Control (GFP knockout)
  13-23:   exp=8.5, trl_eff = TRL_EFF_VALUES

Category 2 — kcat = 0.8 x max (indices 24-35):
  24:      Control (GFP knockout)
  25-35:   exp=8.5, trl_eff = TRL_EFF_VALUES

Category 3 — kcat = 0.6 x max (indices 36-47):
  36:      Control (GFP knockout)
  37-47:   exp=8.5, trl_eff = TRL_EFF_VALUES

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

Modifies (categories 1-3 only, i.e. when a kcat multiplier is applied):
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
	5.0, 4.5, 4.0, 3.5, 3.0,
	2.5, 2.0, 1.5, 1.0, 0.5,
	0.0]
EXPRESSION_FACTOR = 8.5
KCAT_MULTIPLIERS = [None, 1.0, 0.8, 0.6]  # None => no kcat bound
N_TRL_EFF = len(TRL_EFF_VALUES)  # 11
N_PER_CATEGORY = N_TRL_EFF + 1  # 12 (1 control + 11 expression variants)
N_VARIANTS = N_PER_CATEGORY * len(KCAT_MULTIPLIERS)  # 48


def variant_to_category(index: int) -> tuple[int, int]:
	"""Return (cat_idx, local_idx) where cat_idx in 0..3 and local_idx in 0..11."""
	return divmod(index, N_PER_CATEGORY)


def is_control(index: int) -> bool:
	"""Return True if this variant is the control (GFP knockout) for its category.

	Controls are local_idx == 0 within each category (global indices
	{0, 12, 24, 36}).
	"""
	_, local_idx = variant_to_category(index)
	return local_idx == 0


def category_label(cat_idx: int) -> str:
	"""Return a short human label for the category, e.g. 'no_kcat', 'kcat_1.0x'."""
	mult = KCAT_MULTIPLIERS[cat_idx]
	return 'no_kcat' if mult is None else f'kcat_{mult:.1f}x'


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
	Apply variant. Determines the (category, local index) pair for this
	variant, applies the category's kcat multiplier (if any), then sets up
	the internal shift for new gene induction at the appropriate expression
	factor and translation efficiency.
	"""
	cat_idx, local_idx = variant_to_category(index)
	multiplier = KCAT_MULTIPLIERS[cat_idx]

	# Apply kcat constraint if this category has one
	if multiplier is not None:
		assert 'max' in sim_data.process.metabolism.kcat_estimates, \
			"sim_data.process.metabolism.kcat_estimates['max'] missing" \
			" — Parca must populate this"
		metabolism = sim_data.process.metabolism
		metabolism.kcat_estimate_quantile = 'max'
		metabolism.kcat_estimate_multiplier = multiplier
		metabolism.selected_kcat_estimates = {
			key: value * multiplier
			for key, value in metabolism.kcat_estimates['max'].items()
		}

	# Always minimal media
	condition(sim_data, 0)

	# Initialize internal shift dictionary
	setattr(sim_data, 'internal_shift_dict', {})

	# Add the new gene induction to the internal_shift instructions
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, local_idx)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(_induce, 0)]  # knockout = set expression to 0

	# Build description
	cat = category_label(cat_idx)
	if local_idx == 0:
		shortName = f'control_{cat}'
		desc = f'Control (new gene knockout), {cat}.'
	else:
		trl_eff_value = TRL_EFF_VALUES[local_idx - 1]
		shortName = (
			f'trl_eff_{trl_eff_value}_exp_{EXPRESSION_FACTOR}_{cat}')
		desc = (
			f'Expression factor 10^{EXPRESSION_FACTOR - 1}, '
			f'translation efficiency {trl_eff_value}, {cat}.')
	return dict(shortName=shortName, desc=desc), sim_data
