"""
Sweep new-gene (GFP) translation efficiencies at a fixed expression factor
(8.5) across three kcat regimes that differ only in the v2 estimator used,
all applied at multiplier 1.0.  Unlike new_gene_trl_eff_sweep.py (which holds
the estimator at 'max' and sweeps the kcat multiplier), here the multiplier is
fixed and the *estimator* changes: the metabolic reduction comes entirely from
reduced enzyme counts as GFP competes for ribosomes, not from artificial kcat
scaling.

Layout (40 indices total: 1 wildtype + 3 blocks x 13 trl-effs):

  0:       Wildtype control (GFP knockout, no kcat bound)

Block 0 — gfp only, no kcat bound (indices 1-13):
  1-13:    exp=8.5, trl_eff = TRL_EFF_VALUES

Block 1 — gfp + v2_gap kcat bound at 1.0x (indices 14-26):
  14-26:   exp=8.5, trl_eff = TRL_EFF_VALUES

Block 2 — gfp + v2_smoothed_max kcat bound at 1.0x (indices 27-39):
  27-39:   exp=8.5, trl_eff = TRL_EFF_VALUES

trl_eff = 0.0 is included in every block as the within-block "transcription
burden, no translation" control: GFP is transcribed at the fixed expression
factor but not translated, so it imposes no ribosome competition.  Index 0
(GFP knockout) is the absolute baseline.

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Modifies (blocks 1-2 only, i.e. when a kcat estimator is applied):
	sim_data.process.metabolism.kcat_estimate_quantile
	sim_data.process.metabolism.kcat_estimate_multiplier
	sim_data.process.metabolism.selected_kcat_estimates

Modifies (after shift, all GFP variants):
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
)

TRL_EFF_VALUES = [
	5.0, 4.75, 4.5, 4.25, 4.0,
	3.5, 3.0, 2.5, 2.0, 1.5,
	1.0, 0.5, 0.0]
EXPRESSION_FACTOR = 8.5
# Per-block kcat estimator source (None => no kcat bound).  All applied at
# multiplier 1.0; the v2_gap / v2_smoothed_max TSVs are loaded by Parca into
# metabolism.kcat_estimates.
BLOCK_ESTIMATORS = (None, 'v2_gap', 'v2_smoothed_max')
BLOCK_SHORT_NAMES = ('gfp_only', 'gap', 'smoothed_max')
KCAT_MULTIPLIER = 1.0
N_PER_BLOCK = len(TRL_EFF_VALUES)  # 13
N_VARIANTS = 1 + len(BLOCK_ESTIMATORS) * N_PER_BLOCK  # 40


def is_wildtype(index: int) -> bool:
	"""Return True for the index-0 wildtype control (GFP knockout, no kcat)."""
	return index == 0


def variant_to_block(index: int) -> tuple[int, int]:
	"""Return (block_idx, te_idx) for a non-wildtype index.

	block_idx in [0, len(BLOCK_ESTIMATORS)) selects the kcat regime;
	te_idx in [0, N_PER_BLOCK) indexes into TRL_EFF_VALUES.
	"""
	if index < 1 or index >= N_VARIANTS:
		raise ValueError(
			f'new_gene_trl_eff_v2_estimator_sweep: index {index} out of range '
			f'[1, {N_VARIANTS - 1}]; 0 is the wildtype control.')
	return divmod(index - 1, N_PER_BLOCK)


def block_label(block_idx: int) -> str:
	"""Return a short human label for the block, e.g. 'gfp_only', 'gap'."""
	return BLOCK_SHORT_NAMES[block_idx]


def _induce(sim_data, trl_eff_value):
	"""Induce new genes at the fixed EXPRESSION_FACTOR with the given
	translation efficiency.

	trl_eff_value == 0.0 still expresses (transcribes) GFP but blocks its
	translation -- the within-block "transcription burden only" control.
	"""
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	expression_factor = 10 ** (EXPRESSION_FACTOR - 1)
	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		monomer_index = new_gene_monomer_indices[i]

		sim_data.adjust_new_gene_final_expression(
			[gene_index], [expression_factor])
		sim_data.process.translation.translation_efficiencies_by_monomer[
			monomer_index] = trl_eff_value


def _knockout(sim_data, _ignored):
	"""Knock out new gene expression entirely (the index-0 wildtype control)."""
	new_gene_mRNA_ids, new_gene_indices, new_gene_monomer_ids, \
		new_gene_monomer_indices = determine_new_gene_ids_and_indices(
		sim_data)

	for i in range(len(new_gene_indices)):
		gene_index = new_gene_indices[i]
		sim_data.adjust_new_gene_final_expression([gene_index], [0])


def new_gene_trl_eff_v2_estimator_sweep(sim_data, index):
	"""
	Apply variant.  Index 0 is the wildtype control (GFP knockout, no kcat
	bound).  Otherwise decode (block, trl_eff): apply the block's v2 kcat
	estimator at multiplier 1.0 (if any), then induce GFP at the fixed
	expression factor and this variant's translation efficiency via the
	internal shift at NEW_GENE_INDUCTION_GEN.
	"""
	# Always minimal media.
	condition(sim_data, 0)

	# Initialize internal shift dictionary.
	setattr(sim_data, 'internal_shift_dict', {})

	if is_wildtype(index):
		# Wildtype control: knock out GFP, apply no kcat bound.
		if NEW_GENE_INDUCTION_GEN != -1:
			sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
				(_knockout, 0)]
		return dict(
			shortName='wildtype_control',
			desc='Wildtype control (GFP knockout, no kcat bound).',
		), sim_data

	block_idx, te_idx = variant_to_block(index)
	estimator = BLOCK_ESTIMATORS[block_idx]
	trl_eff_value = TRL_EFF_VALUES[te_idx]

	# Apply the block's kcat estimator at multiplier 1.0, if it has one.
	if estimator is not None:
		metabolism = sim_data.process.metabolism
		assert estimator in metabolism.kcat_estimates, (
			f"sim_data.process.metabolism.kcat_estimates['{estimator}'] "
			f"missing -- Parca must populate this from "
			f"reconstruction/ecoli/flat/kcat_estimates/"
			f"kcat_estimates_{estimator}.tsv.")
		base_estimates = metabolism.kcat_estimates[estimator]
		if not base_estimates:
			raise ValueError(
				f"metabolism.kcat_estimates['{estimator}'] is empty -- the v2 "
				f"estimator TSV is still a placeholder.")
		metabolism.kcat_estimate_quantile = estimator
		metabolism.kcat_estimate_multiplier = KCAT_MULTIPLIER
		metabolism.selected_kcat_estimates = {
			key: value * KCAT_MULTIPLIER
			for key, value in base_estimates.items()
		}

	# Induce GFP at NEW_GENE_INDUCTION_GEN with this variant's trl_eff; if a
	# knockout gen is configured, knock GFP back out there.
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, trl_eff_value)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [
			(_knockout, 0)]

	short_block = BLOCK_SHORT_NAMES[block_idx]
	shortName = f'{short_block}_trl_eff_{trl_eff_value}'
	if estimator is None:
		desc = (
			f'GFP exp 10^{EXPRESSION_FACTOR - 1}, translation efficiency '
			f'{trl_eff_value}, no kcat bound.')
	else:
		desc = (
			f'GFP exp 10^{EXPRESSION_FACTOR - 1}, translation efficiency '
			f'{trl_eff_value}, kcat {estimator} at '
			f'{KCAT_MULTIPLIER:.1f}x.')
	return dict(shortName=shortName, desc=desc), sim_data
