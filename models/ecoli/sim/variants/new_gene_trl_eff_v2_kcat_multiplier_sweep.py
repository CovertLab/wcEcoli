"""
Cross GFP burden (translation efficiency) with a kcat MULTIPLIER on the v2
estimators (gap and smoothed_max).  The new_gene_trl_eff_v2_estimator_sweep
showed metabolism is inert at face-value (1.0x) kcat; this sweep lowers the
bound below the estimate via a multiplier -- a proxy for "the true kcats are X x
below this estimate" -- to map where (if anywhere) GFP burden acquires a
metabolic component, and whether the estimator choice then matters.

Layout (61 indices: 1 wildtype + 2 estimators x 6 trl-eff x 5 multipliers):

  0:        Wildtype control (GFP knockout, no kcat bound)
  1 -- 30:  gap          x (trl-eff x multiplier)
  31 -- 60: smoothed_max x (trl-eff x multiplier)

Within each estimator block the cell index decodes as
(trl_idx, mult_idx) = divmod(rem, len(MULTIPLIERS)).  Reference edges per
estimator: multiplier = 1.0 reproduces the GFP-only response; trl_eff = 0.0
reproduces that estimator's pure multiplier sweep (no translation burden).

Reuses the GFP-induction helpers from new_gene_trl_eff_v2_estimator_sweep
(_induce, _knockout, EXPRESSION_FACTOR) and condition / induction gens from
new_gene_internal_shift, so it stays in lockstep with the established sweep.

Modifies (per non-wildtype variant):
	sim_data.process.metabolism.kcat_estimate_quantile
	sim_data.process.metabolism.kcat_estimate_multiplier
	sim_data.process.metabolism.selected_kcat_estimates
	sim_data.condition / external_state.current_timeline_id
	(after shift) new gene expression + translation efficiency
"""

from models.ecoli.sim.variants.new_gene_internal_shift import (
	condition,
	NEW_GENE_INDUCTION_GEN,
	NEW_GENE_KNOCKOUT_GEN,
)
from models.ecoli.sim.variants.new_gene_trl_eff_v2_estimator_sweep import (
	_induce,
	_knockout,
	EXPRESSION_FACTOR,
)


ESTIMATORS = ('v2_gap', 'v2_smoothed_max')  # ('v2_gap',) for the gap-only fallback
ESTIMATOR_SHORT = {'v2_gap': 'gap', 'v2_smoothed_max': 'smoothed_max'}
TRL_EFF_VALUES = [5.0, 4.0, 3.0, 2.0, 1.0, 0.0]   # 6; 0.0 = no-translation ref
MULTIPLIERS = [1.0, 0.8, 0.6, 0.5, 0.4]           # 5; 1.0 = no extra tightening
N_PER_ESTIMATOR = len(TRL_EFF_VALUES) * len(MULTIPLIERS)  # 30
N_VARIANTS = 1 + len(ESTIMATORS) * N_PER_ESTIMATOR        # 61


def is_wildtype(index: int) -> bool:
	"""Return True for the index-0 wildtype control (GFP knockout, no kcat)."""
	return index == 0


def variant_to_cell(index: int) -> tuple[int, int, int]:
	"""Return (estimator_idx, trl_idx, mult_idx) for a non-wildtype index.

	estimator_idx indexes ESTIMATORS; trl_idx indexes TRL_EFF_VALUES;
	mult_idx indexes MULTIPLIERS.
	"""
	if index < 1 or index >= N_VARIANTS:
		raise ValueError(
			f'new_gene_trl_eff_v2_kcat_multiplier_sweep: index {index} out of '
			f'range [1, {N_VARIANTS - 1}]; 0 is the wildtype control.')
	est_idx, rem = divmod(index - 1, N_PER_ESTIMATOR)
	trl_idx, mult_idx = divmod(rem, len(MULTIPLIERS))
	return est_idx, trl_idx, mult_idx


def new_gene_trl_eff_v2_kcat_multiplier_sweep(sim_data, index):
	"""
	Apply variant.  Index 0 is the wildtype control (GFP knockout, no kcat).
	Otherwise decode (estimator, trl_eff, multiplier): apply the estimator's
	kcat bound scaled by the multiplier, then induce GFP at the fixed expression
	factor and this variant's translation efficiency via the internal shift.
	"""
	# Always minimal media.
	condition(sim_data, 0)
	setattr(sim_data, 'internal_shift_dict', {})

	if is_wildtype(index):
		if NEW_GENE_INDUCTION_GEN != -1:
			sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
				(_knockout, 0)]
		return dict(
			shortName='wildtype_control',
			desc='Wildtype control (GFP knockout, no kcat bound).',
		), sim_data

	est_idx, trl_idx, mult_idx = variant_to_cell(index)
	estimator = ESTIMATORS[est_idx]
	trl_eff_value = TRL_EFF_VALUES[trl_idx]
	multiplier = MULTIPLIERS[mult_idx]

	# Apply the estimator's kcat bound scaled by the multiplier.
	metabolism = sim_data.process.metabolism
	assert estimator in metabolism.kcat_estimates, (
		f"sim_data.process.metabolism.kcat_estimates['{estimator}'] missing -- "
		f"Parca must populate this from reconstruction/ecoli/flat/"
		f"kcat_estimates/kcat_estimates_{estimator}.tsv.")
	base_estimates = metabolism.kcat_estimates[estimator]
	if not base_estimates:
		raise ValueError(
			f"metabolism.kcat_estimates['{estimator}'] is empty -- the v2 "
			f"estimator TSV is still a placeholder.")
	metabolism.kcat_estimate_quantile = estimator
	metabolism.kcat_estimate_multiplier = multiplier
	metabolism.selected_kcat_estimates = {
		key: value * multiplier for key, value in base_estimates.items()}

	# Induce GFP at NEW_GENE_INDUCTION_GEN with this variant's trl_eff.
	if NEW_GENE_INDUCTION_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_INDUCTION_GEN] = [
			(_induce, trl_eff_value)]
	if NEW_GENE_KNOCKOUT_GEN != -1:
		sim_data.internal_shift_dict[NEW_GENE_KNOCKOUT_GEN] = [(_knockout, 0)]

	pct = int(round(multiplier * 100))
	return dict(
		shortName=f'{ESTIMATOR_SHORT[estimator]}_trl{trl_eff_value}_x{pct:03d}',
		desc=(
			f'GFP exp 10^{EXPRESSION_FACTOR - 1}, translation efficiency '
			f'{trl_eff_value}, kcat {estimator} at {multiplier:.2f}x.'),
	), sim_data
