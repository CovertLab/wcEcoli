"""
Sweep the three v2 kcat estimators (gap, smoothed_max, drop_top_20) against a
uniform multiplier grid, with a wildtype control.

Estimator sources (loaded from reconstruction/ecoli/flat/kcat_estimates/):
	* v2_gap          -- top-K gap estimator (kcat_estimates_v2_gap.tsv)
	* v2_smoothed_max -- median-filter W=5 estimator
	                     (kcat_estimates_v2_smoothed_max.tsv)
	* v2_drop_top_20  -- drop the 20 largest samples, return the next
	                     (kcat_estimates_v2_drop_top_20.tsv)

Each estimator is swept against the same 10-point multiplier grid for direct
apples-to-apples comparison.

Variant index layout (31 total indices):
	  0       wildtype control (no kcat bounds applied)
	  1 -- 10 v2_gap          x [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
	 11 -- 20 v2_smoothed_max x same grid
	 21 -- 30 v2_drop_top_20  x same grid

Decoding (for index >= 1):
	estimator_idx, mult_idx = divmod(index - 1, N_PER_ESTIMATOR)

The variant function raises ValueError if the chosen estimator's TSV is
empty (i.e. placeholder still in place) so the sweep cannot accidentally
run with zero kcat constraints applied.
"""

ESTIMATOR_LABELS = ('v2_gap', 'v2_smoothed_max', 'v2_drop_top_20')
ESTIMATOR_SHORT_NAMES = ('gap', 'smoothed_max', 'drop_top_20')
MULTIPLIERS = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
N_PER_ESTIMATOR = len(MULTIPLIERS)
N_VARIANTS = 1 + len(ESTIMATOR_LABELS) * N_PER_ESTIMATOR  # 31


def is_wildtype(index: int) -> bool:
	return index == 0


def variant_to_estimator(index: int) -> tuple[int, int]:
	"""Return (estimator_idx, multiplier_idx) for a non-wildtype index.

	estimator_idx in [0, len(ESTIMATOR_LABELS)) selects the source TSV;
	multiplier_idx in [0, N_PER_ESTIMATOR) indexes into MULTIPLIERS.
	"""
	if index < 1 or index >= N_VARIANTS:
		raise ValueError(
			f'kcat_v2_estimator_sweep: index {index} out of range '
			f'[1, {N_VARIANTS - 1}]; 0 is the wildtype control.')
	return divmod(index - 1, N_PER_ESTIMATOR)


def kcat_v2_estimator_sweep(sim_data, index):
	if index == 0:
		return dict(
			shortName='wildtype_control',
			desc='No kcat bounds applied; matches wildtype behavior.',
		), sim_data

	estimator_idx, mult_idx = variant_to_estimator(index)
	label = ESTIMATOR_LABELS[estimator_idx]
	short_name = ESTIMATOR_SHORT_NAMES[estimator_idx]
	multiplier = MULTIPLIERS[mult_idx]
	multiplier_pct = int(round(multiplier * 100))

	metabolism = sim_data.process.metabolism
	base_estimates = metabolism.kcat_estimates.get(label)
	if not base_estimates:
		raise ValueError(
			f"metabolism.kcat_estimates['{label}'] is empty -- the v2 "
			f"estimator TSV at reconstruction/ecoli/flat/kcat_estimates/"
			f"kcat_estimates_{label}.tsv is still a placeholder.  Generate "
			f"a real TSV with models/ecoli/analysis/cohort/"
			f"kcat_estimates_v2_to_tsv.py and drop it in before running "
			f"this sweep.")

	metabolism.kcat_estimate_quantile = label
	metabolism.kcat_estimate_multiplier = multiplier
	metabolism.selected_kcat_estimates = {
		key: value * multiplier for key, value in base_estimates.items()
	}

	return dict(
		shortName=f'{short_name}_x{multiplier_pct:03d}',
		desc=(
			f'kcat estimates from {label} scaled by {multiplier:.2f} '
			f'({multiplier_pct}%).'
		),
	), sim_data
