"""
Tests combinations of kcat estimate quantile choice and a scalar multiplier
applied to all selected kcat values. Use this to sweep how the tightness and
magnitude of kcat-based kinetic targets affect growth.

Index 0 is a wildtype control (no kcat bounds applied). Indices 1-20 are
quantile x multiplier combinations (KCAT_QUANTILES x KCAT_MULTIPLIERS,
row-major), where index i >= 1 maps to:
	quantile   = KCAT_QUANTILES[(i - 1) // len(KCAT_MULTIPLIERS)]
	multiplier = KCAT_MULTIPLIERS[(i - 1)  % len(KCAT_MULTIPLIERS)]

When kcat bounds are active, sim_data.process.metabolism is modified:
	sim_data.process.metabolism.kcat_estimate_quantile  (str)
	sim_data.process.metabolism.kcat_estimate_multiplier  (float)
	sim_data.process.metabolism.selected_kcat_estimates  (dict)

Expected variant indices:
	 0: wildtype control (no kcat bounds)
	 1: p99 x 1.00
	 2: p99 x 0.95
	 3: p99 x 0.90
	 4: p99 x 0.85
	 5: p99 x 0.80
	 6: p95 x 1.00
	 7: p95 x 0.95
	 8: p95 x 0.90
	 9: p95 x 0.85
	10: p95 x 0.80
	11: p90 x 1.00
	12: p90 x 0.95
	13: p90 x 0.90
	14: p90 x 0.85
	15: p90 x 0.80
	16: median x 1.00
	17: median x 0.95
	18: median x 0.90
	19: median x 0.85
	20: median x 0.80
"""

KCAT_QUANTILES = ['p99', 'p95', 'p90', 'median']
KCAT_MULTIPLIERS = [1.0, 0.95, 0.9, 0.85, 0.8]


def kcat_estimate_scale(sim_data, index):
	if index == 0:
		return dict(
			shortName='wildtype_control',
			desc='No kcat bounds applied; matches wildtype behavior.',
		), sim_data

	i = index - 1
	n_multipliers = len(KCAT_MULTIPLIERS)
	quantile = KCAT_QUANTILES[i // n_multipliers]
	multiplier = KCAT_MULTIPLIERS[i % n_multipliers]

	base_estimates = sim_data.process.metabolism.kcat_estimates[quantile]

	sim_data.process.metabolism.kcat_estimate_quantile = quantile
	sim_data.process.metabolism.kcat_estimate_multiplier = multiplier
	sim_data.process.metabolism.selected_kcat_estimates = {
		key: value * multiplier for key, value in base_estimates.items()
	}

	multiplier_pct = int(round(multiplier * 100))
	return dict(
		shortName=f'{quantile}_x{multiplier_pct}pct',
		desc=(
			f'kcat estimates from {quantile} quantile,'
			f' scaled by {multiplier} ({multiplier_pct}%).'
		),
	), sim_data
