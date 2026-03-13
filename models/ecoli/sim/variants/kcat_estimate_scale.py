"""
Tests combinations of kcat estimate quantile choice and a scalar multiplier
applied to all selected kcat values. Use this to sweep how the tightness and
magnitude of kcat-based kinetic targets affect growth.

The active kcat estimates are written to sim_data.process.metabolism so the
metabolism process can apply them as kinetic constraints.

Modifies:
	sim_data.process.metabolism.kcat_estimate_quantile  (str)
	sim_data.process.metabolism.kcat_estimate_multiplier  (float)
	sim_data.process.metabolism.selected_kcat_estimates  (dict)

Expected variant indices (KCAT_QUANTILES x KCAT_MULTIPLIERS, row-major):
	Each index i maps to:
		quantile  = KCAT_QUANTILES[i // len(KCAT_MULTIPLIERS)]
		multiplier = KCAT_MULTIPLIERS[i  % len(KCAT_MULTIPLIERS)]

	 0: p99 x 1.00
	 1: p99 x 0.95
	 2: p99 x 0.90
	 3: p99 x 0.85
	 4: p99 x 0.80
	 5: p95 x 1.00
	 6: p95 x 0.95
	 7: p95 x 0.90
	 8: p95 x 0.85
	 9: p95 x 0.80
	10: p90 x 1.00
	11: p90 x 0.95
	12: p90 x 0.90
	13: p90 x 0.85
	14: p90 x 0.80
	15: median x 1.00
	16: median x 0.95
	17: median x 0.90
	18: median x 0.85
	19: median x 0.80
"""

KCAT_QUANTILES = ['p99', 'p95', 'p90', 'median']
KCAT_MULTIPLIERS = [1.0, 0.95, 0.9, 0.85, 0.8]


def kcat_estimate_scale(sim_data, index):
	n_multipliers = len(KCAT_MULTIPLIERS)
	quantile = KCAT_QUANTILES[index // n_multipliers]
	multiplier = KCAT_MULTIPLIERS[index % n_multipliers]

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
