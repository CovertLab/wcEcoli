"""
Sweeps a scalar multiplier applied to the smoothed_max_buffered kcat estimates
for below-line essential reactions.  Use this to test how the tightness of
kcat-based flux upper bounds affects growth relative to an unconstrained
wildtype control.

Index 0 is a wildtype control (no kcat bounds applied).  Indices 1-21 use the
'smoothed_max_buffered' kcat estimates scaled by a multiplier from 1.0 down to
0.0 in increments of 0.05.  The buffered variant applies a 1.1x multiplier to
low-variance reactions (where smoothed_max/median < 1.05) to provide headroom
for stochastic variation.

When kcat bounds are active, sim_data.process.metabolism is modified:
	sim_data.process.metabolism.kcat_estimate_quantile   (str)
	sim_data.process.metabolism.kcat_estimate_multiplier (float)
	sim_data.process.metabolism.selected_kcat_estimates  (dict)

Expected variant indices:
	 0: wildtype control (no kcat bounds)
	 1: smoothed_max x 100%
	 2: smoothed_max x  95%
	 3: smoothed_max x  90%
	 4: smoothed_max x  85%
	 5: smoothed_max x  80%
	 6: smoothed_max x  75%
	 7: smoothed_max x  70%
	 8: smoothed_max x  65%
	 9: smoothed_max x  60%
	10: smoothed_max x  55%
	11: smoothed_max x  50%
	12: smoothed_max x  45%
	13: smoothed_max x  40%
	14: smoothed_max x  35%
	15: smoothed_max x  30%
	16: smoothed_max x  25%
	17: smoothed_max x  20%
	18: smoothed_max x  15%
	19: smoothed_max x  10%
	20: smoothed_max x   5%
	21: smoothed_max x   0%
"""

import numpy as np

KCAT_QUANTILE = 'smoothed_max_buffered'
KCAT_MULTIPLIERS = np.round(np.arange(1.0, -0.001, -0.05), 2)


def kcat_estimate_scale(sim_data, index):
	if index == 0:
		return dict(
			shortName='wildtype_control',
			desc='No kcat bounds applied; matches wildtype behavior.',
		), sim_data

	multiplier = KCAT_MULTIPLIERS[index - 1]
	base_estimates = sim_data.process.metabolism.kcat_estimates[KCAT_QUANTILE]

	sim_data.process.metabolism.kcat_estimate_quantile = KCAT_QUANTILE
	sim_data.process.metabolism.kcat_estimate_multiplier = multiplier
	sim_data.process.metabolism.selected_kcat_estimates = {
		key: value * multiplier for key, value in base_estimates.items()
	}

	multiplier_pct = int(round(multiplier * 100))
	return dict(
		shortName=f'{KCAT_QUANTILE}_x{multiplier_pct}pct',
		desc=(
			f'kcat estimates from {KCAT_QUANTILE} quantile,'
			f' scaled by {multiplier:.2f} ({multiplier_pct}%).'
		),
	), sim_data
