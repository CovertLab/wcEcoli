"""
Sweeps a scalar multiplier applied to the max kcat estimates for below-line
essential reactions.  Use this to test how the tightness of kcat-based flux
upper bounds affects growth relative to an unconstrained wildtype control.

Index 0 is a wildtype control (no kcat bounds applied).  Indices 1-15 use the
'max' quantile kcat estimates scaled by a multiplier from 1.0 down to 0.8
(15 evenly spaced values via np.linspace).

When kcat bounds are active, sim_data.process.metabolism is modified:
	sim_data.process.metabolism.kcat_estimate_quantile   (str)
	sim_data.process.metabolism.kcat_estimate_multiplier (float)
	sim_data.process.metabolism.selected_kcat_estimates  (dict)

Expected variant indices:
	 0: wildtype control (no kcat bounds)
	 1: max x 100%
	 2: max x  99%
	 3: max x  98%
	 4: max x  97%
	 5: max x  96%
	 6: max x  95%
	 7: max x  94%
	 8: max x  93%
	 9: max x  92%
	10: max x  91%
	11: max x  90%
	12: max x  89%
	13: max x  88%
	14: max x  87%
	15: max x  86%
"""

import numpy as np

KCAT_QUANTILE = 'max'
KCAT_MULTIPLIERS = np.round(np.arange(1.0, 0.855, -0.01), 2)


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
