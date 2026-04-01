"""
Compares smoothed_max_buffered (SMB) vs smoothed_max without buffered
reactions (SM w/o B) at several multiplier levels.

SMB uses all kcat estimates from smoothed_max_buffered, which includes a
1.1x headroom multiplier on low-variance reactions.

SM w/o B uses smoothed_max estimates but *excludes* reactions that were
marked as buffered, leaving them unconstrained.  This tests whether the
buffered reactions are better left unconstrained entirely vs given a
small headroom multiplier.

Expected variant indices:
	0: wildtype control (no kcat bounds)
	1: 1.00 x SMB
	2: 1.00 x SM w/o B
	3: 0.75 x SMB
	4: 0.75 x SM w/o B
	5: 0.50 x SMB
	6: 0.50 x SM w/o B
	7: 0.25 x SMB
	8: 0.25 x SM w/o B
"""

# Multipliers paired as (SMB, SM w/o B) for each tightness level.
MULTIPLIERS = [1.0, 0.75, 0.50, 0.25]


def kcat_estimate_scale(sim_data, index):
	if index == 0:
		return dict(
			shortName='wildtype_control',
			desc='No kcat bounds applied; matches wildtype behavior.',
		), sim_data

	# Odd indices (1, 3, 5, 7) = SMB; even indices (2, 4, 6, 8) = SM w/o B.
	pair_idx = (index - 1) // 2
	use_buffered = (index % 2 == 1)
	multiplier = MULTIPLIERS[pair_idx]
	multiplier_pct = int(round(multiplier * 100))

	metabolism = sim_data.process.metabolism
	buffered_keys = metabolism.kcat_buffered_reactions

	if use_buffered:
		quantile = 'smoothed_max_buffered'
		base_estimates = metabolism.kcat_estimates[quantile]
		label = 'SMB'
	else:
		quantile = 'smoothed_max'
		base_estimates = {
			key: value
			for key, value in metabolism.kcat_estimates[quantile].items()
			if key not in buffered_keys
		}
		label = 'SM_wo_B'

	metabolism.kcat_estimate_quantile = quantile
	metabolism.kcat_estimate_multiplier = multiplier
	metabolism.selected_kcat_estimates = {
		key: value * multiplier for key, value in base_estimates.items()
	}

	return dict(
		shortName=f'{label}_x{multiplier_pct}pct',
		desc=(
			f'kcat estimates from {quantile},'
			f' {"excluding" if not use_buffered else "including"}'
			f' buffered reactions,'
			f' scaled by {multiplier:.2f} ({multiplier_pct}%).'
		),
	), sim_data
