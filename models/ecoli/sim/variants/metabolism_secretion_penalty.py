"""
Variant to compare impact of adjusting the secretion penalty in metabolism.
Penalty varies between 0 and 10 (could go higher if needed).

Modifies:
	sim_data.process.metabolism.secretion_penalty_coeff

Expected variant indices (dependent on SECRETION_PENALTY):
	4: control
	0 (no penalty) - 10 (high penalty)
"""

SECRETION_PENALTY = [0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1, 10]


def metabolism_secretion_penalty(sim_data, index):
	penalty = SECRETION_PENALTY[index]
	sim_data.process.metabolism.secretion_penalty_coeff = penalty

	return dict(
		shortName="penalty={:.0E}".format(penalty),
		desc="Simulation with secretion penalty of {}.".format(penalty)
		), sim_data
