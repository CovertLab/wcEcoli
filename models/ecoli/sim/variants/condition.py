"""
Condition variant for simulations in different environmental conditions

Modifies:
	sim_data.condition
	sim_data.external_state.current_timeline_id

Expected variant indices (dependent on sorted order of sim_data.conditionActiveTfs):
	0: acetate
	1: control
	2: anaerobic
	3: with amino acids
"""
from __future__ import absolute_import, division, print_function


def condition(sim_data, index):
	condition_labels = sorted(sim_data.conditionActiveTfs)
	condition_label = condition_labels[index]
	sim_data.condition = condition_label
	sim_data.external_state.current_timeline_id = condition_label
	sim_data.external_state.saved_timelines[condition_label] = [
		(0, sim_data.conditions[condition_label]["nutrients"])
	]

	return dict(
		shortName = "{}_env".format(condition_label),
		desc = "Simulation of condition {}.".format(condition_label)
		), sim_data
