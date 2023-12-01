"""
Variant to compare the impacts of reversing the orientation of rRNA genes in
multiple media conditions.

Modifies:
	sim_data.process.transcription.rna_data["is_forward"]
	sim_data.process.transcription.rna_data["replication_coordinate"]

Expected variant indices:
	0: control
	1: reverse orientation of all rRNA-encoding genes, minimal media
	2: reverse orientation of all rRNA-encoding genes, rich media
	3: reverse orientation of all rRNA-encoding genes, minimal-to-rich media
		shift
"""

import numpy as np
from wholecell.utils import units

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rrna_orientation(sim_data, index):
	if index > 0:
		is_rrna = sim_data.process.transcription.rna_data['is_rRNA']
		is_forward = sim_data.process.transcription.rna_data['is_forward']
		coordinate = sim_data.process.transcription.rna_data[
			"replication_coordinate"]
		length = sim_data.process.transcription.rna_data["length"].asNumber(units.nt)

		# Reverse orientations of all rRNA genes
		coordinate[is_rrna] = coordinate[is_rrna] + np.multiply(
			length[is_rrna], is_forward[is_rrna])
		is_forward[is_rrna] = ~is_forward[is_rrna]

		# Change media conditions for indexes 2 and 3
		condition_id = sim_data.condition
		if index == 2:
			condition_labels = sim_data.ordered_conditions
			condition_id = condition_labels[1]  # Index for rich media condition
			sim_data.condition = condition_id
			sim_data.external_state.current_timeline_id = condition_id
			sim_data.external_state.saved_timelines[condition_id] = [
				(0, sim_data.conditions[condition_id]["nutrients"])
				]
		elif index == 3:
			saved_timelines = sim_data.external_state.saved_timelines
			timeline_ids = sorted(saved_timelines)
			condition_id = timeline_ids[28]  # Index for minimal-to-rich media shift
			sim_data.external_state.current_timeline_id = condition_id

			# Get possible condition from starting nutrients for proper
			# initialization
			nutrients = saved_timelines[condition_id][0][1]
			conditions = [cond for cond in sim_data.condition_active_tfs
				if sim_data.conditions[cond]['nutrients'] == nutrients]
			sim_data.condition = conditions[0]

		return dict(
			shortName = f"rrna_orientation_reversed_{condition_id}",
			desc = f"Simulation with the orientations of all rRNA genes reversed in {condition_id}"
			), sim_data

	else:
		return CONTROL_OUTPUT, sim_data
