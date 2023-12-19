"""
Variant to compare the impacts of reversing the orientation of rRNA genes in
multiple media conditions.

Modifies:
	sim_data.process.transcription.rna_data["is_forward"]
	sim_data.process.transcription.rna_data["replication_coordinate"]
	sim_data.process.transcription.cistron_data["is_forward"]
	sim_data.process.transcription.cistron_data["replication_coordinate"]

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
		transcription = sim_data.process.transcription
		is_rrna = transcription.rna_data['is_rRNA']
		rna_ids = transcription.rna_data['id']
		rna_is_forward = transcription.rna_data['is_forward']
		rna_coordinates = transcription.rna_data[
			"replication_coordinate"]
		rna_length = transcription.rna_data["length"].asNumber(units.nt)
		cistron_is_forward = transcription.cistron_data['is_forward']
		cistron_coordinates = transcription.cistron_data[
			'replication_coordinate']

		# Reverse the ordering and orientations of genes within all rRNA operons
		for i in np.where(is_rrna)[0]:
			cistron_indexes = transcription.rna_id_to_cistron_indexes(rna_ids[i])
			cistron_pos_within_operon = np.abs(
				cistron_coordinates[cistron_indexes] - rna_coordinates[i])

			cistron_coordinates[cistron_indexes] = (
				rna_coordinates[i]
				+ 2*(rna_is_forward[i] - 0.5)*(rna_length[i] - cistron_pos_within_operon)
				)
			cistron_is_forward[cistron_indexes] = ~cistron_is_forward[cistron_indexes]

		# Reverse orientations of all rRNA transcription units
		rna_coordinates[is_rrna] = rna_coordinates[is_rrna] + np.multiply(
			rna_length[is_rrna], 2*(rna_is_forward[is_rrna] - 0.5))
		rna_is_forward[is_rrna] = ~rna_is_forward[is_rrna]


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
