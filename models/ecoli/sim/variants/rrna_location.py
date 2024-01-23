"""
Variant to compare the impacts of changing the locations of rRNA genes on the
chromosome in multiple media conditions.

Modifies:
	sim_data.process.transcription.rna_data["replication_coordinate"]

Expected variant indices:
	0: control
	1: symmetrically flip the chromosomal positions of all rRNA-encoding genes
	such that most of them lie closer to terC than oriC (Note: there may be
	some genes that physically overlap as a result of this change)
	2: same as 1, in rich media conditions
	3: same as 1, in minimal-to-rich media shift conditions
"""

import numpy as np

from wholecell.utils import units

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rrna_location(sim_data, index):
	if index > 0:
		rna_is_rrna = sim_data.process.transcription.rna_data['is_rRNA']
		rna_coordinates = sim_data.process.transcription.rna_data[
			"replication_coordinate"]
		rna_lengths = sim_data.process.transcription.rna_data['length'].asNumber(
			units.nt)
		gene_is_rrna = sim_data.process.transcription.cistron_data['is_rRNA']
		gene_coordinates = sim_data.process.transcription.cistron_data[
			"replication_coordinate"]
		gene_lengths = sim_data.process.transcription.cistron_data['length'].asNumber(
			units.nt)
		replichore_lengths = sim_data.process.replication.replichore_lengths

		def flip_coordinates(orig_coordinates, length):
			# Note: this assumes the transcription direction of all rRNA genes
			# and transcripts are facing away from the origin
			if orig_coordinates >= 0:
				return replichore_lengths[0] - orig_coordinates - length
			else:
				return -replichore_lengths[1] - orig_coordinates + length

		# Flip coordinates of all rRNA transcription units and genes
		for i in np.where(rna_is_rrna)[0]:
			rna_coordinates[i] = flip_coordinates(
				rna_coordinates[i], rna_lengths[i])

		for i in np.where(gene_is_rrna)[0]:
			gene_coordinates[i] = flip_coordinates(
				gene_coordinates[i], gene_lengths[i])

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
			shortName = f"rrna_relocated_{condition_id}",
			desc = f"Simulation with relocated rRNA genes in {condition_id}",
			), sim_data

	else:
		return CONTROL_OUTPUT, sim_data
