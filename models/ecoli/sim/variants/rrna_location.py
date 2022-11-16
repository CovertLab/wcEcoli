"""
Variant to compare the impacts of changing the locations of rRNA genes on the
chromosome.

Modifies:
	sim_data.process.transcription.rna_data["replication_coordinate"]

Expected variant indices:
	0: control
	1: symmetrically flip the chromosomal positions of all rRNA-encoding geness
	such that most of them lie closer to terC than oriC (Note: there may be
	some genes that physically overlap as a result of this change)
"""

import numpy as np

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rrna_location(sim_data, index):
	if index == 1:
		is_rrna = sim_data.process.transcription.rna_data['is_rRNA']
		coordinate = sim_data.process.transcription.rna_data[
			"replication_coordinate"]
		replichore_lengths = sim_data.process.replication.replichore_lengths

		def flip_coordinates(orig_coordinates):
			if orig_coordinates >= 0:
				return replichore_lengths[0] - orig_coordinates
			else:
				return -replichore_lengths[1] - orig_coordinates

		# Flip coordinates of all rRNA genes
		for i in np.where(is_rrna)[0]:
			coordinate[i] = flip_coordinates(coordinate[i])

		return dict(
			shortName = "rrna_relocated",
			desc = "Simulation with relocated rRNA genes",
			), sim_data

	else:
		return CONTROL_OUTPUT, sim_data
