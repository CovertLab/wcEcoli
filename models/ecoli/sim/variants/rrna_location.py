"""
Variant to compare the impacts of changing the locations of rRNA genes on the
chromosome.

Modifies:
	sim_data.process.transcription.rna_data["replication_coordinate"]

Expected variant indices:
	0: control
	1: symmetrically flip the chromosomal positions of all rRNA-encoding genes
	such that most of them lie closer to terC than oriC (Note: there may be
	some genes that physically overlap as a result of this change)
"""

import numpy as np

from wholecell.utils import units

CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)


def rrna_location(sim_data, index):
	if index == 1:
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

		return dict(
			shortName = "rrna_relocated",
			desc = "Simulation with relocated rRNA genes",
			), sim_data

	else:
		return CONTROL_OUTPUT, sim_data
