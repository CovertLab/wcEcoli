"""
RnaMaturation process

- Converts unprocessed tRNA/rRNA molecules into mature tRNA/rRNAs
"""

import numpy as np

import wholecell.processes.process


class RnaMaturation(wholecell.processes.process.Process):
	""" RnaMaturation """

	_name = "RnaMaturation"

	# Constructor
	def __init__(self):
		super(RnaMaturation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaMaturation, self).initialize(sim, sim_data)

		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data

		# Get matrices and vectors that describe reaction stoichiometries
		self.stoich_matrix = transcription.rna_maturation_stoich_matrix
		self.degraded_nt_counts = transcription.rna_maturation_degraded_nt_counts
		self.n_ppi_added = self.stoich_matrix.toarray().sum(axis=0) - 1

		# Build views
		unprocessed_rna_ids = rna_data['id'][rna_data['is_unprocessed']]
		mature_rna_ids = transcription.mature_rna_data['id']

		self.unprocessed_rnas = self.bulkMoleculesView(unprocessed_rna_ids)
		self.mature_rnas = self.bulkMoleculesView(mature_rna_ids)
		self.fragment_bases = self.bulkMoleculesView(sim_data.molecule_groups.polymerized_ntps)
		self.ppi = self.bulkMoleculeView(sim_data.molecule_ids.ppi)


	def calculateRequest(self):
		unprocessed_rna_counts = self.unprocessed_rnas.total_counts()

		# Request all unprocessed RNAs
		self.unprocessed_rnas.requestAll()

		# Request ppis that need to be added to the 5'-ends of mature RNAs
		self.ppi.requestIs(self.n_ppi_added.dot(unprocessed_rna_counts))


	def evolveState(self):
		# Get counts of unprocessed RNAs
		unprocessed_rna_counts = self.unprocessed_rnas.counts()

		# Calculate numbers of mature RNAs and fragment bases that are generated
		# upon maturation
		n_mature_rnas = self.stoich_matrix.dot(unprocessed_rna_counts)
		n_fragment_bases = np.dot(
			self.degraded_nt_counts.T, unprocessed_rna_counts)

		# Evolve states
		self.mature_rnas.countsInc(n_mature_rnas)
		self.fragment_bases.countsInc(n_fragment_bases)
		self.unprocessed_rnas.countsDec(unprocessed_rna_counts)
		self.ppi.countDec(self.n_ppi_added.dot(unprocessed_rna_counts))

		# Write to listener
		self.writeToListener(
			"RnaMaturationListener", "total_maturation_events",
			unprocessed_rna_counts.sum())
		self.writeToListener(
			"RnaMaturationListener", "total_degraded_ntps",
			n_fragment_bases.sum())
		self.writeToListener(
			"RnaMaturationListener", "unprocessed_rnas_consumed",
			unprocessed_rna_counts)
		self.writeToListener(
			"RnaMaturationListener", "mature_rnas_generated",
			n_mature_rnas)
