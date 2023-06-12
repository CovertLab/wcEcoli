"""
RnaMaturation process

- Converts unprocessed tRNA/rRNA molecules into mature tRNA/rRNAs
- Consolidates the different variants of 23S, 16S, and 5S rRNAs into the single
variant that is used for ribosomal subunits
"""

from itertools import chain

import numpy as np

import wholecell.processes.process
from wholecell.utils import units


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
		mature_rna_data = transcription.mature_rna_data

		# Get matrices and vectors that describe maturation reactions
		self.stoich_matrix = transcription.rna_maturation_stoich_matrix
		self.enzyme_matrix = transcription.rna_maturation_enzyme_matrix.astype(int)
		self.n_required_enzymes = self.enzyme_matrix.sum(axis=1)
		self.degraded_nt_counts = transcription.rna_maturation_degraded_nt_counts
		self.n_ppi_added = self.stoich_matrix.toarray().sum(axis=0) - 1

		# Calculate number of NMPs that should be added when consolidating rRNA
		# molecules
		self.main_23s_rRNA_id = sim_data.molecule_groups.s50_23s_rRNA[0]
		self.main_16s_rRNA_id = sim_data.molecule_groups.s30_16s_rRNA[0]
		self.main_5s_rRNA_id = sim_data.molecule_groups.s50_5s_rRNA[0]

		self.variant_23s_rRNA_ids = sim_data.molecule_groups.s50_23s_rRNA[1:]
		self.variant_16s_rRNA_ids = sim_data.molecule_groups.s30_16s_rRNA[1:]
		self.variant_5s_rRNA_ids = sim_data.molecule_groups.s50_5s_rRNA[1:]

		counts_ACGU = np.vstack((
			rna_data['counts_ACGU'].asNumber(units.nt),
			mature_rna_data['counts_ACGU'].asNumber(units.nt)
			))
		rna_id_to_index = {
			rna_id: i for i, rna_id
			in enumerate(chain(rna_data['id'], mature_rna_data['id']))}
		def calculate_delta_nt_counts(main_id, variant_ids):
			main_index = rna_id_to_index[main_id]
			variant_indexes = np.array([
				rna_id_to_index[rna_id] for rna_id in variant_ids])

			delta_nt_counts = counts_ACGU[variant_indexes, :] - counts_ACGU[main_index, :]
			return delta_nt_counts

		self.delta_nt_counts_23s = calculate_delta_nt_counts(
			self.main_23s_rRNA_id, self.variant_23s_rRNA_ids)
		self.delta_nt_counts_16s = calculate_delta_nt_counts(
			self.main_16s_rRNA_id, self.variant_16s_rRNA_ids)
		self.delta_nt_counts_5s = calculate_delta_nt_counts(
			self.main_5s_rRNA_id, self.variant_5s_rRNA_ids)

		# Build views
		unprocessed_rna_ids = rna_data['id'][rna_data['is_unprocessed']]
		mature_rna_ids = transcription.mature_rna_data['id']
		rna_maturation_enzyme_ids = transcription.rna_maturation_enzymes

		self.unprocessed_rnas = self.bulkMoleculesView(unprocessed_rna_ids)
		self.mature_rnas = self.bulkMoleculesView(mature_rna_ids)
		self.rna_maturation_enzymes = self.bulkMoleculesView(
			rna_maturation_enzyme_ids)
		self.fragment_bases = self.bulkMoleculesView(
			sim_data.molecule_groups.polymerized_ntps)
		self.ppi = self.bulkMoleculeView(sim_data.molecule_ids.ppi)
		self.water = self.bulkMoleculeView(sim_data.molecule_ids.water)
		self.nmps = self.bulkMoleculesView(sim_data.molecule_groups.nmps)
		self.proton = self.bulkMoleculeView(sim_data.molecule_ids.proton)

		self.main_23s_rRNA = self.bulkMoleculeView(self.main_23s_rRNA_id)
		self.main_16s_rRNA = self.bulkMoleculeView(self.main_16s_rRNA_id)
		self.main_5s_rRNA = self.bulkMoleculeView(self.main_5s_rRNA_id)

		self.variant_23s_rRNAs = self.bulkMoleculesView(
			self.variant_23s_rRNA_ids)
		self.variant_16s_rRNAs = self.bulkMoleculesView(
			self.variant_16s_rRNA_ids)
		self.variant_5s_rRNAs = self.bulkMoleculesView(self.variant_5s_rRNA_ids)


	def calculateRequest(self):
		unprocessed_rna_counts = self.unprocessed_rnas.total_counts()
		variant_23s_rRNA_counts = self.variant_23s_rRNAs.total_counts()
		variant_16s_rRNA_counts = self.variant_16s_rRNAs.total_counts()
		variant_5s_rRNA_counts = self.variant_5s_rRNAs.total_counts()
		self.enzyme_availability = self.rna_maturation_enzymes.total_counts().astype(bool)

		# Determine which maturation reactions to turn off based on enzyme
		# availability
		reaction_is_off = (
			self.enzyme_matrix.dot(self.enzyme_availability)
			< self.n_required_enzymes
		)
		unprocessed_rna_counts[reaction_is_off] = 0

		# Request all unprocessed RNAs
		self.unprocessed_rnas.requestIs(unprocessed_rna_counts)

		# Request ppis that need to be added to the 5'-ends of mature RNAs
		self.ppi.requestIs(self.n_ppi_added.dot(unprocessed_rna_counts))

		# Request all variant rRNAs
		self.variant_23s_rRNAs.requestAll()
		self.variant_16s_rRNAs.requestAll()
		self.variant_5s_rRNAs.requestAll()

		# Request NMPs, water, and proton needed to balance mass
		n_added_bases_from_maturation = np.dot(
			self.degraded_nt_counts.T, unprocessed_rna_counts)
		n_added_bases_from_consolidation = (
			self.delta_nt_counts_23s.T.dot(variant_23s_rRNA_counts)
			+ self.delta_nt_counts_16s.T.dot(variant_16s_rRNA_counts)
			+ self.delta_nt_counts_5s.T.dot(variant_5s_rRNA_counts)
		)
		n_added_bases = n_added_bases_from_maturation + n_added_bases_from_consolidation
		n_total_added_bases = n_added_bases.sum()

		self.nmps.requestIs(np.abs(-n_added_bases))
		if n_total_added_bases > 0:
			self.water.requestIs(n_total_added_bases)
		else:
			self.proton.requestIs(-n_total_added_bases)

	def evolveState(self):
		# Get counts of unprocessed RNAs
		unprocessed_rna_counts = self.unprocessed_rnas.counts()

		# Calculate numbers of mature RNAs and fragment bases that are generated
		# upon maturation
		n_mature_rnas = self.stoich_matrix.dot(unprocessed_rna_counts)
		n_added_bases_from_maturation = np.dot(
			self.degraded_nt_counts.T, unprocessed_rna_counts)

		# Evolve states
		self.mature_rnas.countsInc(n_mature_rnas)
		self.unprocessed_rnas.countsDec(unprocessed_rna_counts)
		self.ppi.countDec(self.n_ppi_added.dot(unprocessed_rna_counts))

		# Write to listener
		self.writeToListener(
			"RnaMaturationListener", "total_maturation_events",
			unprocessed_rna_counts.sum())
		self.writeToListener(
			"RnaMaturationListener", "total_degraded_ntps",
			n_added_bases_from_maturation.sum())
		self.writeToListener(
			"RnaMaturationListener", "unprocessed_rnas_consumed",
			unprocessed_rna_counts)
		self.writeToListener(
			"RnaMaturationListener", "mature_rnas_generated",
			n_mature_rnas)
		self.writeToListener(
			"RnaMaturationListener", "maturation_enzyme_counts",
			self.rna_maturation_enzymes.total_counts())

		# Get counts of variant rRNAs
		variant_23s_rRNA_counts = self.variant_23s_rRNAs.counts()
		variant_16s_rRNA_counts = self.variant_16s_rRNAs.counts()
		variant_5s_rRNA_counts = self.variant_5s_rRNAs.counts()

		# Calculate number of NMPs that should be added to balance out the mass
		# difference during the consolidation
		n_added_bases_from_consolidation = (
			self.delta_nt_counts_23s.T.dot(variant_23s_rRNA_counts)
			+ self.delta_nt_counts_16s.T.dot(variant_16s_rRNA_counts)
			+ self.delta_nt_counts_5s.T.dot(variant_5s_rRNA_counts)
			)

		# Evolve states
		self.main_23s_rRNA.countInc(variant_23s_rRNA_counts.sum())
		self.main_16s_rRNA.countInc(variant_16s_rRNA_counts.sum())
		self.main_5s_rRNA.countInc(variant_5s_rRNA_counts.sum())
		self.variant_23s_rRNAs.countsDec(variant_23s_rRNA_counts)
		self.variant_16s_rRNAs.countsDec(variant_16s_rRNA_counts)
		self.variant_5s_rRNAs.countsDec(variant_5s_rRNA_counts)

		# Consume or add NMPs to balance out mass
		n_added_bases = n_added_bases_from_maturation + n_added_bases_from_consolidation
		n_total_added_bases = n_added_bases.sum()

		self.nmps.countsInc(n_added_bases)
		if n_total_added_bases > 0:
			self.water.countDec(n_total_added_bases)
			self.proton.countInc(n_total_added_bases)
		else:
			self.water.countInc(-n_total_added_bases)
			self.proton.countDec(-n_total_added_bases)
