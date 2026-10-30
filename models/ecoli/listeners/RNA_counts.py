#!/usr/bin/env python

"""
RNACounts Listener
"""

import numpy as np
import wholecell.listeners.listener


class RNACounts(wholecell.listeners.listener.Listener):
	"""
	Listener for the counts of each mRNA and rRNA transcription units and
	cistrons. Includes the counts of both partial and full transcripts.

	mRNA counts include only mRNAs that have not been targeted for
	degradation. mRNAs that are being degraded in the 5' to 3' direction (see
	DIRECTIONAL_MRNA_DECAY in rna_degradation.py) are counted separately in
	decaying_mRNA_counts. Cistron counts are calculated from transcription
	unit counts, so all cistrons of a partial transcript are counted.
	available_mRNA_cistron_counts instead counts only the cistrons that can
	bind new ribosomes: cistrons that have been transcribed and have not been
	passed by the degradation front.
	"""
	_name = 'RNACounts'

	def __init__(self, *args, **kwargs):
		super(RNACounts, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(RNACounts, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

		# Get IDs and indexes of all mRNA and rRNA transcription units
		self.all_TU_ids = sim_data.process.transcription.rna_data['id']
		self.mRNA_indexes = np.where(
			sim_data.process.transcription.rna_data['is_mRNA'])[0]
		self.mRNA_TU_ids = self.all_TU_ids[self.mRNA_indexes]
		self.rRNA_indexes = np.where(
			sim_data.process.transcription.rna_data['is_rRNA'])[0]
		self.rRNA_TU_ids = self.all_TU_ids[self.rRNA_indexes]

		# Get IDs and indexes of all mRNA and rRNA cistrons
		self.all_cistron_ids = sim_data.process.transcription.cistron_data['id']
		self.cistron_is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		self.mRNA_cistron_ids = self.all_cistron_ids[self.cistron_is_mRNA]
		self.cistron_is_rRNA = sim_data.process.transcription.cistron_data['is_rRNA']
		self.rRNA_cistron_ids = self.all_cistron_ids[self.cistron_is_rRNA]

		# Get mapping matrix between TUs and cistrons
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

		# Get start positions of each cistron within each TU, sorted by TU,
		# to count the cistrons available on partial and decaying mRNAs
		cistron_TU_pairs = sorted(
			(TU_index, cistron_index, start_position)
			for (cistron_index, TU_index), (start_position, _)
			in sim_data.process.transcription.cistron_start_end_pos_in_tu.items())
		pair_TU_indexes = np.array([pair[0] for pair in cistron_TU_pairs])
		self.pair_cistron_indexes = np.array([pair[1] for pair in cistron_TU_pairs])
		self.pair_start_positions = np.array([pair[2] for pair in cistron_TU_pairs])
		self.n_cistrons_per_TU = np.bincount(
			pair_TU_indexes, minlength=len(self.all_TU_ids))
		self.TU_pair_offsets = np.concatenate((
			[0], np.cumsum(self.n_cistrons_per_TU)[:-1]))

	def allocate(self):
		super(RNACounts, self).allocate()

		self.mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.full_mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.partial_mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.decaying_mRNA_counts = np.zeros(len(self.mRNA_TU_ids), dtype=np.int64)
		self.mRNA_cistron_counts = np.zeros(len(self.mRNA_cistron_ids), dtype=np.int64)
		self.full_mRNA_cistron_counts = np.zeros(len(self.mRNA_cistron_ids), dtype=np.int64)
		self.partial_mRNA_cistron_counts = np.zeros(len(self.mRNA_cistron_ids), dtype=np.int64)
		self.available_mRNA_cistron_counts = np.zeros(len(self.mRNA_cistron_ids), dtype=np.int64)
		self.partial_rRNA_counts = np.zeros(len(self.rRNA_TU_ids), dtype=np.int64)
		self.partial_rRNA_cistron_counts = np.zeros(len(self.rRNA_cistron_ids), dtype=np.int64)

	def update(self):
		# Get attributes of RNAs
		RNAs = self.uniqueMolecules.container.objectsInCollection('RNA')
		TU_indexes, can_translate, is_full_transcript, is_decaying, transcript_lengths, decay_front = RNAs.attrs(
			'TU_index', 'can_translate', 'is_full_transcript', 'is_decaying',
			'transcript_length', 'decay_front')
		# mRNAs being degraded in the 5' to 3' direction are counted separately
		is_intact_mRNA = can_translate & ~is_decaying
		is_rRNA = np.isin(TU_indexes, self.rRNA_indexes)

		# Get counts of mRNA and rRNA transcription units
		all_TU_counts = np.bincount(
			TU_indexes[np.logical_or(is_intact_mRNA, is_rRNA)], minlength=len(self.all_TU_ids))
		self.mRNA_counts = all_TU_counts[self.mRNA_indexes]
		full_TU_counts = np.bincount(
			TU_indexes[np.logical_and(is_intact_mRNA, is_full_transcript)],
			minlength=len(self.all_TU_ids))
		self.full_mRNA_counts = full_TU_counts[self.mRNA_indexes]
		partial_TU_counts = all_TU_counts - full_TU_counts
		self.partial_mRNA_counts = self.mRNA_counts - self.full_mRNA_counts
		self.decaying_mRNA_counts = np.bincount(
			TU_indexes[can_translate & is_decaying],
			minlength=len(self.all_TU_ids))[self.mRNA_indexes]
		# All unique rRNAs are partially transcribed
		self.partial_rRNA_counts = all_TU_counts[self.rRNA_indexes]

		# Calculate counts of mRNA cistrons from transcription unit counts
		# TODO (ggsun): Partial RNA cistron counts should take into account
		# 	the lengths of each RNA transcript.
		cistron_counts = self.cistron_tu_mapping_matrix.dot(
			all_TU_counts)
		self.mRNA_cistron_counts = cistron_counts[self.cistron_is_mRNA]
		self.full_mRNA_cistron_counts = self.cistron_tu_mapping_matrix.dot(
			full_TU_counts)[self.cistron_is_mRNA]
		self.partial_mRNA_cistron_counts = self.cistron_tu_mapping_matrix.dot(
			partial_TU_counts)[self.cistron_is_mRNA]
		self.partial_rRNA_cistron_counts = cistron_counts[self.cistron_is_rRNA]

		# Calculate counts of mRNA cistrons that can bind new ribosomes. All
		# cistrons on full transcripts that are not being degraded are
		# available. Cistrons on partial transcripts and on transcripts being
		# degraded are available if they have been transcribed and have not
		# been passed by the degradation front.
		is_other_mRNA = can_translate & ~(is_full_transcript & ~is_decaying)
		other_decay_front = np.where(
			is_decaying[is_other_mRNA], decay_front[is_other_mRNA], -1)
		self.available_mRNA_cistron_counts = self.full_mRNA_cistron_counts + \
			self._count_available_cistrons(
				TU_indexes[is_other_mRNA], transcript_lengths[is_other_mRNA],
				other_decay_front)[self.cistron_is_mRNA]

	def _count_available_cistrons(self, TU_indexes, transcript_lengths, decay_front):
		"""
		Counts the cistrons on a set of mRNAs that have been transcribed and
		have not been passed by the degradation front.

		Args:
			TU_indexes: TU index of each mRNA
			transcript_lengths: Current length of each mRNA
			decay_front: Position of the degradation front of each mRNA (-1
				for mRNAs that are not being degraded)

		Returns:
			Counts of each cistron
		"""
		n_cistrons = self.n_cistrons_per_TU[TU_indexes]
		mRNA_rep = np.repeat(np.arange(len(TU_indexes)), n_cistrons)
		pair_indexes = (np.repeat(self.TU_pair_offsets[TU_indexes], n_cistrons)
			+ np.arange(n_cistrons.sum())
			- np.repeat(np.cumsum(n_cistrons) - n_cistrons, n_cistrons))
		start_positions = self.pair_start_positions[pair_indexes]
		is_available = (
			(transcript_lengths[mRNA_rep] > start_positions)
			& (start_positions >= decay_front[mRNA_rep]))
		return np.bincount(
			self.pair_cistron_indexes[pair_indexes][is_available],
			minlength=len(self.all_cistron_ids))

	def tableCreate(self, tableWriter):
		subcolumns = {
			'mRNA_counts': 'mRNA_ids',
			'full_mRNA_counts': 'mRNA_ids',
			'partial_mRNA_counts': 'mRNA_ids',
			'decaying_mRNA_counts': 'mRNA_ids',
			'mRNA_cistron_counts': 'mRNA_cistron_ids',
			'full_mRNA_cistron_counts': 'mRNA_cistron_ids',
			'partial_mRNA_cistron_counts': 'mRNA_cistron_ids',
			'available_mRNA_cistron_counts': 'mRNA_cistron_ids',
			'partial_rRNA_counts': 'rRNA_ids',
			'partial_rRNA_cistron_counts': 'rRNA_cistron_ids',
			}

		tableWriter.writeAttributes(
			mRNA_ids = self.mRNA_TU_ids.tolist(),
			mRNA_cistron_ids = self.mRNA_cistron_ids.tolist(),
			rRNA_ids = self.rRNA_TU_ids.tolist(),
			rRNA_cistron_ids = self.rRNA_cistron_ids.tolist(),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			mRNA_counts = self.mRNA_counts,
			full_mRNA_counts = self.full_mRNA_counts,
			partial_mRNA_counts = self.partial_mRNA_counts,
			decaying_mRNA_counts = self.decaying_mRNA_counts,
			mRNA_cistron_counts = self.mRNA_cistron_counts,
			full_mRNA_cistron_counts = self.full_mRNA_cistron_counts,
			partial_mRNA_cistron_counts = self.partial_mRNA_cistron_counts,
			available_mRNA_cistron_counts = self.available_mRNA_cistron_counts,
			partial_rRNA_counts = self.partial_rRNA_counts,
			partial_rRNA_cistron_counts = self.partial_rRNA_cistron_counts,
			)
