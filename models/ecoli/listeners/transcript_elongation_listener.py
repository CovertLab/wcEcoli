"""
TranscriptElongationListener
"""

import numpy as np

import wholecell.listeners.listener

class TranscriptElongationListener(wholecell.listeners.listener.Listener):
	""" TranscriptElongationListener """

	_name = 'TranscriptElongationListener'

	def __init__(self, *args, **kwargs):
		super(TranscriptElongationListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TranscriptElongationListener, self).initialize(sim, sim_data)

		# Attributes
		self.rnaIds = sim_data.process.transcription.rna_data['id']
		self.attenuated_rnas = self.rnaIds[sim_data.process.transcription.attenuated_rna_indices]
		n_attenuated = len(self.attenuated_rnas)
		self.cistron_ids = sim_data.process.transcription.cistron_data['id']
		self.n_cistrons = self.cistron_ids.size
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

		# Columns
		self.countRnaSynthesized = np.zeros(sim_data.process.transcription.rna_data.fullArray().size, np.int64)
		self.countRnaCistronSynthesized = np.zeros(self.n_cistrons, np.int64)
		self.countNTPsUSed = 0
		self.attenuation_probability = np.zeros(n_attenuated)
		self.counts_attenuated = np.zeros(n_attenuated, np.int64)

	def update(self):
		# Map completed transcripts from transcription units to cistrons. For
		# polycistronic operons, each completed TU is counted once for every
		# constituent cistron (mirrors RnapData.rna_init_event_per_cistron).
		# Unlike initiation events, these counts exclude transcripts lost to
		# tRNA attenuation, since attenuated transcripts never complete.
		self.countRnaCistronSynthesized = self.cistron_tu_mapping_matrix.dot(
			self.countRnaSynthesized)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'countRnaSynthesized': 'rnaIds',
			'countRnaCistronSynthesized': 'cistron_ids',
			'attenuation_probability': 'attenuated_rnas',
			'counts_attenuated': 'attenuated_rnas',
		}

		tableWriter.writeAttributes( # TODO: reconsider attribute names
			countRnaSynthesized = self.countUnits,
			countRnaCistronSynthesized = self.countUnits,
			countNTPsUSed = self.countUnits,
			rnaIds = list(self.rnaIds),
			cistron_ids = list(self.cistron_ids),
			attenuated_rnas = list(self.attenuated_rnas),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			countRnaSynthesized = self.countRnaSynthesized,
			countRnaCistronSynthesized = self.countRnaCistronSynthesized,
			countNTPsUSed = self.countNTPsUSed,
			attenuation_probability = self.attenuation_probability,
			counts_attenuated = self.counts_attenuated,
			)
