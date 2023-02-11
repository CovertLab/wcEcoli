"""
RnaDegradationListener
"""

import numpy as np

import wholecell.listeners.listener

class RnaDegradationListener(wholecell.listeners.listener.Listener):
	""" RnaDegradationListener """

	_name = 'RnaDegradationListener'

	def __init__(self, *args, **kwargs):
		super(RnaDegradationListener, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaDegradationListener, self).initialize(sim, sim_data)
		self.rna_ids = sim_data.process.transcription.rna_data['id']
		self.mature_rna_ids = sim_data.process.transcription.mature_rna_data['id']
		self.n_transcribed_rnas = len(self.rna_ids)
		self.mature_rna_exists = (len(self.mature_rna_ids) > 0)
		self.cistron_ids = sim_data.process.transcription.cistron_data['id']
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

		cistron_id_to_index = {
			cistron_id: i for (i, cistron_id) in enumerate(self.cistron_ids)
			}
		self.mature_rna_cistron_indexes = np.array([
			cistron_id_to_index[rna_id[:-3]] for rna_id in self.mature_rna_ids
			])

	def allocate(self):
		super(RnaDegradationListener, self).allocate()

		self.count_RNA_degraded = np.zeros(len(self.rna_ids) + len(self.mature_rna_ids), np.int64)
		self.count_RNA_degraded_per_cistron = np.zeros(len(self.cistron_ids), np.int64)
		self.nucleotidesFromDegradation = 0
		self.FractionActiveEndoRNases = 0.
		self.DiffRelativeFirstOrderDecay = 0.
		self.FractEndoRRnaCounts = 0.
		self.fragmentBasesDigested = 0

	def update(self):
		count_RNA_degraded_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.count_RNA_degraded[:self.n_transcribed_rnas])

		# Add degraded counts from mature RNAs
		if self.mature_rna_exists:
			count_RNA_degraded_per_cistron[self.mature_rna_cistron_indexes] += self.count_RNA_degraded[self.n_transcribed_rnas:]

		self.count_RNA_degraded_per_cistron = count_RNA_degraded_per_cistron

	def tableCreate(self, tableWriter):
		subcolumns = {
			'count_RNA_degraded': 'rna_ids',
			'count_RNA_degraded_per_cistron': 'cistron_ids'}

		tableWriter.writeAttributes(
			rna_ids = list(self.rna_ids) + list(self.mature_rna_ids),
			cistron_ids = list(self.cistron_ids),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			count_RNA_degraded = self.count_RNA_degraded,
			count_RNA_degraded_per_cistron = self.count_RNA_degraded_per_cistron,
			nucleotidesFromDegradation = self.nucleotidesFromDegradation,
			FractionActiveEndoRNases = self.FractionActiveEndoRNases,
			DiffRelativeFirstOrderDecay = self.DiffRelativeFirstOrderDecay,
			FractEndoRRnaCounts = self.FractEndoRRnaCounts,
			fragmentBasesDigested = self.fragmentBasesDigested,
			)
