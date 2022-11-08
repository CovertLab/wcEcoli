"""
RnaMaturation listener
"""

import numpy as np

import wholecell.listeners.listener

class RnaMaturationListener(wholecell.listeners.listener.Listener):
	""" RnaMaturationListener """

	_name = 'RnaMaturationListener'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RnaMaturationListener, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaMaturationListener, self).initialize(sim, sim_data)

		transcription = sim_data.process.transcription
		self.unprocessed_rna_ids = transcription.rna_data['id'][
			transcription.rna_data['is_unprocessed']]
		self.mature_rna_ids = transcription.mature_rna_data['id']
		self.enzyme_ids = transcription.rna_maturation_enzymes

	# Allocate memory
	def allocate(self):
		super(RnaMaturationListener, self).allocate()

		self.total_maturation_events = 0
		self.total_degraded_ntps = 0
		self.unprocessed_rnas_consumed = np.zeros(
			len(self.unprocessed_rna_ids), dtype=np.float64)
		self.mature_rnas_generated = np.zeros(
			len(self.mature_rna_ids), dtype=np.float64)
		self.maturation_enzyme_counts = np.zeros(
			len(self.enzyme_ids), dtype=np.float64)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		subcolumns = {
			'unprocessed_rnas_consumed': 'unprocessed_rna_ids',
			'mature_rnas_generated': 'mature_rna_ids',
			'maturation_enzyme_counts': 'enzyme_ids',
			}

		tableWriter.writeAttributes(
			unprocessed_rna_ids=list(self.unprocessed_rna_ids),
			mature_rna_ids=list(self.mature_rna_ids),
			enzyme_ids=self.enzyme_ids,
			subcolumns=subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			total_maturation_events = self.total_maturation_events,
			total_degraded_ntps = self.total_degraded_ntps,
			unprocessed_rnas_consumed = self.unprocessed_rnas_consumed,
			mature_rnas_generated = self.mature_rnas_generated,
			maturation_enzyme_counts = self.maturation_enzyme_counts,
			)
