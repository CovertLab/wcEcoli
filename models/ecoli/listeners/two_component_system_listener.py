"""
Two Component Systems Listener
"""

import numpy as np
import wholecell.listeners.listener

class TwoComponentSystems(wholecell.listeners.listener.Listener):
	"""
	Listener for the reaction events for two component system reactions.
	"""
	_name = 'TwoComponentSystems'

	def __init__(self, *args, **kwargs):
		super(TwoComponentSystems, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(TwoComponentSystems, self).initialize(sim, sim_data)

		# Get the IDs of two component systems reactions:
		self.reactionIDs = list(sim_data.process.two_component_system.rxn_ids)


	def allocate(self):
		super(TwoComponentSystems, self).allocate()

		self.complexationEvents = np.zeros(
			len(self.reactionIDs),
			np.int64)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'complexationEvents': 'reactionIDs'}

		tableWriter.writeAttributes(
			twoComponentSystemRxnIds = self.reactionIDs,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			complexationEvents = self.complexationEvents,
			)
