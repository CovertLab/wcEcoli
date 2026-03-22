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

		# Get the IDs of molecules involved in two component system reactions:
		self.two_component_system_molecule_ids = (
			list(sim_data.process.two_component_system.molecule_names))


	def allocate(self):
		super(TwoComponentSystems, self).allocate()

		self.delta2CSMolecules = np.zeros(
			len(self.two_component_system_molecule_ids),
			np.int64)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'delta2CSMolecules': 'twoComponentSystemMoleculeIds'
			}

		tableWriter.writeAttributes(
			twoComponentSystemMoleculeIds = self.two_component_system_molecule_ids,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			delta2CSMolecules = self.delta2CSMolecules,
			)
