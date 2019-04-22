
from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.states.internal_state
import wholecell.views.view
from wholecell.utils import units

from wholecell.containers.bulk_objects_container import BulkObjectsContainer

ASSERT_POSITIVE_CONCENTRATIONS = True

class NegativeConcentrationError(Exception):
	pass


class CellProperties(wholecell.states.internal_state.InternalState):
	_name = 'CellProperties'

	def __init__(self, *args, **kwargs):
		self.container = None
		self.property_ids = None
		self.property_values = None

		super(CellProperties, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(CellProperties, self).initialize(sim, sim_data)

		# self._processIDs = sim.processes.keys()

		self.internal_states = sim.internal_states

		# load constants
		self.nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		self.cell_density = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		# initialize properties
		self.cell_mass = 0.0
		self.volume = self.cell_mass / self.cell_density
		self.counts_to_molar = 0.0  # (units.g * units.mol / units.L)
		# self.counts_to_molar = 1.0 / (self.nAvogadro * self.volume)

		# initialize properties IDs and values
		self.property_ids = [
			'cell_mass',
			'volume',
			'counts_to_molar',
			 ]
		self.property_values = np.array([
			self.cell_mass,
			self.volume,
			self.counts_to_molar,
		])

		# create bulk container for property values.
		self.container = BulkObjectsContainer(self.property_ids, dtype=np.float64)
		self.container.countsIs(self.property_values)


	def update(self):

		# update cell_mass, volume, countsToMolar
		masses = sum(state.mass() for state in self.internal_states.itervalues())
		postEvolveMasses = masses[1, ...]
		self.cell_mass = postEvolveMasses.sum() # sum over all dimensions
		self.volume = self.cell_mass / self.cell_density
		self.counts_to_molar = 1 / (self.nAvogadro * self.volume)

		self.property_values = np.array([
			self.cell_mass,
			self.volume,
			self.counts_to_molar,
		])

		# update property values in container
		self.container.countsIs(self.property_values)


	def calculatePreEvolveStateMass(self):
		pass


	def calculatePostEvolveStateMass(self):
		pass


	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		tableWriter.writeAttributes()


	def tableAppend(self, tableWriter):
		tableWriter.append()


class CellPropertiesViewBase(object):
	_stateID = 'CellProperties'

	def __init__(self, *args, **kwargs):
		super(CellPropertiesViewBase, self).__init__(*args, **kwargs)
		self._containerIndexes = None  # subclasses must set this

	def _updateQuery(self):
		self._totalIs(self._state.container._counts[self._containerIndexes])

	def _counts(self):
		return self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex].copy()

	def _countsIs(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] = values

	def _countsInc(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		values = np.asarray(values, dtype=self._containerIndexes.dtype)
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] += values

	def _countsDec(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		values = np.asarray(values, dtype=self._containerIndexes.dtype)
		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] -= values

	# Request
	def requestIs(self, value):
		self._requestedCount[:] = value

	def requestAll(self):
		self._requestedCount[:] = self._totalCount


class CellPropertiesView(CellPropertiesViewBase):
	def __init__(self, *args, **kwargs):
		super(CellPropertiesView, self).__init__(*args, **kwargs)

		# State references
		self._containerIndexes = self._state.container._namesToIndexes((self._query,))

	def count(self):
		return self._counts()[0]

	def countIs(self, value):
		self._countsIs(value)

	def countInc(self, value):
		self._countsInc(value)

	def countDec(self, value):
		self._countsDec(value)
