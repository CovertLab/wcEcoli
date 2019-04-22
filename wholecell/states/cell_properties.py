
from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.states.internal_state
import wholecell.views.view

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

		self._processIDs = sim.processes.keys()

		# initialize molecule IDs and concentrations
		self.property_ids = []
		self.property_values = np.array([])

		# create bulk container for molecule concentrations. This uses concentrations instead of counts.
		self.container = BulkObjectsContainer(self.property_ids, dtype=np.float64)
		self.container.countsIs(self.property_values)


	def update(self):
		'''update cell properties'''

		# TODO -- update cellMass, volume, countsToMolar


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
