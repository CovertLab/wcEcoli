#!/usr/bin/env python

"""
Environment.py

External State which represents the class of environmental molecules.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

from itertools import izip

import numpy as np

import wholecell.states.external_state
import wholecell.views.view

from wholecell.utils import units

from wholecell.utils.constants import REQUEST_PRIORITY_DEFAULT

ASSERT_POSITIVE_COUNTS = True

class NegativeCountsError(Exception):
	pass

class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):

		self._moleculeIDs = None

		super(Environment, self).__init__(*args, **kwargs)


	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		# Load constants
		self._moleculeIDs = sim_data.externalState.environment.environmentData['id']
		self._nutrientData = sim_data.externalState.environment.nutrientData
		self._condition = sim_data.externalState.environment.condition
		self._nutrientsTimeSeriesLabel = sim_data.externalState.environment.nutrientsTimeSeriesLabel


	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		tableWriter.writeAttributes(
			processNames = self._processIDs,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			counts = self.container._counts,
			)

	def tableLoad(self, tableReader, tableIndex):
		self.container.tableLoad(tableReader, tableIndex)


class EnvironmentViewBase(wholecell.views.view.View):
	_stateID = 'Environment'

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


class EnvironmentView(EnvironmentViewBase):
	def __init__(self, *args, **kwargs):
		super(EnvironmentView, self).__init__(*args, **kwargs)

		# State references
		assert len(set(self._query)) == len(self._query), "Environment molecules views cannot contain duplicate entries"
		self._containerIndexes = self._state.container._namesToIndexes(self._query)


	def _dataSize(self):
		return len(self._query)


	def counts(self):
		return self._counts()


	def mass(self):
		return self._mass()


	def countsIs(self, values):
		self._countsIs(values)


	def countsInc(self, values):
		self._countsInc(values)


	def countsDec(self, values):
		self._countsDec(values)


class EnvironmentView(EnvironmentViewBase):
	def __init__(self, *args, **kwargs):
		super(EnvironmentView, self).__init__(*args, **kwargs)

		# State references
		self._containerIndexes = self._state.container._namesToIndexes((self._query,))


	def count(self):
		return self._counts()[0]


	def mass(self):
		return self._mass()


	def countIs(self, value):
		self._countsIs(value)


	def countInc(self, value):
		self._countsInc(value)


	def countDec(self, value):
		self._countsDec(value)

