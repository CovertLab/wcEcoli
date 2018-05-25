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


class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):

		self._moleculeIDs = None

		super(Environment, self).__init__(*args, **kwargs)


	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self.nutrientsTimeSeries = sim_data.nutrientsTimeSeries
		self.nutrientData = sim_data.externalState.environment.nutrientData

		# current condition and timeseries
		self.nutrientsTimeSeriesLabel = sim_data.externalState.environment.nutrientsTimeSeriesLabel
		self.condition = self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel][0][1]


	def update(self):

		while len(self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel]) and \
			self.time() > self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel][0][0]:
				_, self.condition = self.nutrientsTimeSeries[self.nutrientsTimeSeriesLabel].popleft()


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

