#!/usr/bin/env python

"""
environment.py

The Environment State handles the identity and dynamic properties of environmental
molecules in the simulation.

@author: Eran Agmon
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/14/2018
"""

from __future__ import division

import re
from itertools import izip

import numpy as np

import wholecell.internal_states.state
import wholecell.views.view

#TODO (Eran) create environment objects container
from wholecell.containers.environment_objects_container import EnvironmentObjectsContainer

from wholecell.utils import units

from wholecell.utils.constants import REQUEST_PRIORITY_DEFAULT

class Environment(wholecell.internal_states.state.State):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):
		self.container = None
		self._moleculeMass = None
		self._moleculeIDs = None

		super(Environment, self).__init__(*args, **kwargs)


	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self._processIDs = sim.processes.keys()

		# Load constants
		# TODO (Eran) create environment objects for sim_data
		self._moleculeIDs = sim_data.state.environment.environmentData['id']
		self._moleculeMass = sim_data.state.environment.environmentData['mass'].asNumber(units.fg / units.mol) / sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		self._submassNameToIndex = sim_data.submassNameToIndex

		# Create the container for molecule counts
		self.container = EnvironmenObjectsContainer(self._moleculeIDs)


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


	def _setCounts(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] = values


	def _increaseCounts(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

		self._state._countsAllocatedFinal[self._containerIndexes, self._processIndex] += values


	def _decreaseCounts(self, values):
		assert (np.size(values) == np.size(self._containerIndexes)) or np.size(values) == 1, 'Inappropriately sized values'

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


	def setCount(self, values):
		self._countsIs(values)


	def increaseCount(self, values):
		self._countsInc(values)


	def decreaseCount(self, values):
		self._countsDec(values)

