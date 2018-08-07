#!/usr/bin/env python

"""
External state that represents environmental molecules and conditions.

	- nutrients_time_series: a list of tuples that include time and nutrients in
		which shifts occur.
	- nutrients: a string specifying the current nutrient condition.
	- times: a list of all times at which the nutrients shift.

	Functions:
	----------
	- update: updates nutrients according to nutrients_time_series

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy as np

import wholecell.states.external_state
import wholecell.views.view

from wholecell.utils import units
from wholecell.containers.environment_objects_container import EnvironmentObjectsContainer

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L

ASSERT_POSITIVE_CONCENTRATIONS = True

class NegativeConcentrationError(Exception):
	pass

class Environment(wholecell.states.external_state.ExternalState):
	_name = 'Environment'

	def __init__(self, *args, **kwargs):
		self.container = None
		self._moleculeIDs = None
		self._concentrations = None
		self._volume = None
		self._infinite_environment = None

		super(Environment, self).__init__(*args, **kwargs)


	def initialize(self, sim, sim_data):
		super(Environment, self).initialize(sim, sim_data)

		self._processIDs = sim.processes.keys()

		# load constants
		self._nAvogadro = sim_data.constants.nAvogadro

		# get molecule IDs and initial concentrations
		moleculeIDs = [id for id, value in sim_data.external_state.environment.nutrients.iteritems()]
		concentrations = np.array([value.asNumber() for id, value in sim_data.external_state.environment.nutrients.iteritems()])

		self.setLocalEnvironment(moleculeIDs, concentrations)

		# environment time series data
		self.environment_dict = sim_data.external_state.environment.environment_dict
		self.nutrients_time_series_label = sim_data.external_state.environment.nutrients_time_series_label
		self.nutrients_time_series = sim_data.external_state.environment.nutrients_time_series[
			self.nutrients_time_series_label
			]

		self.nutrients = self.nutrients_time_series[0][1]
		self._times = [t[0] for t in self.nutrients_time_series]

		# get volume if volume is infinite (default), countsInc is skipped
		self._volume = self.nutrients_time_series[0][2]
		if np.isinf(self._volume):
			self._infinite_environment = True
		else:
			self._infinite_environment = False
			self._volume = self._volume * (units.L)

		# create container for molecule concentrations
		self.container = EnvironmentObjectsContainer(self._moleculeIDs)
		self.container.concentrationsIs(self._concentrations)
		self.container.volumeIs(self._volume)

		# the length of the longest nutrients name, for padding in nutrients listener
		self._nutrients_name_max_length = len(max([t[1] for t in self.nutrients_time_series], key=len))


	def setLocalEnvironment(self, moleculeIDs, concentrations):
		self._moleculeIDs = moleculeIDs
		self._concentrations = concentrations


	def update(self):
		current_index = [i for i, t in enumerate(self._times) if self.time()>=t][-1]

		if self.nutrients != self.nutrients_time_series[current_index][1]:
			self.nutrients = self.nutrients_time_series[current_index][1]
			self._concentrations = np.array([value.asNumber() for id, value in self.environment_dict[self.nutrients].iteritems()])
			self.container.concentrationsIs(self._concentrations)

			self._volume = self.nutrients_time_series[current_index][2]
			if np.isinf(self._volume):
				self._infinite_environment = True
			else:
				self._infinite_environment = False
				self._volume = float(self._volume) * (units.L)

		if ASSERT_POSITIVE_CONCENTRATIONS and not (self._concentrations >= 0).all():
			raise NegativeConcentrationError(
					"Negative environment concentration(s) in self._concentrations:\n"
					+ "\n".join(
					"{}".format(
						self._moleculeIDs[molIndex],
						)
					for molIndex in np.where(self._concentrations < 0)[0]
					)
				)


	def _counts_to_concentration(self, counts):
		concentrations = counts / (self._volume * self._nAvogadro).asNumber(VOLUME_UNITS / COUNTS_UNITS)
		return concentrations


	def tableCreate(self, tableWriter):
		self.container.tableCreate(tableWriter)
		tableWriter.writeAttributes(
			nutrientTimeSeriesLabel = self.nutrients_time_series_label,
			)

	def tableAppend(self, tableWriter):

		tableWriter.append(
			nutrientCondition = self.nutrients.ljust(self._nutrients_name_max_length),
			nutrientConcentrations=self._concentrations,
			volume=self._volume,
			)


class EnvironmentViewBase(object):
	_stateID = 'Environment'

	def __init__(self, state, process, query): # weight, priority, coupling id, option to not evaluate the query
		self._state = state
		self._state.viewAdd(self)
		self._processId = process.name()
		self._processIndex = process._processIndex
		self._query = query
		self._concentrations = np.zeros(self._dataSize(), np.float64) # number of objects that satisfy the query


	# Interface to State
	def _updateQuery(self):
		self._totalIs(self._state.container._concentrations[self._containerIndexes])


	def _totalIs(self, value):
		self._concentrations[:] = value


	def _countsInc(self, counts):
		assert (np.size(counts) == np.size(self._containerIndexes)) or np.size(counts) == 1, 'Inappropriately sized values'
		change_concentrations = self._state._counts_to_concentration(counts)

		self._state._concentrations[self._containerIndexes] += change_concentrations


	# Interface to Process
	def _totalConcentrations(self):
		return np.array(self._state._concentrations)[self._containerIndexes].copy()



class EnvironmentView(EnvironmentViewBase):
	def __init__(self, *args, **kwargs):
		super(EnvironmentView, self).__init__(*args, **kwargs)

		# State references
		assert len(set(self._query)) == len(self._query), "Environment views cannot contain duplicate entries"
		self._containerIndexes = self._state.container._namesToIndexes(self._query)


	def _dataSize(self):
		return len(self._query)


	def totalConcentrations(self):
		return self._totalConcentrations()


	def countsInc(self, counts):
		if self._state._infinite_environment:
			return

		self._countsInc(counts)
