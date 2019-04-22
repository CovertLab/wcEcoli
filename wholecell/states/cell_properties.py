
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
		self.waterIndex = sim_data.submassNameToIndex["water"]
		self.cell_density = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		# initialize properties
		# TODO -- these need to be initialized to the correct value
		self.cell_mass = 1339.0
		self.dry_mass = 403.0
		self.volume = self.cell_mass / self.cell_density
		# self.counts_to_molar = 1.0  # (units.g * units.mol / units.L)
		self.counts_to_molar = 1.0 / (self.nAvogadro * self.volume)

		# initialize properties IDs and values
		self.property_ids = [
			'cell_mass',
			'dry_mass',
			'volume',
			'counts_to_molar',
			 ]
		self.property_values = np.array([
			self.cell_mass,
			self.dry_mass,
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
		submasses = postEvolveMasses.sum(axis=0)  # sum over the processes

		self.water_mass = submasses[self.waterIndex]
		self.dry_mass = self.cell_mass - self.water_mass

		self.volume = self.cell_mass / self.cell_density
		self.counts_to_molar = 1 / (self.nAvogadro * self.volume)

		self.property_values = np.array([
			self.cell_mass,
			self.dry_mass,
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

	def __init__(self, state, process, query): # weight, priority, coupling id, option to not evaluate the query
		self._state = state
		self._state.viewAdd(self)
		self._processId = process.name()
		self._processIndex = process._processIndex
		self._query = query
		self._values = np.zeros(self._dataSize(), np.float64) # number of objects that satisfy the query

	# Interface to State
	def _updateQuery(self):
		self._valueIs(self._state.container._counts[self._containerIndexes])

	def _valueIs(self, value):
		self._values[:] = value

	def _countsInc(self, counts):
		return

	# Interface to Process
	def _allValues(self):
		return np.array(self._state.property_values)[self._containerIndexes].copy()


class CellPropertiesView(CellPropertiesViewBase):
	def __init__(self, *args, **kwargs):
		super(CellPropertiesView, self).__init__(*args, **kwargs)

		# State references
		assert len(set(self._query)) == len(self._query), "cell property views cannot contain duplicate entries"
		self._containerIndexes = self._state.container._namesToIndexes(self._query)

	def _dataSize(self):
		return len(self._query)

	def allValues(self):
		return self._allValues()

	# def countsInc(self, molecule_ids, counts):
	# 	# self._state._environment_deltas = counts
	# 	# self._state.accumulate_deltas(molecule_ids, counts)
	# 	return
