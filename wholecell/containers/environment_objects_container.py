'''

environment_objects_container.py

'''

from __future__ import absolute_import

import numpy as np


class EnvironmentObjectsContainer(object):
	'''
	EnvironmentObjectsContainer

	A wrapper around a NumPy array that tracks the concentrations of environment objects.
	Environment objects are bulk molecules in the environment.
	'''

	def __init__(self, objectNames, dtype = np.float64):
		self._objectNames = list(objectNames)

		self._nObjects = len(self._objectNames)

		self._objectIndex = {objectName:index for index, objectName in enumerate(self._objectNames)}

		self._concentrations = np.zeros(len(self._objectNames), dtype)

		self._volume = 'infinite'


	def volume(self):
		return self._volume


	def volumeIs(self, value):
		self._volume = value


	def concentrations(self, names = None):
		if names is None:
			return self._concentrations.copy()

		else:
			return self._concentrations[self._namesToIndexes(names)]


	def concentrationsIs(self, values, names = None):
		if names is None:
			self._concentrations[:] = values

		else:
			self._concentrations[self._namesToIndexes(names)] = values


	def concentrationsView(self, names = None):
		if names is None:
			return _EnvironmentObjectsView(self, np.arange(self._nObjects))

		else:
			return _EnvironmentObjectsView(self, self._namesToIndexes(names))


	def concentration(self, name):
		return self._concentrations[self._objectIndex[name]]


	def concentrationIs(self, value, name):
		self._concentrations[self._objectIndex[name]] = value


	def concentrationView(self, name):
		return _EnvironmentObjectView(self, self._objectIndex[name])

	def objectNames(self):
		return tuple(self._objectNames)

	def emptyLike(self):
		names = self.objectNames()
		new_copy = EnvironmentObjectsContainer(names)
		return new_copy

	def _namesToIndexes(self, names):
		return np.array([self._objectIndex[name] for name in names])


	def __eq__(self, other):
		return (self._concentrations == other._concentrations).all()


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			objectNames = self._objectNames
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			concentrations = self._concentrations
			)


	def tableLoad(self, tableReader, tableIndex):
		assert self._objectNames == tableReader.readAttribute("objectNames")

		self._concentrations = tableReader.readRow(tableIndex)["concentrations"]



class _EnvironmentObjectsView(object):
	'''
	_EnvironmentObjectsView

	An accessor for a subset of objects in a EnvironmentObjectsContainer.
	'''

	def __init__(self, container, indexes):
		self._container = container
		self._indexes = indexes


	def concentrations(self):
		return self._container._concentrations[self._indexes]


	def concentrationsIs(self, values):
		self._container._concentrations[self._indexes] = values



class _EnvironmentObjectView(object):
	'''
	_EnvironmentObjectView

	An accessor for a single object in a EnvironmentObjectsContainer.
	'''

	def __init__(self, container, index):
		self._container = container
		self._index = index


	def concentration(self):
		return self._container._concentrations[self._index]


	def concentrationIs(self, values):
		self._container._concentrations[self._index] = values
