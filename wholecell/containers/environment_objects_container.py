'''

environment_objects_container.py

'''

from __future__ import division

import numpy as np


class EnvironmentObjectsContainer(object):
	'''
	EnvironmentObjectsContainer

	A wrapper around a NumPy array that tracks the counts of environmental objects.
	Environmental objects are those external to the whole cells.
	'''

	def __init__(self, objectNames, dtype = np.int64):
		self._objectNames = list(objectNames)

		self._nObjects = len(self._objectNames)

		self._objectIndex = {objectName:index for index, objectName in enumerate(self._objectNames)}

		self._counts = np.zeros(len(self._objectNames), dtype)


	def counts(self, names = None):
		if names is None:
			return self._counts.copy()

		else:
			return self._counts[self._namesToIndexes(names)]


	def setCounts(self, values, names = None):
		if names is None:
			self._counts[:] = values

		else:
			self._counts[self._namesToIndexes(names)] = values


	def increaseCounts(self, values, names = None):
		if names is None:
			self._counts[:] += values

		else:
			self._counts[self._namesToIndexes(names)] += values


	def decreaseCounts(self, values, names = None): # TODO: raise exception if > max?
		if names is None:
			self._counts[:] -= values

		else:
			self._counts[self._namesToIndexes(names)] -= values


	def countsView(self, names = None):
		if names is None:
			return _EnvironmentObjectsView(self, np.arange(self._nObjects))

		else:
			return _EnvironmentObjectsView(self, self._namesToIndexes(names))


	def count(self, name):
		return self._counts[self._objectIndex[name]]


	def setCount(self, value, name):
		self._counts[self._objectIndex[name]] = value


	def increaseCount(self, value, name):
		self._counts[self._objectIndex[name]] += value


	def decreaseCount(self, value, name):
		self._counts[self._objectIndex[name]] -= value


	def countView(self, name):
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
		return (self._counts == other._counts).all()


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			objectNames = self._objectNames
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			counts = self._counts
			)


	def tableLoad(self, tableReader, tableIndex):
		assert self._objectNames == tableReader.readAttribute("objectNames")

		self._counts = tableReader.readRow(tableIndex)["counts"]


class _EnvironmentObjectsView(object):
	'''
	_EnvironmentObjectsView

	An accessor for a subset of objects in a EnvironmentObjectsContainer.
	'''

	def __init__(self, container, indexes):
		self._container = container
		self._indexes = indexes


	def counts(self):
		return self._container._counts[self._indexes]


	def setCounts(self, values):
		self._container._counts[self._indexes] = values


	def increaseCounts(self, values):
		self._container._counts[self._indexes] += values


	def decreaseCounts(self, values):
		self._container._counts[self._indexes] -= values


class _EnvironmentObjectView(object):
	'''
	_EnvironmentObjectView

	An accessor for a single object in a EnvironmentObjectsContainer.
	'''

	def __init__(self, container, index):
		self._container = container
		self._index = index


	def count(self):
		return self._container._counts[self._index]


	def setCount(self, values):
		self._container._counts[self._index] = values


	def increaseCount(self, values):
		self._container._counts[self._index] += values


	def decreaseCount(self, values):
		self._container._counts[self._index] -= values