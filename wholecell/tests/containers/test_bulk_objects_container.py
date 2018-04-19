'''
test_bulk_objects_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.containers.bulk_objects_container import BulkObjectsContainer

OBJECT_NAMES = ('ATP', 'glucose', 'glycine')
OBJECT_COUNTS = [100, 20, 10]


class Test_BulkObjectsContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = createContainer()


	def tearDown(self):
		pass

	# Interface methods

	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_counts(self):
		self.assertEqual(
			self.container.counts().tolist(),
			OBJECT_COUNTS
			)

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			OBJECT_COUNTS[1:]
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsIs(self):
		newCounts = [10, 20.0, 30]
		self.container.countsIs(newCounts)

		self.assertEqual(
			self.container.counts().tolist(),
			newCounts
			)

		newCounts = [15, 5.0]
		self.container.countsIs(newCounts, OBJECT_NAMES[1:])

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			newCounts
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsInc(self):
		incCounts = [10, 20.0, 30]
		newCounts = (self.container.counts() + np.array(incCounts)).tolist()
		self.container.countsInc(incCounts)

		self.assertEqual(
			self.container.counts().tolist(),
			newCounts
			)

		incCounts = [15, 5.0]
		newCounts = (self.container.counts(OBJECT_NAMES[1:]) + np.array(incCounts)).tolist()
		self.container.countsInc(incCounts, OBJECT_NAMES[1:])

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			newCounts
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsDec(self):
		decCounts = [30, 10.0, 5]
		newCounts = (np.array(OBJECT_COUNTS) - np.array(decCounts)).tolist()
		self.container.countsDec(decCounts)

		self.assertEqual(
			self.container.counts().tolist(),
			newCounts
			)

		decCounts = [5.0, 2]
		newCounts = (self.container.counts(OBJECT_NAMES[1:]) - np.array(decCounts)).tolist()
		self.container.countsDec(decCounts, OBJECT_NAMES[1:])

		self.assertEqual(
			self.container.counts(OBJECT_NAMES[1:]).tolist(),
			newCounts
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsView_countsInc(self):
		countsView = self.container.countsView()
		incCounts = [10, 20.0, 30]
		newCounts = (self.container.counts() + np.array(incCounts)).tolist()
		countsView.countsInc(incCounts)

		self.assertEqual(countsView.counts().tolist(), newCounts)

		countsView = self.container.countsView(OBJECT_NAMES[1:])
		incCounts = [15, 5.0]
		newCounts = (self.container.counts(OBJECT_NAMES[1:]) + np.array(incCounts)).tolist()
		countsView.countsInc(incCounts)

		self.assertEqual(countsView.counts().tolist(), newCounts)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countsView_countsDec(self):
		countsView = self.container.countsView()
		decCounts = [30, 10.0, 5]
		newCounts = (np.array(OBJECT_COUNTS) - np.array(decCounts)).tolist()
		countsView.countsDec(decCounts)

		self.assertEqual(countsView.counts().tolist(), newCounts)

		countsView = self.container.countsView(OBJECT_NAMES[1:])
		decCounts = [5.0, 2]
		newCounts = (self.container.counts(OBJECT_NAMES[1:]) - np.array(decCounts)).tolist()
		countsView.countsDec(decCounts)

		self.assertEqual(countsView.counts().tolist(), newCounts)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countView(self):
		countView = self.container.countView(OBJECT_NAMES[2])
		incCount = 333.0
		newCount = OBJECT_COUNTS[2] + incCount
		countView.countInc(incCount)

		self.assertEqual(countView.count(), newCount)

		newCount -= incCount
		countView.countDec(incCount)
		self.assertEqual(countView.count(), newCount)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_count(self):
		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			OBJECT_COUNTS[0]
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countIs(self):
		newCount = 40
		self.container.countIs(newCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countInc(self):
		incCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)

		incCount = 40.0
		newCount = self.container.count(OBJECT_NAMES[0]) + incCount
		self.container.countInc(incCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_countDec(self):
		decCount = 40
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)

		decCount = 40.0
		newCount = self.container.count(OBJECT_NAMES[0]) - decCount
		self.container.countDec(decCount, OBJECT_NAMES[0])

		self.assertEqual(
			self.container.count(OBJECT_NAMES[0]),
			newCount
			)


	# Internal methods


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_namesToIndexes(self):
		# Test normal ordering
		names = OBJECT_NAMES
		self.assertEqual(
			self.container._namesToIndexes(names).tolist(),
			np.arange(len(OBJECT_NAMES)).tolist()
			)

		# Test reverse ordering
		names = OBJECT_NAMES[::-1]
		self.assertEqual(
			self.container._namesToIndexes(names).tolist(),
			np.arange(len(OBJECT_NAMES))[::-1].tolist()
			)

		# Test subset
		names = OBJECT_NAMES[::2]
		self.assertEqual(
			self.container._namesToIndexes(names).tolist(),
			np.arange(len(OBJECT_NAMES))[::2].tolist()
			)


	@noseAttrib.attr('smalltest', 'bulkObjects')
	def test_eq(self):
		newContainer = createContainer()

		self.assertEqual(
			self.container.counts().tolist(),
			newContainer.counts().tolist()
			)

# TODO: view tests


def createContainer():
	container = BulkObjectsContainer(OBJECT_NAMES)

	container.countsIs(OBJECT_COUNTS)

	return container
