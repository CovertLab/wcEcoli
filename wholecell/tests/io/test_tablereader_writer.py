from __future__ import absolute_import, division, print_function

import numpy as np
import numpy.testing as npt
import shutil
import tempfile
import unittest

import nose.plugins.attrib as noseAttrib


def noop_decorator():
	def noop(fcn):
		return fcn
	return noop

# Cope if this is not running under kernprof.
__builtins__.setdefault('profile', noop_decorator())


from wholecell.io.tablereader import TableReader, DoesNotExistError
from wholecell.io.tablewriter import TableWriter, MissingFieldError, UnrecognizedFieldError


COLUMNS = 'x y z theta'.split()
DATA = {key: np.arange(10.0) + ord(key[0]) for key in COLUMNS}


class Test_BulkObjectsContainer(unittest.TestCase):
	def setUp(self):
		self.test_dir = None

	def tearDown(self):
		if self.test_dir:
			shutil.rmtree(self.test_dir)

	def make_test_dir(self):
		if not self.test_dir:
			self.test_dir = tempfile.mkdtemp()

	@noseAttrib.attr('smalltest', 'table')
	def test_basic(self):
		'''Test a table with some float arrays and no attributes.'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.test_dir)
		writer.append(**DATA)

		d2 = {key: -10 * value for key, value in DATA.iteritems()}
		writer.append(**d2)

		no_theta = dict(d2)
		del no_theta['theta']
		with self.assertRaises(MissingFieldError):
			writer.append(**no_theta)

		with self.assertRaises(UnrecognizedFieldError):
			writer.append(JUNK=np.arange(4), **d2)

		writer.close()

		with self.assertRaises(ValueError):
			writer.append(**DATA)  # writer is closed

		# --- Read ---
		reader = TableReader(self.test_dir)
		self.assertEqual([], reader.attributeNames())
		self.assertEqual(set(COLUMNS), set(reader.columnNames()))

		with self.assertRaises(DoesNotExistError):
			reader.readAttribute('x')

		column_name = COLUMNS[0]
		actual = reader.readColumn(column_name)  # the basic readColumn() case
		expected = np.vstack((DATA[column_name], d2[column_name]))
		npt.assert_array_equal(expected, actual)

		actual = reader.readRow(1)
		self.assertEqual(set(d2.keys()), set(actual.keys()))
		for key in d2:
			npt.assert_array_equal(d2[key], actual[key])

		actual = reader.readRow(0)
		self.assertEqual(set(DATA.keys()), set(actual.keys()))
		for key in DATA:
			npt.assert_array_equal(DATA[key], actual[key])

		column_name = COLUMNS[2]
		values = list(reader.iterColumn(column_name))
		npt.assert_array_equal((DATA[column_name], d2[column_name]), values)

		with self.assertRaises(DoesNotExistError):
			reader.readColumn('JUNK')

		reader.close()


	@noseAttrib.attr('smalltest', 'table')
	def test_attributes(self):
		'''Test a table with attributes and no columns.'''
		def check_attributes(attribute_dict):
			for k, v in attribute_dict.iteritems():
				self.assertEqual(attribute_dict[k], reader.readAttribute(k))

		self.make_test_dir()

		d1 = dict(mercury=1, venus=2, earth=3, mars=4)
		d2 = dict(jupiter=[50, 60, 70], saturn='Saturn')
		d3 = dict(uranus=700.0, neptune=800.5)
		keys = set(d1.keys() + d2.keys() + d3.keys())

		# --- Write ---
		writer = TableWriter(self.test_dir)
		writer.writeAttributes(**d1)
		writer.writeAttributes(**d2)
		writer.writeAttributes(**d3)
		writer.close()

		# --- Read ---
		reader = TableReader(self.test_dir)
		self.assertEqual([], reader.columnNames())
		self.assertEqual(keys, set(reader.attributeNames()))
		self.assertEqual(len(keys), len(reader.attributeNames()))

		check_attributes(d1)
		check_attributes(d3)
		check_attributes(d2)
		check_attributes(d1)
		check_attributes(d2)

		with self.assertRaises(DoesNotExistError):
			reader.readColumn('JUNK')

		reader.close()

	# TODO(jerry): Test readColumn() w/non-default indices and block_read.
