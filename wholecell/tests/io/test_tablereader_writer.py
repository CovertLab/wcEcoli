from __future__ import absolute_import, division, print_function

import numpy as np
import numpy.testing as npt
import os
import shutil
import tempfile
import unittest

import nose.plugins.attrib as noseAttrib


def noop_decorator(fcn):
	return fcn

# TODO(jerry): Remove before flight: Workaround for @profile defined only by kernprof.
__builtins__.setdefault('profile', noop_decorator)


from wholecell.io.tablereader import TableReader, DoesNotExistError
from wholecell.io.tablewriter import (TableWriter, MissingFieldError,
	UnrecognizedFieldError, VariableEntrySize)


COLUMNS = 'x y z theta'.split()
DATA = {key: np.arange(10.0) + ord(key[0]) for key in COLUMNS}


# TODO(jerry): Test readColumn() w/non-default indices.
# TODO(jerry): Test structured dtypes.

class Test_TableReader_Writer(unittest.TestCase):
	def setUp(self):
		self.test_dir = None
		self.table_path = None

	def tearDown(self):
		if self.test_dir:
			shutil.rmtree(self.test_dir)

	def make_test_dir(self):
		'''Create a temp test dir. TableWriter must make its own subdir.'''
		if not self.test_dir:
			self.test_dir = tempfile.mkdtemp()
			self.table_path = os.path.join(self.test_dir, 'Main')

	@noseAttrib.attr('smalltest', 'table')
	def test_basic(self):
		'''Test a table with some float arrays and no attributes.'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
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
		reader = TableReader(self.table_path)
		self.assertEqual([], reader.attributeNames())
		self.assertEqual({'_version'}, set(reader.allAttributeNames()))
		self.assertEqual(set(COLUMNS), set(reader.columnNames()))

		with self.assertRaises(DoesNotExistError):
			reader.readAttribute('x')

		column_name = COLUMNS[0]
		actual = reader.readColumn(column_name)  # the basic readColumn() case
		expected = np.vstack((DATA[column_name], d2[column_name]))
		npt.assert_array_equal(expected, actual)
		self.assertEqual(2, actual.ndim)

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
		writer = TableWriter(self.table_path)
		writer.writeAttributes(**d1)
		writer.writeAttributes(**d2)
		writer.writeAttributes(**d3)
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
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

	@noseAttrib.attr('smalltest', 'table')
	def test_array_dimensions(self):
		'''Test 0D, 2D, and non-uniform array lengths, also automatic int to
		float value conversion.
		'''
		self.make_test_dir()
		d0 = {key: 19 for key in DATA.iterkeys()}
		d2 = {key: value.reshape(2, -1) for key, value in DATA.iteritems()}
		d3 = {key: value[2:] for key, value in DATA.iteritems()}

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(**DATA)  # 1-D float arrays in table row 0

		with self.assertRaises(VariableEntrySize):
			writer.append(**d0)  # 0-D arrays: inconsistent entry sizes

		writer.append(**d2)  # 2-D arrays in table row 2

		with self.assertRaises(VariableEntrySize):
			writer.append(**d3)  # narrower 1-D arrays than row 0

		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		self.assertEqual(set(COLUMNS), set(reader.columnNames()))

		column_name = COLUMNS[0]
		actual = reader.readColumn(column_name)
		self.assertEqual(1, actual[1].ndim)  # 2-D d2 values reshaped to 1-D
		expected = np.vstack((DATA[column_name], DATA[column_name]))
		npt.assert_array_equal(expected, actual)

	@noseAttrib.attr('smalltest', 'table')
	def test_scalars(self):
		'''Test the case where all rows contain scalars, where readColumn()
		returns a 1-D array instead of a 2-D array!
		'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(x=20)  # scalar int
		writer.append(x=21)
		writer.append(x=22)
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual = reader.readColumn('x')
		self.assertEqual(1, actual.ndim)
		self.assertEqual((3,), actual.shape)
		npt.assert_array_equal([20, 21, 22], actual)

	@noseAttrib.attr('smalltest', 'table')
	def test_1_element_arrays(self):
		'''Test the case where all rows contain 1-element arrays, where
		readColumn() returns a 1-D array instead of a 2-D array!
		'''
		self.make_test_dir()

		# --- Write ---
		writer = TableWriter(self.table_path)
		writer.append(x=[20])
		writer.append(x=[21])
		writer.append(x=[22])
		writer.close()

		# --- Read ---
		reader = TableReader(self.table_path)
		actual = reader.readColumn('x')
		self.assertEqual(1, actual.ndim)
		self.assertEqual((3,), actual.shape)
		npt.assert_array_equal([20, 21, 22], actual)

	@unittest.skip('sensitive to implementation details; maybe not portable')
	@noseAttrib.attr('smalltest', 'table')
	def test_path_clash(self):
		'''Test what happens if two TableWriters are writing to the same directory.'''
		self.make_test_dir()
		d1 = {key: np.arange(5) + ord(key[0]) for key in COLUMNS}

		writer1 = TableWriter(self.table_path)
		writer1.append(**d1)
		writer1.append(**d1)
		writer1.append(**d1)

		writer2 = TableWriter(self.table_path)
		writer2.append(**d1)

		d2 = {key: -value for key, value in d1.iteritems()}
		writer2.append(**d2)

		writer1.writeAttributes(x=1, y=2, z=3)
		writer2.writeAttributes(x=10, a='a', b='b')

		# WARNING: The results depend on the closing order, the longer files,
		# and probably I/O buffer sizes!
		writer1.close()
		writer2.close()

		reader = TableReader(self.table_path)
		column_name = COLUMNS[0]
		actual = reader.readColumn(column_name)
		expected = np.vstack((d1[column_name], d2[column_name], d1[column_name]))
		npt.assert_array_equal(expected, actual)

		# WARNING: The last table to write any attributes (including automatic
		# writes of internal Table metadata) takes all. This esp. matters if
		# the column names are stored in an attribute.
		self.assertEqual(set('x a b'.split()), set(reader.attributeNames()))
