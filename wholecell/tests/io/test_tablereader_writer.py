from __future__ import absolute_import, division, print_function

import numpy as np
import numpy.testing as npt
import shutil
import tempfile
import unittest

import nose.plugins.attrib as noseAttrib


def noop_decorator():
	def same(fcn):
		return fcn
	return same

# TODO(jerry): Remove before flight: Workaround for @profile defined only by kernprof.
__builtins__.setdefault('profile', noop_decorator())


from wholecell.io.tablereader import TableReader, DoesNotExistError, VariableWidthError
from wholecell.io.tablewriter import TableWriter, MissingFieldError, UnrecognizedFieldError


COLUMNS = 'x y z theta'.split()
DATA = {key: np.arange(10.0) + ord(key[0]) for key in COLUMNS}


# TODO(jerry): Test readColumn() w/non-default indices and block_read.
# TODO(jerry): Test structured dtypes.
# TODO(jerry): Test or delete iterColumn().

class Test_TableReader_Writer(unittest.TestCase):
	def setUp(self):
		self.test_dir = None

	def tearDown(self):
		if self.test_dir:
			shutil.rmtree(self.test_dir)

	def make_test_dir(self):
		if not self.test_dir:
			self.test_dir = tempfile.mkdtemp()

	def assert_equal_row(self, written_row, read_row):
		'''Assert that a read-in row (dict of cell values) matches the written
		row, adjusting for the fact that **TableWriter stores a 1-D array in
		each table cell even if the written cell value was a scalar or an
		n-D array.**
		'''
		self.assertEqual(set(written_row.keys()), set(read_row.keys()))

		for key, expected_value in written_row.iteritems():
			actual_value = read_row[key]
			if isinstance(expected_value, np.ndarray):
				expected_value = expected_value.reshape(expected_value.size)
			self.assertIsInstance(actual_value, np.ndarray)
			npt.assert_array_equal(expected_value, actual_value)

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
		self.assert_equal_row(d2, actual)

		actual = reader.readRow(0)
		self.assert_equal_row(DATA, actual)

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
		writer = TableWriter(self.test_dir)
		writer.append(**DATA)
		writer.append(**d0)  # row 1
		writer.append(**d2)
		writer.append(**d3)
		writer.close()

		# --- Read ---
		reader = TableReader(self.test_dir)
		self.assertEqual(set(COLUMNS), set(reader.columnNames()))

		column_name = COLUMNS[0]
		with self.assertRaises(VariableWidthError):
			reader.readColumn(column_name)

		actual = reader.readRow(1)
		self.assertEqual(d0, actual)
		self.assert_equal_row(d0, actual)
		for key in DATA:  # expect 1-length arrays
			self.assertEqual(np.float64, actual[key].dtype)
			self.assertEqual((1,), actual[key].shape)

		actual = reader.readRow(2)
		self.assert_equal_row(d2, actual)

		actual = reader.readRow(3)
		self.assert_equal_row(d3, actual)

		actual = reader.readRow(0)
		self.assert_equal_row(DATA, actual)
