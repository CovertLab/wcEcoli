
from __future__ import absolute_import, division, print_function

import os
import json
import struct

import numpy as np

from wholecell.utils import filepath


__all__ = [
	"TableWriter",
	# "TableWriterError",
	# "MissingFieldError",
	# "UnrecognizedFieldError",
	# "AttributeAlreadyExistsError",
	# "AttributeTypeError",
	# "VariableEntrySizeError",
	# "DtypeTooComplexError",
	]

VERSION = 3  # should update this any time there is a spec-breaking change

FILE_ATTRIBUTES = "attributes.json"
COLUMN_SIGNATURE = 0xDECAFF

# Column data file's header struct. See HEADER.pack(), below.
HEADER, DTYPE_BYTE_LEN = struct.Struct('<2I 3H 64p'), 64


class TableWriterError(Exception):
	"""
	Base exception class for TableWriter-associated exceptions.
	"""
	pass

class MissingFieldError(TableWriterError):
	"""
	An error raised when TableWriter.append is called without providing all
	field names.
	"""
	pass

class UnrecognizedFieldError(TableWriterError):
	"""
	An error raised when TableWriter.append is called with an unexpected field
	name.
	"""
	pass

class AttributeAlreadyExistsError(TableWriterError):
	"""
	An error raised when TableWriter.writeAttributes is called with an
	attribute name that has already been used.
	"""
	pass

class AttributeTypeError(TableWriterError):
	"""
	An error raised when TableWriter.writeAttributes is called with a type that
	does not appear to be JSON-serializable.
	"""
	pass

class VariableEntrySizeError(TableWriterError):
	"""Error raised on attempt to write an entry that's not the same size as
	the previous entries in the same Column.
	"""
	pass

class DtypeTooComplexError(TableWriterError):
	"""Error raised on attempt to write an entry whose NumPy dtype description
	doesn't fit in the relevant HEADER field. This is unexpected since all the
	current descriptions fit in 6 bytes including the p-string length byte:
	`"<i8"`, `"<f8"`, `"|b1"`, `"|S7"`.
	"""
	pass


class _Column(object):
	"""
	Manages the written data for a specific TableWriter field.

	Each field in a 'table' corresponds to a 'column' that is written to on
	each append operation.  This private class encapsulates the logic and data
	for a particular 'column'.

	Parameters:
		path (str): The path for this particular column's data file.

	Notes
	-----
	See TableWriter for more information about output file and directory
	structure.  Each column creates one file containing the packed HEADER
	followed by the array entries.

	TODO (John): With some adjustment this class could be made public, and used
		as a lightweight alternative to TableWriter in addition to part of
		TableWriter's internal implementation.
	"""

	def __init__(self, path):
		self._path = path
		self._data = open(path, "wb")
		self._dtype = None
		self._bytes_per_entry = None
		self._elements_per_entry = None  # aka subcolumn count


	def append(self, value):
		"""
		Appends an array-like entry to the end of a column, converting it
		to a 1-D array.

		The first call to this method will define the column's NumPy dtype,
		element array size (subcolumns), and element size in bytes.
		Subsequent entries must be consistent.

		Parameters:
			value (array-like): A NumPy ndarray or anything that can be cast to
				an array via np.asarray, including a scalar (i.e. 0D array).
		"""

		value = np.asarray(value, self._dtype)
		data_bytes = value.tobytes()

		if self._dtype is None:
			self._dtype = value.dtype

			descr = self._dtype.descr
			if len(descr) == 1 and descr[0][0] == "":
				descr = descr[0][1]
			descr_json = json.dumps(descr, separators=(',', ':'))

			self._bytes_per_entry = len(data_bytes)
			self._elements_per_entry = value.size

			if len(descr_json) >= DTYPE_BYTE_LEN:
				raise DtypeTooComplexError(descr_json)
			header = HEADER.pack(
				COLUMN_SIGNATURE,          # I: magic signature
				self._bytes_per_entry,     # I: bytes/entry
				self._elements_per_entry,  # H: subcolumns
				1,                         # H: entries/written block
				0,                         # H: compression type code (0 = none)
				descr_json,                # 64p: element dtype description
				)

			self._data.write(header)

		elif self._bytes_per_entry != len(data_bytes):
			raise VariableEntrySizeError(
				'Entry size in bytes, elements {} is inconsistent with {}'
				' for Table column {}'.format(
					(len(data_bytes), value.size),
					(self._bytes_per_entry, self._elements_per_entry),
					self._path))

		self._data.write(data_bytes)


	def close(self):
		"""
		Close the column's data file.

		Notes
		-----
		Trying to append after closing will raise an error.
		"""

		self._data.close()


	def __del__(self):
		"""
		Explicitly closes the output file once the instance is totally
		dereferenced.

		Notes
		-----
		This will lead to errors consequent of operating on a closed file if
		references to the output file (which ought to be private)
		exist.  This is desirable, because such references should not be made.
		"""
		self.close()


class TableWriter(object):
	"""
	Generic live streaming output writer for NumPy ndarrays.

	NumPy can save and load arrays to disk.  This class provides a convenient
	interface to repeated appending of NumPy array data to an output file,
	which can be loaded one column at a time via the TableReader class.

	A Table has one or more named columns. Each (row x column) entry is a
	1-D NumPy array.

	Output file structure:

	<root directory> : Root path, provided by during instantiation.
		/attributes.json : A JSON file containing the attributes and metadata.
		/<column name> : A file per column containing a HEADER followed
				by the entries which are binary data from NumPy ndarrays.

	Parameters:
		path (str): Path to the directory to create.  All data will be saved
			within this directory.  If the directory already exists, it won't
			raise an error, but don't put multiple Tables in the same directory.

	See also
	--------
	wholecell.io.tablereader.TableReader

	Notes
	-----
	The terms used in this class are adapted from the analogy of a spreadsheet
	or table.

	Each append() operation adds a new 'row' to the table containing one
	'entry' (or table cell) per 'column'.  Each 'column' corresponds to a
	'field' and holds a fixed array size and data type for all of its entries.
	append() will convert scalar values and higher dimension arrays to 1D
	arrays.  The 1D array elements are called 'subcolumns'.

	Both simple and structured ndarray dtypes are supported.  Structured arrays
	are a good way to work around the dimensional limitations.

	NumPy object arrays (e.g. arrays not of one or more standard dtypes) are
	not supported.  They might work, but under the hood will rely on pickling
	for saving and loading, which will be terribly slow.

	A Table also stores named 'attributes' meant for user-provided annotations
	that may be useful for downstream analysis or portability.  E.g. a list of
	element names for a vector-column.  Each attribute can be written only
	once -- attributes do not support time-series data.

	TODO (John): Test portability across machines (particularly, different
		operating systems).

	TODO (John): Consider writing all fields simultaneously
		as part of a structured array (i.e. a hybrid data type).
	"""

	def __init__(self, path):
		self._path = filepath.makedirs(path)
		self._columns = None

		self._attributes = {}
		self._attributes_filename = os.path.join(path, FILE_ATTRIBUTES)
		self.writeAttributes(_version=VERSION)


	def append(self, **namesAndValues):
		"""
		Write a new row of values, writing a 1-D NumPy array to each named
		column.

		The first call to this method will define the column names and dtypes.
		Subsequent calls will validate the names and types for consistency.

		Parameters:
			**namesAndValues (dict[str, array-like]):  The column names (fields)
				and associated values to append to the end of the columns.

		Notes
		-----
		All fields must be provided every time this method is called.
		"""

		if self._columns is None:
			self._columns = {
				name:_Column(os.path.join(self._path, name))
				for name in namesAndValues.viewkeys()
				}

		else:
			missingFields = self._columns.viewkeys() - namesAndValues.viewkeys()
			unrecognizedFields = namesAndValues.viewkeys() - self._columns.viewkeys()

			if missingFields:
				raise MissingFieldError(
					"Missing fields: {}".format(", ".join(missingFields))
					)

			if unrecognizedFields:
				raise UnrecognizedFieldError(
					"Unrecognized fields: {}".format(", ".join(unrecognizedFields))
					)

		for name, value in namesAndValues.viewitems():
			self._columns[name].append(value)


	def writeAttributes(self, **namesAndValues):
		"""
		Writes JSON-serializable data.

		Oftentimes some additional data is needed to contextualize the data in
		the columns.  This method can be called to write JSON-serializable
		data (e.g. a list of strings) alongside the column data.

		Parameters:
			**namesAndValues (dict[str, JSON-serializable]): The named
				attribute values.

		NOTE: TableWriter uses attribute names starting with "_" to hold its
		internal metadata.

		Notes
		-----
		This method can be called at any time, so long as an attribute name is
		not reused. That's because attributes are not designed for time-series
		data. [We could almost require all writeAttributes() to happen before
		append() except 'lengthSec' and 'endTime' attributes get written to
		'Main' at the end of each generation.]
		"""

		for name, value in namesAndValues.viewitems():
			if name in self._attributes:
				raise AttributeAlreadyExistsError(
					"An attribute named '{}' already exists.".format(name)
					)

			try:
				if isinstance(value, np.ndarray):
					print("Warning - converting '{}' attribute from ndarray to list for JSON serialization.".format(name))
					value = value.tolist()

				json.dumps(value)  # test that it's JSON serializable

			except TypeError:
				raise AttributeTypeError(
					"Attribute '{}' value ({!r}) was not JSON serializable.".format(name, value)
					)

			self._attributes[name] = value

		filepath.write_json_file(self._attributes_filename, self._attributes, indent=1)

	def close(self):
		"""
		Close the output files (columns).

		Notes
		-----
		Trying to append after closing will raise an error.

		"""

		if self._columns is not None:
			for column in self._columns.viewvalues():
				column.close()


	def __del__(self):
		"""
		Close the output files once the instance is totally dereferenced.
		"""
		self.close()
