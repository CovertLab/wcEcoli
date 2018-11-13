
from __future__ import absolute_import, division, print_function

import os
import json
import numpy as np
import struct
import zlib

from wholecell.utils import filepath


__all__ = [
	"TableWriter",
	# "TableWriterError",
	# "MissingFieldError",
	# "UnrecognizedFieldError",
	# "AttributeAlreadyExistsError",
	# "AttributeTypeError",
	# "VariableEntrySizeError",
	]

VERSION = 3  # should update this any time there is a spec-breaking change

FILE_ATTRIBUTES = "attributes.json"

# Chunk type and size.
CHUNK_HEADER = struct.Struct('>4s I')
COLUMN_CHUNK_TYPE = 'COLM'  # column file's header chunk
BLOCK_CHUNK_TYPE = 'BLOC'   # data block chunk

# Column header struct. See the pack() calls for field details.
COLUMN_STRUCT = struct.Struct('>2I 2H')

COMPRESSION_TYPE_NONE = 0
COMPRESSION_TYPE_ZLIB = 1

# zlib's default compression level is 6.
#
# Measuring zlib compression of Mass/time, BulkMolecules/atpRequested, and
# BulkMolecules/counts, higher levels yielded a fraction of a percent space
# savings at substantial time costs. For BulkMolecules/counts, level 7 adds
# ~40% time, level 8 takes 5x, and level 9 takes 15x time.
ZLIB_LEVEL = 6

# Pack enough entries into each data block to total about this many bytes
# before compression.
#
# Compressing multiple entries as a block saves considerable space for small
# and medium size entries. E.g. compressing an 8-byte Mass/time entry doubles
# the size, and adding a chunk header makes it 3x, while compressing 512 of
# them into a block gets the entire chunk down to 68% of the input size.
#
# BulkMolecules/atpRequested has 96-byte entries which compress individually to
# 30% including header. It goes down to 10% after packing 16 entries together,
# 9% when packing 32 entries together, and about 7% when packing 512 entries
# together.
#
# Packing also saves compression time and should save some I/O time.
# TODO: Measure I/O time at different block sizes.
#
# The column file format isn't designed for random access but a reader could
# skip chunk to chunk without decompressing them to get to a desired block, and
# adding a chunk index would enable block-level random access.
#
# zlib supports incrementally compressing a bytestring at a time to save RAM
# and use the cumulative data history to improve the overall compression ratio,
# but measurements show it can't obviate data blocks. For readers to decompress
# an entry at a time, the writer has to do a Z_SYNC_FLUSH which reduces the
# compression ratio (even if the writer discards the '\x00\x00\xff\xff' suffix
# and the reader restores it). So the writer still needs to group entries into
# blocks. Then cumulative incremental compression has a tiny +/- impact on
# compression size, and the ability to read an entry without decompressing
# everything before it depends on saving & reloading the compression state,
# which takes space and that feature is not in the Python zlib library.
BLOCK_BYTES_GOAL = 4096

V2_DIR_COLUMNS = "columns"  # format v2's directory of column files


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

class TableExitsError(TableWriterError):
	"""Error raised on attempt to create a Table in a directory that already
	has a table (completed or in progress).
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
	structure.  Each column creates one file containing a COLUMN_CHUNK_TYPE
	chunk followed by the data BLOCK_CHUNK_TYPE chunks. Each data block
	contains one or more NumPy array 'entries', optionally compressed. One
	can read it using the Python chunk library Chunk(file, False).

	TODO (John): With some adjustment this class could be made public, and used
		as a lightweight alternative to TableWriter in addition to part of
		TableWriter's internal implementation.
	"""

	def __init__(self, path, compression_type=COMPRESSION_TYPE_ZLIB):
		if compression_type not in (COMPRESSION_TYPE_NONE, COMPRESSION_TYPE_ZLIB):
			raise ValueError('Unknown compression type {}'.format(compression_type))

		self._path = path
		self._data = open(path, "wb")
		self._dtype = None
		self._entries_per_block = None
		self._bytes_per_entry = None
		self._elements_per_entry = None  # aka subcolumn count
		self._compression_type = compression_type
		self._current_block = []


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

		# First entry: Write the column header.
		if self._dtype is None:
			self._dtype = value.dtype
			self._bytes_per_entry = value.nbytes
			self._entries_per_block = int(np.ceil(BLOCK_BYTES_GOAL / value.nbytes))
			self._elements_per_entry = value.size

			descr = self._dtype.descr
			if len(descr) == 1 and descr[0][0] == "":
				descr = descr[0][1]
			descr_json = json.dumps(descr, separators=(',', ':'))

			chunk_header = CHUNK_HEADER.pack(
				COLUMN_CHUNK_TYPE,
				COLUMN_STRUCT.size + len(descr_json),
				)
			column_struct = COLUMN_STRUCT.pack(
				self._bytes_per_entry,     # I: bytes/entry, before packing, compression, and CHUNK_HEADER
				self._elements_per_entry,  # I: subcolumns
				self._entries_per_block,   # H: packed entries/block; the last block may be smaller
				self._compression_type,    # H: compression type code
				)
			self._data.write(chunk_header + column_struct + descr_json)

		# Later entry: Check consistency.
		elif self._bytes_per_entry != len(data_bytes):
			raise VariableEntrySizeError(
				'Entry size in bytes, elements {} is inconsistent with {}'
				' for Table column {}'.format(
					(len(data_bytes), value.size),
					(self._bytes_per_entry, self._elements_per_entry),
					self._path))

		# Collect up an I/O block. Compress and write it when it's full.
		self._current_block.append(data_bytes)
		if len(self._current_block) >= self._entries_per_block:
			self._write_block()
		elif self._data.closed:
			raise ValueError('I/O operation on closed file')


	def _write_block(self):
		'''Compress and write the current block, if any.'''
		if self._current_block:
			block_data = b''.join(self._current_block)

			if self._compression_type == COMPRESSION_TYPE_ZLIB:
				block_data = zlib.compress(block_data, ZLIB_LEVEL)

			block_header = CHUNK_HEADER.pack(BLOCK_CHUNK_TYPE, len(block_data))
			self._data.write(block_header + block_data)

			del self._current_block[:]


	def close(self):
		"""
		Finish writing and close the column's data file. Idempotent.

		Notes
		-----
		Trying to append after closing will raise an error.
		"""

		if not self._data.closed:
			try:
				self._write_block()
				self._data.truncate()
			finally:
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
		/<column name> : A file per column containing a COLUMN_CHUNK_TYPE of
				header chunk followed by data chunks containing compressed
				NumPy ndarray data. It's extensible.

	Parameters:
		path (str): Path to the directory to create.  All data will be saved
			within this directory.  It's OK if the directory already exists but
			not OK if it contains Table files. This will raise TableExitsError
			if its attributes.json file already exists. Existing or concurrent
			Table files would confuse each other.

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

		if os.path.exists(self._attributes_filename) or os.path.exists(V2_DIR_COLUMNS):
			raise TableExitsError('In {}'.format(self._path))

		# The column file's magic signature mostly obviates the '_version'
		# attribute but writing the attributes file now lets the above check
		# also prevent competing TableWriters.
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
