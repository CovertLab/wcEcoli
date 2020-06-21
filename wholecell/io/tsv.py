"""
CSV reader and writer that handles Python 2/3 compatibility and defaults to TAB
delimiters.
"""

from __future__ import absolute_import, division, print_function

import csv
from io import TextIOWrapper
from typing import Any, IO, Iterator, Union

import six


def reader(csvfile, dialect='excel', delimiter='\t', **fmtparams):
	# type: (IO[bytes], Union[str, csv.Dialect], str, **Any) -> Iterator[list]
	"""Open a csv.reader(), handling Python 2/3 I/O compatibility, and defaults
	to TAB delimiters.

	REQUIRES: csvfile must be a buffered byte reader, e.g. from
	io.open(filename, 'rb') or io.BytesIO(buffer).

	NOTE: To support Unicode in PY2, this could be extended to decode the
	csv.reader's output bytes.
	"""
	if six.PY2:
		input_file = csvfile
	else:
		input_file = TextIOWrapper(csvfile, encoding='utf-8', newline='')
	return csv.reader(input_file, dialect=dialect, delimiter=delimiter, **fmtparams)


def writer(csvfile, dialect='excel', delimiter='\t', **fmtparams):
	# type: (IO[bytes], Union[str, csv.Dialect], str, **Any) -> Any
	"""Open a csv.writer(), handling Python 2/3 I/O compatibility and defaults
	to TAB delimiters.

	REQUIRES: csvfile must be a buffered byte writer, e.g. from
	io.open(filename, 'wb') or io.BytesIO(buffer).

	NOTE: To support Unicode in PY2, this could be extended to encode the
	csv.writer's input bytes.
	"""
	if six.PY2:
		input_file = csvfile
	else:
		input_file = TextIOWrapper(
			csvfile, encoding='utf-8', newline='', line_buffering=True)
	return csv.writer(input_file, dialect=dialect, delimiter=delimiter, **fmtparams)
