"""Unit test for the tsv module."""
from __future__ import absolute_import, division, print_function

from io import BytesIO
import six
from six.moves import zip_longest
import unittest

from wholecell.io import tsv


FIELD_NAMES = ['id', 'ourLocation', u'\u20ac:xyz', 'mass (units.g)']
INPUT_ROWS = [
	b'id\tourLocation\t\xE2\x82\xAC:xyz\tmass (units.g)',
	b'G6660-MONOMER\t[x]\tLocation information from Lopez Campistrous 2005.\t98.6',
	b'2.71828\t[c]\tLocation from \xe2\x8a\x972011.\t12']


class Test_Tsv(unittest.TestCase):
	def test_reader(self):
		byte_stream = BytesIO(b'\n'.join(INPUT_ROWS))
		reader = tsv.reader(byte_stream)

		for row, expected in zip_longest(reader, INPUT_ROWS, fillvalue=404):
			if six.PY2:
				# tsv.reader doesn't support unicode in PY2.
				assert row == expected.split(b'\t')
			else:
				assert row == expected.decode('utf-8').split('\t')

	def test_writer(self):
		byte_stream = BytesIO()
		writer = tsv.writer(byte_stream)
		for row in INPUT_ROWS:
			if six.PY2:
				# tsv.writer doesn't support unicode in PY2.
				writer.writerow(row.split(b'\t'))
			else:
				writer.writerow(row.decode('utf-8').split('\t'))

		data = byte_stream.getvalue()
		assert data == b'\r\n'.join(INPUT_ROWS + [b''])
