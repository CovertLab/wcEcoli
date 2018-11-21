'''
Compare whether two simulation runs produced identical output and, if not,
print info on where they differ. This also returns a shell exit status code.

This aim is to check that a code change such as a simulation speedup didn't
accidentally change the output. This does not support a tolerance between nor
check for semantic equivalence.

Example command line:

    diff_simouts.py out/manual/wildtype_000000/000000/generation_000000/000000/simOut \
    	out/experiment/wildtype_000000/000000/generation_000000/000000/simOut
'''

from __future__ import absolute_import, division, print_function

import numpy as np
import os
from pprint import pprint
import sys

from wholecell.io.tablereader import TableReader, VersionError


def array_equal(array1, array2):
	'''Returns True if the two ndarrays are equal, checking the shape and all
	elements, allowing for NaN values.

	Args:
		array1 (np.ndarray): array to compare
		array2 (np.ndarray): array to compare
	'''
	try:
		np.testing.assert_array_equal(array1, array2)
		return True
	except AssertionError:
		return False

def diff_name_lists(diffs, kind, list1, list2):
	'''Diff name list `list2` against `list1`.

	Args:
		diffs (dict): a dict to accumulate info about the differences
		kind (str): the kind of things named in the lists, e.g. 'Table'
		list1 (list[str]): a list of things from one simOut dir
		list2 (list[str]): a list of things from the other simOut dir
	Returns:
		set[str]: a set of the names in common
	Side effects:
		adds into about the differing names to `diffs`
	'''
	set1 = set(list1)
	set2 = set(list2)
	missing = set1 - set2
	extra = set2 - set1

	if missing:
		diffs['-Missing: ' + kind] = sorted(missing)
	if extra:
		diffs['+Extra: ' + kind] = sorted(extra)

	return set1 & set2

def diff_tables(diffs, table_name, simout_dir1, simout_dir2):
	'''Diff the named Table in `simout_dir2` against the one in `simout_dir1`.

	Args:
		diffs (dict): a dict to accumulate info about the differences
		table_name (str): the name of a table in both simOut dirs
		simout_dir1 (str): the reference simOut directory
		simout_dir2 (str): the other simOut directory
	'''
	table_errors = []
	table1 = table2 = None
	try:
		table1 = TableReader(os.path.join(simout_dir1, table_name))
	except VersionError as e:
		table_errors.append("simOut dir 1 doesn't have a Table: " + str(e))

	try:
		table2 = TableReader(os.path.join(simout_dir2, table_name))
	except VersionError as e:
		table_errors.append("simOut dir 2 doesn't have a Table: " + str(e))

	if table_errors:
		diffs['Subdir: ' + table_name] = table_errors
		return

	# Compare Column names.
	column_names1 = table1.columnNames()
	column_names2 = table2.columnNames()
	table_diffs = {}
	common_columns = diff_name_lists(table_diffs, 'Column file', column_names1, column_names2)

	# Diff the Column contents in common_columns and add to table_diffs.
	diff_contents = [
		name for name in common_columns
		if not array_equal(table1.readColumn2D(name), table2.readColumn2D(name))
		]
	if diff_contents:
		table_diffs['Unequal Columns'] = sorted(diff_contents)

	# Diff the Table Attributes and add to table_diffs.
	attribute_names1 = table1.attributeNames()
	attribute_names2 = table2.attributeNames()
	common_attributes = diff_name_lists(table_diffs, 'Attribute', attribute_names1, attribute_names2)

	diff_attributes = [
		key for key in common_attributes
		if table1.readAttribute(key) != table2.readAttribute(key)]
	if diff_attributes:
		table_diffs['Unequal Attributes'] = sorted(diff_attributes)

	if table_diffs:
		diffs['Table: ' + table_name] = table_diffs

def diff_simout(simout_dir1, simout_dir2):
	'''Diff two simOut dirs. Return a dict describing the differences.'''
	subdir_names1 = os.listdir(simout_dir1)
	subdir_names2 = os.listdir(simout_dir2)

	diffs = {}
	common_tables = diff_name_lists(diffs, 'Subdir', subdir_names1, subdir_names2)

	for subdir_name in common_tables:
		if subdir_name in {'Daughter1', 'Daughter2'}:
			# TODO(jerry): Compare (instead of ignoring) the inherited_state?
			pass
		else:
			diff_tables(diffs, subdir_name, simout_dir1, simout_dir2)

	return diffs

def cmd_diff_simout(simout_dir1, simout_dir2):
	'''Command line diff simout_dir2 against reference simout_dir1, with
	messages and returning an exit status code.
	'''
	if simout_dir1 == simout_dir2:
		return "diff_simouts: Don't diff a simOut directory against itself"

	print('Using simOut {} to check {}\n'.format(simout_dir1, simout_dir2))

	if not os.path.isdir(simout_dir1):
		return 'diff_simouts: No simOut directory 1: {}'.format(simout_dir1)
	if not os.path.isdir(simout_dir2):
		return 'diff_simouts: No simOut directory 2: {}'.format(simout_dir2)

	diffs = diff_simout(simout_dir1, simout_dir2)

	if diffs:
		pprint(diffs)
		return 1
	else:
		print('The simOut dirs match')
		return 0


if __name__ == '__main__':
	if len(sys.argv) == 2 and sys.argv[1] in {'-h', '--help'}:
		print('diff_simouts <simOut1> <simOut2> -- diffs two WCM simOut directories')
		sys.exit()

	if len(sys.argv) < 3:
		sys.exit('diff_simouts: Need 2 simOut directory arguments')

	exit_status = cmd_diff_simout(sys.argv[1], sys.argv[2])
	sys.exit(exit_status)
