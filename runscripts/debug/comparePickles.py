#! /usr/bin/env python
"""
Compare all the .cPickle files in a pair of directories. Show any differences.

Usage (DIR is a path like 'out/manual/intermediates'):
	runscripts/debug/comparePickles.py DIR1 DIR2
"""

import argparse
import os
import sys

from runscripts.reflect.object_tree import diff_dirs, diff_files


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Compare two .cPickle files"
					" or all the .cPickle files in two directories (in"
					" modification-time order)."
					" Print a count and optionally a summary of the differences.")
	parser.add_argument('-c', '--count', action='store_true',
		help="Print just the diff line count for each file, skipping the"
			 " detailed diff lines.")
	parser.add_argument('path', metavar='PATH', nargs=2,
		help="The two pickle files or directories to compare.")

	args = parser.parse_args()
	path1, path2 = args.path

	if os.path.isfile(path1):
		diff_count = diff_files(path1, path2, print_diff_lines=not args.count)
	else:
		diff_count = diff_dirs(path1, path2, print_diff_lines=not args.count)

	sys.exit(3 if diff_count else 0)
