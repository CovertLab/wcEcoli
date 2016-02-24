#!/usr/bin/env python

"""
Analysis firetask tools
"""

import os

def run_function(f, args, name):
	try:
		print "%s: Running %s" % (time.ctime(), name)
		f(*args)
	except KeyboardInterrupt:
		import sys; sys.exit(1)
		

def run_specific_order(directory, fileList, fileName, position=-1, verbose=False):
	"""
	Moves files mentioned in fileName to end of fileList.
	"""
	# Need to access last index directly, instead of with -1
	if position == -1:
		position = len(fileList)

	with open(os.path.join(directory, fileName), 'r') as f:
		runlast = []
		for line in f:
			line = line.strip("\n")
			if line in fileList:
				runlast.append(line)
				idx = fileList.index(line)
				fileList.insert(position, fileList.pop(idx))
			else:
				if verbose:
					print "\nUnknown file in %s: %s" % (fileName, line)
	if verbose:
		print "Running files in %s after others: \n%s\n" % (fileName, ", ".join(runlast))
	return fileList