"""
Test cases for new tablereader methods

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import
from __future__ import division

import os
import time

import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

ITERS = 10
BLOCK_SIZE = 5000  # roughly number of proteins or RNA


def test_old(reader, indices):
	start = time.time()
	for i in range(ITERS):
		counts = reader.readColumn('counts')[:, indices]
	end = time.time()

	print('\tOld method: {:.3f}'.format(end - start))

	return counts

def test_new_block(reader, indices):
	start = time.time()
	for i in range(ITERS):
		counts = reader.readColumn('counts', indices)
	end = time.time()

	print('\tNew method, block read: {:.3f}'.format(end - start))

	return counts

def test_new_multiple(reader, indices):
	start = time.time()
	for i in range(ITERS):
		counts = reader.readColumn('counts', indices, False)
	end = time.time()

	print('\tNew method, multiple reads: {:.3f}'.format(end - start))

	return counts

def test_all_functions(text, reader, indices):
	functions = [test_old, test_new_block, test_new_multiple]

	print(text)
	results = [f(reader, indices) for f in functions]
	if np.all(results[0] == results[1]) and np.all(results[0] == results[2]):
		print('\tResults match!')
	else:
		print('\tResults do not match')

def test_two_functions(text, reader, indices):
	functions = [test_old, test_new_block]

	print(text)
	results = [f(reader, indices) for f in functions]
	if np.all(results[0] == results[1]):
		print('\tResults match!')
	else:
		print('\tResults do not match')

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		# Listeners used
		bulk_molecules = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
		bulk_ids = bulk_molecules.readAttribute('objectNames')
		n_mols = len(bulk_ids)

		# Test reads
		## Single index
		indices = np.array([0])
		test_all_functions('One index', bulk_molecules, indices)

		## First and last index
		indices = np.array([0, n_mols-1])
		test_all_functions('First and last indices', bulk_molecules, indices)

		## Large block
		indices = np.array(range(BLOCK_SIZE))
		test_all_functions('Block indices', bulk_molecules, indices)

		## 2 Large blocks
		indices = np.array(range(BLOCK_SIZE) + range(n_mols)[-BLOCK_SIZE:])
		test_all_functions('Two blocks of indices', bulk_molecules, indices)

		## Dispersed reads
		indices = np.linspace(0, n_mols-1, BLOCK_SIZE, dtype=np.int64)
		test_all_functions('Dispersed indices', bulk_molecules, indices)

		## Random reads
		indices = np.array(range(n_mols))
		np.random.shuffle(indices)
		indices = indices[:BLOCK_SIZE]
		test_all_functions('Random indices', bulk_molecules, indices)

		## All indices
		indices = np.array(range(n_mols))
		test_all_functions('All indices', bulk_molecules, indices)


if __name__ == '__main__':
	Plot().cli()
