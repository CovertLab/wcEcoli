"""
Generates a .json file containing the dynamic locations of molecules bound to
the chromosome.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/2/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os
import json

import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read replichore lengths from sim_data
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		replication_data_reader = TableReader(
			os.path.join(simOutDir, "ReplicationData"))

		# Load time data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		# Load replisome attributes
		fork_coordinates = replication_data_reader.readColumn("fork_coordinates")
		domain_indexes = replication_data_reader.readColumn("fork_domains")
		fork_unique_index = replication_data_reader.readColumn("fork_unique_index")

		# Get unique indexes of each fork
		all_unique_indexes = np.unique(fork_unique_index[~np.isnan(fork_unique_index)])
		n_replisomes = len(all_unique_indexes)

		# Parse data such that one column corresponds to one unique replisome
		fork_coordinates_parsed = np.zeros((fork_coordinates.shape[0], n_replisomes))
		domain_indexes_parsed = np.zeros(n_replisomes, dtype=np.int64)

		for mol_idx, unique_idx in enumerate(all_unique_indexes):
			rows, cols = np.where(fork_unique_index == unique_idx)
			fork_coordinates_parsed[rows, mol_idx] = fork_coordinates[rows, cols]
			domain_indexes_parsed[mol_idx] = domain_indexes[rows[0], cols[0]]

		# Build dictionary of chromosome data
		chromosome_data = {
			"metadata": metadata,
			"time": time.tolist(),
			"right_replichore_length": replichore_lengths[0],
			"left_replichore_length": replichore_lengths[1],
			"replisomes": {
					"counts": n_replisomes,
					"coordinates": fork_coordinates_parsed.tolist(),
					"domain_indexes": domain_indexes_parsed.tolist()
				},
			}

		# Output dictionary to json file
		with open(os.path.join(plotOutDir, plotOutFileName + ".json"), 'w') as f:
			f.write(json.dumps(chromosome_data, indent=4))


if __name__ == '__main__':
	Plot().cli()
