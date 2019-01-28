"""
Template for single analysis plots

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
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

		# Build dictionary of chromosome data
		chromosome_data = {
			"metadata": metadata,
			"time": time.tolist(),
			"right_replichore_length": replichore_lengths[0],
			"left_replichore_length": replichore_lengths[1],
			"replisomes": {
					"coordinates": fork_coordinates.tolist(),
					"domain_indexes": domain_indexes.tolist()
				}
			}

		# Output dictionary to json file
		out_file = open(os.path.join(plotOutDir, plotOutFileName + ".json"), 'w')
		out_file.write(json.dumps(chromosome_data, indent=4))


if __name__ == '__main__':
	Plot().cli()
