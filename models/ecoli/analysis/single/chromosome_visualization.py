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
		rnap_data_reader = TableReader(
			os.path.join(simOutDir, "RnapData"))

		# Load time data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		n_timesteps = len(time)

		# Load replisome attributes
		fork_coordinates = replication_data_reader.readColumn("fork_coordinates")
		fork_domain_indexes = replication_data_reader.readColumn("fork_domains")
		fork_unique_indexes = replication_data_reader.readColumn("fork_unique_index")

		# Load active RNAP attributes
		rnap_coordinates = rnap_data_reader.readColumn("active_rnap_coordinates")
		rnap_domain_indexes = rnap_data_reader.readColumn("active_rnap_domain_indexes")
		rnap_unique_indexes = rnap_data_reader.readColumn("active_rnap_unique_indexes")

		# Get list of unique indexes of each replisome
		fork_unique_index_list = np.unique(
			fork_unique_indexes[~np.isnan(fork_unique_indexes)])
		n_unique_replisomes = len(fork_unique_index_list)

		# Parse data such that one row of column corresponds to one unique
		# replisome. The active mask boolean array is set to True if the
		# corresponding replisome (column) is active at the given timestep
		# (row).
		fork_coordinates_parsed = np.zeros((n_timesteps, n_unique_replisomes))
		fork_active_mask_parsed = np.zeros(
			(n_timesteps, n_unique_replisomes), dtype=np.bool)
		fork_domain_indexes_parsed = np.zeros(n_unique_replisomes, dtype=np.int64)

		for mol_idx, unique_idx in enumerate(fork_unique_index_list):
			time_index, col_index = np.where(fork_unique_indexes == unique_idx)
			fork_coordinates_parsed[time_index, mol_idx] = fork_coordinates[time_index, col_index]
			fork_active_mask_parsed[time_index, mol_idx] = True

			# Domain indexes are static - just take the first value
			fork_domain_indexes_parsed[mol_idx] = fork_domain_indexes[
				time_index[0], col_index[0]]

		# Crop out full columns of NaNs and replace NaNs to zeros for RNAP data
		rnap_isnan = np.isnan(rnap_coordinates)
		n_nan_columns = (rnap_isnan.sum(axis=0) == n_timesteps).sum()

		rnap_coordinates_cropped = np.nan_to_num(rnap_coordinates[:, :-n_nan_columns])
		rnap_domain_indexes_cropped = np.nan_to_num(rnap_domain_indexes[:, :-n_nan_columns])
		rnap_unique_indexes_cropped = np.nan_to_num(rnap_unique_indexes[:, :-n_nan_columns])

		# Build dictionary of chromosome data
		chromosome_data = {
			"metadata": metadata,
			"time": [round(t, 2) for t in time],
			"right_replichore_len": replichore_lengths[0],
			"left_replichore_len": replichore_lengths[1],
			"replisomes": {
				"active": fork_active_mask_parsed.tolist(),
				"counts": n_unique_replisomes,
				"coordinates": fork_coordinates_parsed.tolist(),
				"domain_indexes": fork_domain_indexes_parsed.tolist()
				},
			"active_RNAPs": {
				"active": np.logical_not(rnap_isnan).tolist(),
				"coordinates": rnap_coordinates_cropped.tolist(),
				"domain_indexes": rnap_domain_indexes_cropped.tolist(),
				"unique_indexes": rnap_unique_indexes_cropped.tolist()
				}
			}

		# Output dictionary to json file
		with open(os.path.join(plotOutDir, plotOutFileName + ".json"), 'w') as f:
			f.write(json.dumps(chromosome_data))


if __name__ == '__main__':
	Plot().cli()
