"""
Plots a violin/swarm plot of the relative amounts of the sum of excess protein
monomers produced for heterogeneous protein complexes.
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		complex_ids = sim_data.process.complexation.ids_complexes
		subunit_ids = sim_data.process.complexation.molecule_names
		stoich_matrix_monomers = sim_data.process.complexation.stoich_matrix_monomers()

		# Listeners used
		monomer_counts_reader =	TableReader(os.path.join(simOutDir, 'MonomerCounts'))

		# Load data
		(all_complex_counts, ) = read_bulk_molecule_counts(simOutDir, (complex_ids, ))
		all_monomer_counts = monomer_counts_reader.readColumn('monomerCounts')
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)}

		# Take averages across timepoints
		all_complex_counts_mean = all_complex_counts.mean(axis=0)
		all_monomer_counts_mean = all_monomer_counts.mean(axis=0)

		excess_monomer_fractions = []

		for i, complex_id in enumerate(complex_ids):
			subunit_mask = stoich_matrix_monomers[:, i] < 0
			monomer_stoichs = -stoich_matrix_monomers[:, i][subunit_mask]
			monomer_indexes = np.where(subunit_mask)[0]

			# Skip homogeneous protein complexes
			if len(monomer_indexes) == 1:
				continue

			try:
				monomer_indexes = np.array([
					monomer_id_to_index[subunit_ids[subunit_index]]
					for subunit_index in monomer_indexes
					])
			except KeyError:
				# Skip complexes with non-protein monomers or monomers that
				# were excluded from the model
				continue

			monomer_counts = all_monomer_counts_mean[monomer_indexes]
			complex_count = all_complex_counts_mean[i]

			if np.all(monomer_counts > 0):
				excess_monomer_fractions.append(
					((monomer_counts - monomer_stoichs*complex_count) / monomer_counts).mean()
					)

		plt.figure()

		# TODO: Build violin/swarm plot, highlight complexes with operon-bound
		# 	monomers?

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
