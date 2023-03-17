"""
Plot to identify candidate sub generational genes
"""
from __future__ import absolute_import, division, print_function

import pickle
import os

from wholecell.utils.sparkline import whitePadSparklineAxis
import numpy as np
from matplotlib import pyplot as plt
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		simOutDir = os.path.join(cell_paths[0], "simOut")

		time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()

		cell_cycle_length = read_stacked_columns(cell_paths, 'Main', 'time',fun=lambda x: (x[-1] - x[0])).squeeze()

		generation_time = np.cumsum(cell_cycle_length)

		generation_index = np.searchsorted(time, generation_time)

		# Determine all subgen genes
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		all_cistron_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts')

		cistron_expressed_bool = (all_cistron_counts == 0)*1

		cistron_expressed_bool_by_generation = np.split(cistron_expressed_bool, generation_index, axis=0)

		cistron_expressed_by_generation = np.array([
			np.sum(cistron_expressed_bool_by_generation[i], axis = 0)
			for i in range(len(cell_cycle_length))
		])

		frequency_of_cistron_per_generation = cistron_expressed_by_generation/cell_cycle_length[:, None]

		cistron_expressed_per_generation_bool = (cistron_expressed_by_generation >= 1) * 1

		frequency_of_cistron_over_all_generations = np.sum(cistron_expressed_per_generation_bool, axis = 0)/len(cell_cycle_length)

		subgenerational_genes_mask = (
				(frequency_of_cistron_over_all_generations > 0)
				& (frequency_of_cistron_over_all_generations < 1)
		)

		subgenerational_gene_ids = np.array(mRNA_cistron_ids)[subgenerational_genes_mask]
		not_expressed_gene_ids = np.array(mRNA_cistron_ids)[frequency_of_cistron_over_all_generations == 0]
		always_expressed_gene_ids = np.array(mRNA_cistron_ids)[frequency_of_cistron_over_all_generations == 1]

		plt.figure(figsize=(8.5, 11))
		x = ['always_expressed','subgenerational', 'not_expressed']
		y = [len(always_expressed_gene_ids), len(subgenerational_gene_ids), len(not_expressed_gene_ids)]

		plt.bar(x, y, color = 'b')

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close()


if __name__ == "__main__":
	Plot().cli()
