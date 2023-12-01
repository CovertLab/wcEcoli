"""
Comparison plot to compare the average number of collisions that occur between
replisomes and RNAPs.
"""

import os
from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure,
    read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader



class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		def read_sims(ap, sim_data):
			# Get total counts of collisions from all successful sims
			cell_paths = ap.get_cells(only_successful=True)
			n_cells = len(cell_paths)
			all_headon_collisions = read_stacked_columns(
				cell_paths, 'RnapData', 'n_headon_collisions',
				ignore_exception=True, fun=lambda x: x.sum())
			all_codirectional_collisions = read_stacked_columns(
				cell_paths, 'RnapData', 'n_codirectional_collisions',
				ignore_exception=True, fun=lambda x: x.sum())
			total_headon_collisions = all_headon_collisions.sum()
			total_codirectional_collisions = all_codirectional_collisions.sum()

			# Divide by number of cells
			return {
				'n_headon_collisions': total_headon_collisions / n_cells,
				'n_codirectional_collisions': total_codirectional_collisions / n_cells,
			}

		data1 = read_sims(ap1, sim_data1)
		data2 = read_sims(ap2, sim_data2)

		fig = plt.figure(figsize=(4, 5))

		# Plot bar plots for number of collisions
		ax = fig.add_subplot(1, 1, 1)
		ax.bar(0, data1['n_codirectional_collisions'], width=0.7, alpha=0.5,
			   color='#555555', label='codirectional')
		ax.bar(0, data1['n_headon_collisions'], width=0.7, alpha=0.5,
			   color='C3', bottom=data1['n_codirectional_collisions'], label='headon')
		ax.bar(1, data2['n_codirectional_collisions'], width=0.7, alpha=0.5,
			   color='#555555')
		ax.bar(1, data2['n_headon_collisions'], width=0.7, alpha=0.5,
			   color='C3', bottom=data2['n_codirectional_collisions'])

		ax.set_xticks([0, 1])
		ax.set_xticklabels(['reference', 'input'])
		ax.set_xlim([-0.8, 1.8])
		ax.set_ylabel('Average number of replisome-RNAP collisions')
		ax.legend(loc='upper left', prop={'size': 8}, bbox_to_anchor=(1.04, 1))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
