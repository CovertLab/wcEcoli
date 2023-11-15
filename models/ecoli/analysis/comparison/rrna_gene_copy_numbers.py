"""
Comparison plot to compare the copy numbers of genes encoding rRNAs.
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


GENE_ID_TO_RRNA_OPERON_ID = {
	'EG30084': 'rrnA',
	'EG30085': 'rrnB',
	'EG30086': 'rrnC',
	'EG30087': 'rrnD',
	'EG30088': 'rrnE',
	'EG30089': 'rrnG',
    'EG30090': 'rrnH',
	}

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if ap1.n_generation <= 4 or ap2.n_generation <= 4:
			print('Not enough generations to run analysis.')
			return

		def read_sims(ap):
			# Get cell paths
			cell_paths = ap.get_cells(
				generation=np.arange(4, ap.n_generation),
				only_successful=True)

			# Get gene_ids attribute from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			gene_ids = rna_synth_prob_reader.readAttribute('gene_ids')

			# Get indexes of 16S genes (first gene in each operon)
			rrna_gene_indexes = np.array([
				gene_ids.index(key) for key in GENE_ID_TO_RRNA_OPERON_ID.keys()
				])

			# Get copy numbers of 16S genes
			rrna_gene_copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'gene_copy_number',
				ignore_exception=True, fun=lambda x: x[:, rrna_gene_indexes])
			avg_rrna_gene_copy_numbers = rrna_gene_copy_numbers.mean(axis=0)

			return avg_rrna_gene_copy_numbers

		copy_numbers1 = read_sims(ap1)
		copy_numbers2 = read_sims(ap2)

		fig = plt.figure(figsize=(9, 3))
		gs = fig.add_gridspec(1, 2, width_ratios=(6.5, 1))

		# Plot bar plots for each rRNA operon
		ax0 = fig.add_subplot(gs[0, 0])
		for i, (c1, c2) in enumerate(zip(copy_numbers1, copy_numbers2)):
			ax0.bar(
				1.5*i - 0.25, c1, width=0.5, alpha=0.5, color='C0', label='reference')
			ax0.bar(
				1.5*i + 0.25, c2, width=0.5, alpha=0.5, color='C1', label='input')

		ax0.set_xticks(1.5*np.arange(7))
		ax0.set_xticklabels([v for v in GENE_ID_TO_RRNA_OPERON_ID.values()])
		ax0.set_xlim([-0.8, 9.8])
		ax0.set_ylabel('copy numbers')
		ax0.spines['top'].set_visible(False)
		ax0.spines['right'].set_visible(False)

		# Plot bar plots for total copy numbers
		ax1 = fig.add_subplot(gs[0, 1])

		ax1.bar(
			-0.25, copy_numbers1.sum(), width=0.5, alpha=0.5, color='C0', label='reference')
		ax1.bar(
			0.25, copy_numbers2.sum(), width=0.5, alpha=0.5, color='C1', label='input')

		ax1.set_xticks([0])
		ax1.set_xticklabels([
			'all operons'])
		ax1.set_xlim([-0.8, 0.8])
		ax1.spines['top'].set_visible(False)
		ax1.spines['right'].set_visible(False)

		ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})

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
