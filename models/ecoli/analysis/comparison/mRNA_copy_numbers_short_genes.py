"""
Generates a scatter plot comparing mRNA copy numbers for short genes whose
RNA-Seq read counts are misreported as zero.
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
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


FIGSIZE = (4, 3.9)
BOUNDS = [0.5, 2.5]

SHORT_GENE_IDS_TO_ACTUAL_READ_COUNTS = {
	'EG11269': 13.016887039317732,
	'EG11271': 4.579771048659526,
	'EG11272': 0.3918528586244949,
	'EG11274': 10.090414561057367,
	'EG11276': 0.0,
	'EG11277': 31.420200804608612,
	'EG11280': 22.84349367645335,
	'G0-10439': 0.0,
	'G0-10580': 8.51021823855618,
	'G0-10588': 0.19765505228974092,
	'G0-10641': 25.171596929773155,
	'G0-10645': 0.0,
	'G0-10651': 1.934269898213109,
	'G0-10689': 5.781810967296733,
	}


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		bounds = BOUNDS

		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, validation_data2 = self.setup(input_sim_dir)

		cell_paths = ap2.get_cells()
		simOutDir = os.path.join(cell_paths[0], 'simOut')
		RNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_ids = RNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Get mask for mRNA genes that encode for ribosomal proteins or RNAPs
		all_cistron_ids = sim_data2.process.transcription.cistron_data['id']
		is_mRNA = sim_data2.process.transcription.cistron_data['is_mRNA']
		assert np.all(
			mRNA_cistron_ids == all_cistron_ids[is_mRNA])

		mRNA_gene_ids = sim_data2.process.transcription.cistron_data['gene_id'][is_mRNA]
		gene_id_to_index = {
			gene_id: i for i, gene_id in enumerate(mRNA_gene_ids)
			}
		short_gene_indexes = np.array([
			gene_id_to_index[gene_id] for gene_id in SHORT_GENE_IDS_TO_ACTUAL_READ_COUNTS])
		actual_read_counts = np.array([
			v for v in SHORT_GENE_IDS_TO_ACTUAL_READ_COUNTS.values()])

		def read_data(ap):
			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			# Sample initial mRNA counts from each cell
			all_initial_counts = read_stacked_columns(
				cell_paths, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True, fun=lambda x: x[0])

			return all_initial_counts

		c1 = read_data(ap1)
		c2 = read_data(ap2)

		if len(c1) == 0 or len(c2) == 0:
			print('Skipping analysis -- not enough sims run.')
			return

		m1 = c1.mean(axis=0)
		m2 = c2.mean(axis=0)

		# Normalize counts from two conditions
		ratio = m1[all_plotted_mRNA_mask].sum()/m2[all_plotted_mRNA_mask].sum()
		m2 = ratio * m2.astype(float)

		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)
		ax.plot(bounds, bounds, ls='--', lw=2, c='k', alpha=0.05)
		ax.scatter(
			np.log10(actual_read_counts + 1),
			np.log10(m1[short_gene_indexes] + 1),
			c='#aaaaaa', edgecolor='none', s=12, alpha=0.7,
			label=f'ribosomal subunits',
			clip_on=False)
		ax.scatter(
			np.log10(actual_read_counts + 1),
			np.log10(m2[short_gene_indexes] + 1),
			c='#222222', edgecolor='none', s=12, alpha=0.7,
			label='RNAP subunits',
			clip_on=False)

		ax.set_xlabel('$\log_{10}$(normalized mRNA copies + 1), from manual alignment')
		ax.set_ylabel('$\log_{10}$(normalized mRNA copies + 1), from simulated cells')
		ax.set_xticks(np.arange(bounds[0], bounds[1] + 0.5, 0.5))
		ax.set_yticks(np.arange(bounds[0], bounds[1] + 0.5, 0.5))
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_position(("outward", 15))
		ax.spines["left"].set_position(("outward", 15))
		ax.set_xlim(bounds)
		ax.set_ylim(bounds)
		ax.legend(loc=2, prop={'size': 8})

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


if __name__ == '__main__':
	Plot().cli()
