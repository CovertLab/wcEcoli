"""
Compare histograms of mRNA lengths.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError
from wholecell.utils import units


FIGSIZE = (4, 4)
BOUNDS = [1, 5]
N_BINS = 20

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if ap1.n_generation <= 2 or ap2.n_generation <= 2:
			print('Skipping analysis -- not enough sims run.')
			return

		def read_sims(ap, sim_data):
			rna_lengths = sim_data.process.transcription.rna_data['length'].asNumber(units.nt)
			is_mRNA = sim_data.process.transcription.rna_data['is_mRNA']

			# Ignore data from first two gens
			cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

			mRNA_counts = read_stacked_columns(
				cell_paths, 'RNACounts', 'full_mRNA_counts', ignore_exception=True)

			mean_mRNA_counts = mRNA_counts.mean(axis=0)
			mRNA_lengths = rna_lengths[is_mRNA]

			return mRNA_lengths, mean_mRNA_counts

		lengths1, mean_counts1 = read_sims(ap1, sim_data1)
		lengths2, mean_counts2 = read_sims(ap2, sim_data2)

		weighted_mean1 = lengths1.dot(mean_counts1) / mean_counts1.sum()
		weighted_mean2 = lengths2.dot(mean_counts2) / mean_counts2.sum()

		m1 = (mean_counts1 > 0).sum()
		m2 = (mean_counts2 > 0).sum()

		weighted_std1 = np.sqrt(
			mean_counts1.dot((lengths1 - weighted_mean1)**2) /
			((m1 - 1) / m1 * mean_counts1.sum())
			)
		weighted_std2 = np.sqrt(
			mean_counts2.dot((lengths2 - weighted_mean2)**2) /
			((m2 - 1) / m2 * mean_counts2.sum())
			)

		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)

		bins = np.logspace(BOUNDS[0], BOUNDS[1], N_BINS + 1)

		ax.hist(
			lengths1, bins=bins, weights=mean_counts1, alpha=0.5,
			label=f'reference ({weighted_mean1:.1f} $\pm$ {weighted_std1:.1f}, n={len(lengths1)})')
		ax.axvline(weighted_mean1, ls='--', lw=2, c='C0')
		ax.hist(
			lengths2, bins=bins, weights=mean_counts2, alpha=0.5,
			label=f'input ({weighted_mean2:.1f}  $\pm$ {weighted_std2:.1f}, n={len(lengths2)})')
		ax.axvline(weighted_mean2, ls='--', lw=2, c='C1')
		ax.legend(prop={'size': 6})

		ax.set_xlim([10**1.5, 10**4.5])
		ax.set_xscale('log')
		ax.set_xlabel('mRNA lengths (nt)')
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)

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
