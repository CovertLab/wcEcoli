"""
Compare histograms of active ribosome counts between two sets of simulations.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError


FIGSIZE = (4, 4)
BOUNDS = [0, 30000]
N_BINS = 30

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if ap1.n_generation <= 4 or ap2.n_generation <= 4:
			print('Not enough generations to run analysis.')
			return

		def read_sims(ap):
			# Ignore data from first four gens
			cell_paths = ap.get_cells(generation=np.arange(4, ap.n_generation))

			# Get index of active ribosomes in the unique molecule counts reader
			unique_molecule_counts_reader = TableReader(
				os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
			unique_molecule_ids = unique_molecule_counts_reader.readAttribute(
				'uniqueMoleculeIds')
			active_ribosome_idx = unique_molecule_ids.index('active_ribosome')

			# Get counts of active ribosome at first timestep from every cell
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts','uniqueMoleculeCounts',
				ignore_exception=True, fun=lambda x: x[0, active_ribosome_idx])

			return active_ribosome_counts

		# Get initial active ribosome counts from each set of sims
		active_ribosome_counts1 = read_sims(ap1)
		active_ribosome_counts2 = read_sims(ap2)

		# Plot histogram
		fig = plt.figure(figsize=FIGSIZE)
		ax = fig.add_subplot(1, 1, 1)

		bins = np.linspace(BOUNDS[0], BOUNDS[1], N_BINS + 1)

		ax.hist(
			active_ribosome_counts1, bins=bins, alpha=0.5,
			label=f'reference ({np.mean(active_ribosome_counts1):.2f} $\pm$ {np.std(active_ribosome_counts1):.2f}, n={len(active_ribosome_counts1)})')
		ax.axvline(np.mean(active_ribosome_counts1), ls='--', lw=2, c='C0')
		ax.hist(
			active_ribosome_counts2, bins=bins, alpha=0.5,
			label=f'input ({np.mean(active_ribosome_counts2):.2f} $\pm$ {np.std(active_ribosome_counts2):.2f}, n={len(active_ribosome_counts2)})')
		ax.axvline(np.mean(active_ribosome_counts2), ls='--', lw=2, c='C1')

		ax.legend(prop={'size': 6})

		ax.set_xlim(BOUNDS)
		ax.set_xlabel('Active ribosome counts')
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
