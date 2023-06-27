"""
Comparison of target synthesis probabilities vs actual synthesis
probabilities for transcription units whose synthesis probabilities exceeded
the limit set by the physical size and the elongation rates of RNA polymerases.
"""

import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Listeners used
		rna_synth_prob_reader =	TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

		# Load data
		time_min = read_stacked_columns(cell_paths, 'Main', 'time') / 60
		tu_ids = rna_synth_prob_reader.readAttribute('rnaIds')
		target_rna_synth_prob = read_stacked_columns(cell_paths,
			'RnaSynthProb', 'target_rna_synth_prob')
		actual_rna_synth_prob = read_stacked_columns(cell_paths,
			'RnaSynthProb', 'actual_rna_synth_prob')
		tu_is_overcrowded = read_stacked_columns(cell_paths,
			'RnaSynthProb', 'tu_is_overcrowded')

		# Get indexes of TUs that were overcrowded at some point in the sim
		overcrowded_tu_indexes = np.where(tu_is_overcrowded.sum(axis=0))[0]
		n_overcrowded_tus = len(overcrowded_tu_indexes)

		# Plot the target vs actual rna synthesis probabilites of these TUs
		plt.figure(figsize=(6, 1.5*n_overcrowded_tus))

		for i, tu_index in enumerate(overcrowded_tu_indexes):
			ax = plt.subplot(n_overcrowded_tus, 1, i + 1)
			ax.plot(time_min, target_rna_synth_prob[:, tu_index], label='target')
			ax.plot(time_min, actual_rna_synth_prob[:, tu_index], label='actual')
			ax.set_ylabel(f'{tu_ids[i][:-3]}\nsynthesis probs')

			if i == 0:
				ax.set_title(f'Total number of overcrowded TUs: {n_overcrowded_tus}')
				ax.legend(loc=1)

			if i == n_overcrowded_tus - 1:
				ax.set_xlabel('Time [min]')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == '__main__':
	Plot().cli()
