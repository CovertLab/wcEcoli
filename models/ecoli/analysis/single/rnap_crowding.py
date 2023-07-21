"""
Comparison of target synthesis probabilities vs actual synthesis
probabilities for transcription units whose synthesis probabilities exceeded
the limit set by the physical size and the elongation rates of RNA polymerases.
"""

import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rna_synth_prob_reader =	TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time_min = (main_reader.readColumn('time') - initial_time) / 60
		tu_ids = rna_synth_prob_reader.readAttribute('rnaIds')
		target_rna_synth_prob = rna_synth_prob_reader.readColumn(
			'target_rna_synth_prob')
		actual_rna_synth_prob = rna_synth_prob_reader.readColumn(
			'actual_rna_synth_prob')
		tu_is_overcrowded = rna_synth_prob_reader.readColumn(
			'tu_is_overcrowded')

		# Get indexes of TUs that were overcrowded at some point in the sim
		overcrowded_tu_indexes = np.where(tu_is_overcrowded.sum(axis=0))[0]
		n_overcrowded_tus = len(overcrowded_tu_indexes)

		if n_overcrowded_tus > 0:
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
		else:
			# Generate empty plot if no overcrowding occurred
			plt.figure(figsize=(6, 1.5))
			ax = plt.subplot(1, 1, 1)
			ax.set_title('No TUs were overcrowded.')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
