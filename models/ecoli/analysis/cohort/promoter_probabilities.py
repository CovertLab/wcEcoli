"""
Compare expected TF binding probabilities with observed probabilities.
"""

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


def plot(ax, expected, actual, label):
	min_prob = min(expected[expected > 0].min(), actual[actual > 0].min())
	max_prob = np.max((expected, actual))

	# Plot data
	ax.loglog([min_prob, max_prob], [min_prob, max_prob], '--k')
	ax.loglog(expected, actual, 'o', alpha=0.5)

	# Axes formatting
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel(f'Expected probability')
	ax.set_ylabel(f'Observed {label} probability')


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(variantDir, cohort_plot=True)

		expected = []
		bound = []
		active_tfs = []
		total_tfs = []
		promoters = []
		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			rna_synth_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

			# Load data
			expected.append(rna_synth_reader.readColumn('pPromoterBound'))
			bound.append(rna_synth_reader.readColumn('nActualBound'))
			active_tfs.append(rna_synth_reader.readColumn('n_active_tfs'))
			total_tfs.append(rna_synth_reader.readColumn('total_tf_counts'))
			promoters.append(rna_synth_reader.readColumn('n_available_promoters'))
		expected = np.vstack(expected).mean(0)
		bound = np.vstack(bound)
		active_tf_prob = (bound / np.vstack(active_tfs))
		active_tf_prob[bound == 0] = 0
		active_tf_prob = active_tf_prob.mean(0)
		total_tf_prob = (bound / np.vstack(total_tfs))
		total_tf_prob[bound == 0] = 0
		total_tf_prob = total_tf_prob.mean(0)
		promoter_prob = (bound / np.vstack(promoters))
		promoter_prob[bound == 0] = 0
		promoter_prob = promoter_prob.mean(0)

		plt.figure(figsize=(5, 12))

		plot(plt.subplot(3, 1, 1), expected, active_tf_prob, 'active TF')
		plot(plt.subplot(3, 1, 2), expected, total_tf_prob, 'total TF')
		plot(plt.subplot(3, 1, 3), expected, promoter_prob, 'promoter')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
