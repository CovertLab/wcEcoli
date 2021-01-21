"""
Compare expected promoter binding probabilities with observed probabilities.
"""

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(variantDir, cohort_plot=True)

		expected = []
		bound = []
		promoters = []
		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			rna_synth_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

			# Load data
			expected.append(rna_synth_reader.readColumn('pPromoterBound'))
			bound.append(rna_synth_reader.readColumn('nActualBound'))
			promoters.append(rna_synth_reader.readColumn('n_available_promoters'))

		expected = np.vstack(expected).mean(0)
		bound = np.vstack(bound)
		promoter_prob = (bound / np.vstack(promoters))
		promoter_prob[bound == 0] = 0
		actual = promoter_prob.mean(0)

		min_prob = min(expected[expected > 0].min(), actual[actual > 0].min())
		max_prob = np.max((expected, actual))

		plt.figure()
		ax = plt.gca()

		# Plot data
		ax.plot([0, max_prob], [0, max_prob], '--k')
		ax.plot(expected, actual, 'o', alpha=0.5)
		ax.set_xscale('symlog', linthreshx=min_prob)
		ax.set_yscale('symlog', linthreshy=min_prob)

		# Axes formatting
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.set_xlabel(f'Expected promoter bound probability')
		ax.set_ylabel(f'Observed promoter bound probability')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
