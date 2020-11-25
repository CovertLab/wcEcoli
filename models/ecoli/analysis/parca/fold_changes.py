"""
Compare fold changes from different sources.
"""

from six.moves import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = cPickle.load(f)

		# Load original fold change data
		original_fold_changes = {}
		consistent = {}
		for fc in raw_data.fold_changes:
			pair = (fc['TF'], fc['Target'])
			original_fold_changes[pair] = fc['log2 FC mean']
			consistent[pair] = fc['Regulation_direct'] <= 2

		# Load NCA fold change data
		nca_fold_changes = {}
		for fc in raw_data.fold_changes_nca:
			nca_fold_changes[(fc['TF'], fc['Target'])] = fc['log2 FC mean']

		# Compare regulation pairs in both datasets
		original = []
		nca = []
		for pair in original_fold_changes:
			if pair in nca_fold_changes and consistent[pair]:
				original.append(original_fold_changes[pair])
				nca.append(nca_fold_changes[pair])
		original = np.array(original)
		nca = np.array(nca)
		pearson = stats.pearsonr(original, nca)

		# Get range for y=x line
		min_val = min(original.min(), nca.min())
		max_val = max(original.max(), nca.max())

		plt.figure()

		# Plot scatter
		plt.plot(original, nca, 'x', label='log2 fold changes')
		xlim = plt.xlim()
		ylim = plt.ylim()

		# Plot y=x and revert to previous limits
		plt.plot([min_val, max_val], [min_val, max_val], 'k--', label='y=x')
		plt.xlim(xlim)
		plt.ylim(ylim)

		# Format axes
		plt.xlabel('Original fold change')
		plt.ylabel('NCA predicted fold change')
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.legend(fontsize=8, frameon=False)

		# Display stats
		plt.title(f'Pearson r={pearson[0]:.3f} (p={pearson[1]:.0e}, n={len(original)})', fontsize=10)

		# Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
