"""
Template for variant analysis plots

TODO:
	add labels for variants
	add labels for bar plot
	add x labels
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


def density_plot(x, y, n_bins=200):
	k = stats.kde.gaussian_kde([x, y])
	xi, yi = np.mgrid[x.min():x.max():n_bins * 1j, y.min():y.max():n_bins * 1j]
	zi = k(np.vstack([xi.flatten(), yi.flatten()]))
	plt.contour(xi, yi, zi.reshape(xi.shape))


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		n_aas = len(sim_data.molecule_groups.amino_acids)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		n_bar = 1
		n_others = 0  # TODO: add plots for ppGpp, ribosomes, RNAP etc
		n_cols = n_aas + n_others + n_bar
		scaling = 3
		plt.figure(figsize=(scaling*n_cols, scaling*n_variants))
		gs = gridspec.GridSpec(nrows=n_variants, ncols=n_cols)

		for row, variant in enumerate(variants):
			cell_paths = ap.get_cells(variant=[variant])

			# Load data
			# TODO: load AA, synthesis enzymes, synthetases and overlay on same plots?
			growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True).squeeze()
			uncharged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'uncharged_trna_conc', remove_first=True).T

			# Plot scatter plot and format axes
			r = []
			for col, conc in enumerate(uncharged_trna):
				ax = plt.subplot(gs[row, col])
				density_plot(conc, growth)

				self.remove_border(ax)
				ax.tick_params(labelsize=6)
				if col != 0:
					ax.set_yticklabels([])

				r.append(stats.pearsonr(conc, growth)[0])

			# Plot bar plot of correlations
			ax = plt.subplot(gs[row, -1])
			plt.bar(range(len(r)), r)

			self.remove_border(ax)
			ax.tick_params(labelsize=6)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
