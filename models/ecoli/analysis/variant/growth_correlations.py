"""
Template for variant analysis plots

TODO:
	add labels for variants
	add labels for bar plot
	add x labels
	2D density plot instead of points
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


MA_WINDOW = 151  # 5 min moving average if timestep is 2 sec


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_bar(self, ax, cell_property, growth_rates, xlabel, tick_labels):
		r = np.zeros((cell_property.shape[0], len(growth_rates)))
		for i, property in enumerate(cell_property):
			for j, growth in enumerate(growth_rates):
				r[i, j] = stats.pearsonr(property[:len(growth)], growth)[0]

		ax.bar(range(r.shape[0]), r[:, 0])
		self.remove_border(ax)

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		n_bar = 1
		n_others = 0  # TODO: add plots for ppGpp, ribosomes, RNAP etc
		n_cols = n_aas + n_others + n_bar
		scaling = 3
		n_cols = 7
		plt.figure(figsize=(scaling*n_cols, scaling*n_variants))
		gs = gridspec.GridSpec(nrows=n_variants, ncols=n_cols)

		for row, variant in enumerate(variants):
			cell_paths = ap.get_cells(variant=[variant])

			# Load data
			# TODO: load AA, synthesis enzymes, synthetases and overlay on same plots?
			counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True).squeeze()
			growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
				remove_first=True).squeeze() * 3600
			uncharged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'uncharged_trna_conc',
				remove_first=True).T
			charged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'charged_trna_conc',
				remove_first=True).T
			aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc',
				remove_first=True).T
			synthetase_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'synthetase_conc',
				remove_first=True).T
			aa_supply_enzymes_fwd = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_fwd',
				remove_first=True).T * counts_to_molar
			aa_supply_enzymes_rev = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_rev',
				remove_first=True).T * counts_to_molar

			growth_ma = np.convolve(growth, np.ones(MA_WINDOW) / MA_WINDOW, mode='valid')
			fraction_charged = charged_trna / (uncharged_trna + charged_trna)

			all_growth = [growth, growth_ma]
			self.plot_bar(plt.subplot(gs[row, 0]), uncharged_trna, all_growth, '', aa_ids)
			self.plot_bar(plt.subplot(gs[row, 1]), charged_trna, all_growth, '', aa_ids)
			self.plot_bar(plt.subplot(gs[row, 2]), fraction_charged, all_growth, '', aa_ids)
			self.plot_bar(plt.subplot(gs[row, 3]), aa_conc, all_growth, '', aa_ids)
			self.plot_bar(plt.subplot(gs[row, 4]), synthetase_conc, all_growth, '', aa_ids)
			self.plot_bar(plt.subplot(gs[row, 5]), aa_supply_enzymes_fwd, all_growth, '', aa_ids)
			self.plot_bar(plt.subplot(gs[row, 6]), aa_supply_enzymes_rev, all_growth, '', aa_ids)
			# # Plot scatter plot and format axes
			# r = []
			# for col, conc in enumerate(aa_supply_enzymes_rev):
			# 	ax = plt.subplot(gs[row, col])
			# 	ax.plot(conc, growth, 'o', alpha=0.2, markersize=2)
			#
			# 	self.remove_border(ax)
			# 	ax.tick_params(labelsize=6)
			# 	if col != 0:
			# 		ax.set_yticklabels([])
			#
			# 	# TODO: r for multiple window sizes - should see effect decrease as window increases
			# 	r.append(stats.pearsonr(conc[:len(growth_ma)], growth_ma)[0])
			#
			# # Plot bar plot of correlations
			# ax = plt.subplot(gs[row, -1])
			# plt.bar(range(len(r)), r)
			#
			# self.remove_border(ax)
			# ax.tick_params(labelsize=6)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata, extension='.png')
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
