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

		x = np.arange(r.shape[0])
		ax.bar(x, r[:, 0])

		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xticks(x)
		ax.set_xticklabels(tick_labels, fontsize=6, rotation=45, ha='right')
		ax.set_xlabel(xlabel, fontsize=8)
		ax.set_ylabel('Growth correlation', fontsize=8)

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids]
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
			aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc',
				remove_first=True).T
			uncharged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'uncharged_trna_conc',
				remove_first=True).T
			charged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'charged_trna_conc',
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

			plot_data = [
				(aa_conc, 'Amino acid conc'),
				(uncharged_trna, 'Uncharged tRNA conc'),
				(charged_trna, 'Charged tRNA conc'),
				(fraction_charged, 'Fraction charged'),
				(synthetase_conc, 'Synthetase conc'),
				(aa_supply_enzymes_fwd, 'Forward enzyme conc'),
				(aa_supply_enzymes_rev, 'Reverse enzyme conc'),
				]
			for col, (property, label) in enumerate(plot_data):
				self.plot_bar(plt.subplot(gs[row, col]), property, all_growth, label, aa_ids)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
