"""
Compare cell properties at varying levels of ppGpp concentration.  Useful with
with the ppgpp_conc variant and in comparison with data presented in Zhu et al.
2019. https://academic.oup.com/nar/article/47/9/4684/5420536.
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.ppgpp_conc import BASE_FACTOR, CONDITIONS, split_index
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import COLORS_SMALL


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_data(self, axes, ppgpp, y, yerr, ylabel, condition_labels, conditions, factors):
		condition_base_factor = dict(zip(CONDITIONS, BASE_FACTOR))

		raw_ax, norm_ax = axes
		for condition, color in zip(np.unique(conditions), COLORS_SMALL):
			mask = conditions == condition
			raw_ax.errorbar(ppgpp[mask], y[mask], yerr=yerr[mask], fmt='o', color=color,
				label=condition_labels[condition])

			if condition_base_factor[condition] in factors[mask]:
				ref_idx = factors[mask] == condition_base_factor[condition]
				ref_val = y[mask][ref_idx]
				norm_ax.errorbar(ppgpp[mask], y[mask] / ref_val, yerr=yerr[mask] / ref_val,
					fmt='o', color=color, label=condition_labels[condition])
				norm_ax.axhline(1, linestyle='--', color='k', linewidth=0.5)

		raw_ax.set_ylabel(ylabel, fontsize=8)
		norm_ax.set_ylabel(f'Normalized {ylabel}', fontsize=8)
		for ax in axes:
			ax.set_xlabel('ppGpp conc', fontsize=8)
			ax.tick_params(labelsize=6)
			self.remove_border(ax)

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		conditions = np.zeros(n_variants, int)
		factors = np.zeros(n_variants)
		ppgpp_mean = np.zeros(n_variants)
		ppgpp_std = np.zeros(n_variants)
		growth_rate_mean = np.zeros(n_variants)
		growth_rate_std = np.zeros(n_variants)
		rna_to_protein_mean = np.zeros(n_variants)
		rna_to_protein_std = np.zeros(n_variants)
		elong_rate_mean = np.zeros(n_variants)
		elong_rate_std = np.zeros(n_variants)
		for i, variant in enumerate(variants):
			all_cells = ap.get_cells(variant=[variant], only_successful=True)
			conditions[i], factors[i] = split_index(variant)

			# Read data from listeners
			ppgpp = read_stacked_columns(all_cells, 'GrowthLimits', 'ppgpp_conc')
			elong_rate = read_stacked_columns(all_cells, 'RibosomeData', 'effectiveElongationRate')
			growth_rate = read_stacked_columns(all_cells, 'Mass', 'instantaneous_growth_rate', remove_first=True) * 3600
			rna_mass = read_stacked_columns(all_cells, 'Mass', 'rnaMass')
			protein_mass = read_stacked_columns(all_cells, 'Mass', 'proteinMass')
			rna_to_protein = (rna_mass / protein_mass)

			# Calculate mean and std for each value
			ppgpp_mean[i] = ppgpp.mean()
			ppgpp_std[i] = ppgpp.std()
			growth_rate_mean[i] = growth_rate.mean()
			growth_rate_std[i] = growth_rate.std()
			rna_to_protein_mean[i] = rna_to_protein.mean()
			rna_to_protein_std[i] = rna_to_protein.std()
			elong_rate_mean[i] = elong_rate.mean()
			elong_rate_std[i] = elong_rate.std()

		condition_labels = sim_data.ordered_conditions

		# Create plots
		_, axes = plt.subplots(3, 2, figsize=(10, 10))

		## Bar plots of cell properties
		self.plot_data(axes[0, :], ppgpp_mean, growth_rate_mean, growth_rate_std,
			'Growth rate (1/hr)', condition_labels, conditions, factors)
		self.plot_data(axes[1, :], ppgpp_mean, elong_rate_mean, elong_rate_std,
			'Elongation rate (AA/s)', condition_labels, conditions, factors)
		self.plot_data(axes[2, :], ppgpp_mean, rna_to_protein_mean, rna_to_protein_std,
			'RNA/protein', condition_labels, conditions, factors)

		plt.legend(fontsize=6)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
