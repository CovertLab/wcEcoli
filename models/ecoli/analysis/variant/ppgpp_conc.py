"""
Compare cell properties at varying levels of ppGpp concentration.  Useful with
with the ppgpp_conc variant and in comparison with data presented in Zhu et al.
2019. https://academic.oup.com/nar/article/47/9/4684/5420536.
"""

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_data(self, ax, ppgpp, y, yerr, ylabel):
		ax.errorbar(ppgpp, y, yerr=yerr, fmt='o')

		plt.xlabel('ppGpp conc', fontsize=8)
		plt.ylabel(ylabel, fontsize=8)

		ax.tick_params(labelsize=6)
		self.remove_border(ax)

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		ppgpp_mean = np.zeros(n_variants)
		ppgpp_std = np.zeros(n_variants)
		growth_rate_mean = np.zeros(n_variants)
		growth_rate_std = np.zeros(n_variants)
		rna_to_protein_mean = np.zeros(n_variants)
		rna_to_protein_std = np.zeros(n_variants)
		elong_rate_mean = np.zeros(n_variants)
		elong_rate_std = np.zeros(n_variants)
		for i, variant in enumerate(variants):
			all_cells = ap.get_cells(variant=[variant])

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

		# Create plots
		plt.figure(figsize=(5, 10))
		gs = gridspec.GridSpec(3, 1)

		## Bar plots of cell properties
		self.plot_data(plt.subplot(gs[0, 0]), ppgpp_mean, growth_rate_mean, growth_rate_std, 'Growth rate (1/hr)')
		self.plot_data(plt.subplot(gs[1, 0]), ppgpp_mean, elong_rate_mean, elong_rate_std, 'Elongation rate (AA/s)')
		self.plot_data(plt.subplot(gs[2, 0]), ppgpp_mean, rna_to_protein_mean, rna_to_protein_std, 'RNA/protein')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
