import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import COLORS_COLORBLIND as COLORS


FONT_SIZE=9
MAX_CELL_LENGTH = 180  # filter sims that reach the max time of 180 min


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def hist(self, ax, data, xlabel, bin_width=1., xlim=None, sf=1):
		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			bins = max(1, int(np.ceil((variant_data.max() - variant_data.min()) / bin_width)))
			mean = variant_data.mean()
			std = variant_data.std()
			ax.hist(variant_data, bins, color=color, alpha=0.5,
				label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
			ax.axvline(mean, color=color, linestyle='--', linewidth=1)

		if xlim:
			ax.set_xlim(xlim)
		self.remove_border(ax)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		doubling_times = {}
		growth_rates = {}
		new_gene_monomer_counts = {}
		new_gene_mRNA_counts = {}

		def downsample(x):
			"""Average every n_downsample points to one value to smooth and downsample"""
			n_downsample = 100
			if (extra_points := x.shape[0] % n_downsample) != 0:
				x = x[:-extra_points]
			return x.reshape(-1, n_downsample).mean(1).reshape(-1, 1)

		for variant in self.ap.get_variants():
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt < MAX_CELL_LENGTH]
			growth_rates[variant] = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'Mass', 'instantaneous_growth_rate',
				remove_first=True, fun=downsample).squeeze() * 3600.

			# TODO get indexes for New Gene mRNA and Protein

			all_mRNA_counts = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'mRNACounts', 'mRNA_counts')
			new_gene_mRNA_counts_var = all_mRNA_counts[:,-1]
			new_gene_mRNA_counts_var = downsample(new_gene_mRNA_counts_var)
			new_gene_mRNA_counts_var[new_gene_mRNA_counts_var == 0] = 10**(-10)
			new_gene_mRNA_counts[variant] = np.log10(new_gene_mRNA_counts_var)

			all_monomer_counts = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'MonomerCounts', 'monomerCounts')
			new_gene_monomer_counts_var = all_monomer_counts[:, -1]
			new_gene_monomer_counts_var = downsample(new_gene_monomer_counts_var)
			new_gene_monomer_counts_var[new_gene_monomer_counts_var == 0] = 10 ** (-10)
			new_gene_monomer_counts[variant] = np.log10(new_gene_monomer_counts_var)

		_, axes = plt.subplots(2, 1, figsize=(10, 10))

		self.hist(axes[0], new_gene_monomer_counts, 'Log10(New Gene Protein Counts)')
		self.hist(axes[1], new_gene_mRNA_counts, 'Log10(New Gene mRNA Counts)')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		axes[0].set_xlim([-10, 8])
		axes[1].set_xlim([-10, 8])
		exportFigure(plt, plotOutDir, plotOutFileName + '_trimmed', metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
