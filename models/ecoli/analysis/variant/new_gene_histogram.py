"""
Plot mRNA and protein counts for new genes
"""

import numpy as np
from matplotlib import pyplot as plt

import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS


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

		def downsample(x):
			"""Average every n_downsample points to one value to smooth and downsample"""
			n_downsample = 100
			if (extra_points := x.shape[0] % n_downsample) != 0:
				x = x[:-extra_points]
			return x.reshape(-1, n_downsample).mean(1).reshape(-1, 1)

		# Data extraction
		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt < MAX_CELL_LENGTH]
			growth_rates[variant] = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'Mass', 'instantaneous_growth_rate',
				remove_first=True, fun=downsample).squeeze() * 3600.

			all_mRNA_counts = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'mRNACounts', 'mRNA_counts')
			all_monomer_counts = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'MonomerCounts', 'monomerCounts')
			if variant == min_variant: ### TODO flag new gene mRNAs and proteins more efficiently
				# Extract mRNA indexes for each new gene
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')
				mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
				mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
				new_gene_mRNA_ids = [k for k, v in mRNA_idx.items() if k.startswith('NG')]
				new_gene_mRNA_indexes = [v for k, v in mRNA_idx.items() if k.startswith('NG')]
				assert len(new_gene_mRNA_ids) != 0, 'no new gene mRNAs found'

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				monomer_idx = {monomer: i for i, monomer in
							   enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_ids = [k for k, v in monomer_idx.items() if k.startswith('NG')]
				new_gene_monomer_indexes = [v for k, v in monomer_idx.items() if k.startswith('NG')]
				assert len(new_gene_monomer_ids) != 0, 'no new gene proteins found'

				assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), 'number of new gene monomers and mRNAs should be equal'

				new_gene_mRNA_counts = [{} for id in new_gene_mRNA_ids]
				new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]

			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_counts_var = all_mRNA_counts[:,new_gene_mRNA_indexes[i]]
				new_gene_mRNA_counts_var = downsample(new_gene_mRNA_counts_var)
				new_gene_mRNA_counts_var[new_gene_mRNA_counts_var == 0] = 10**(-10)
				new_gene_mRNA_counts[i][variant] = np.log10(new_gene_mRNA_counts_var)

				new_gene_monomer_counts_var = all_monomer_counts[:, new_gene_monomer_indexes[i]]
				new_gene_monomer_counts_var = downsample(new_gene_monomer_counts_var)
				new_gene_monomer_counts_var[new_gene_monomer_counts_var == 0] = 10 ** (-10)
				new_gene_monomer_counts[i][variant] = np.log10(new_gene_monomer_counts_var)

		# Plotting
		for i in range(len(new_gene_mRNA_ids)):
			_, axes = plt.subplots(2, 1, figsize=(10, 10))
			self.hist(axes[0], new_gene_monomer_counts[i], 'Log10(' + new_gene_monomer_ids[i][:-3] +' Counts)')
			self.hist(axes[1], new_gene_mRNA_counts[i], 'Log10(' + new_gene_mRNA_ids[i][:-3] +' Counts)')
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName+'_'+new_gene_monomer_ids[i][:-3], metadata)

			axes[0].set_xlim([-11, 8])
			axes[1].set_xlim([-11, 8])
			exportFigure(plt, plotOutDir, plotOutFileName+'_'+new_gene_monomer_ids[i][:-3] +'_trimmed', metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
