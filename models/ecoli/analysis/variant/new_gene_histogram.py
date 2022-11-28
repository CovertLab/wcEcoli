"""
Plot histogram of mRNA and protein counts for new genes, colored by variant
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
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def hist(self, ax, data, xlabel, bin_width=1., xlim=None, sf=1):
		if xlim:
			bins = np.histogram(range(xlim[0],xlim[1]+1), bins=int(np.ceil((xlim[1]-xlim[0])/bin_width)))[1]

		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			if not xlim:
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
		# Data extraction
		print("---Data Extraction---")
		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			if len(all_cells) >= MIN_LATE_CELL_INDEX:
				early_cell_index = list(range(MIN_LATE_CELL_INDEX))
				late_cell_index = list(range(MIN_LATE_CELL_INDEX, len(all_cells)))

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

				new_gene_mRNA_counts_early_gens = [{} for id in new_gene_mRNA_ids]
				new_gene_monomer_counts_early_gens = [{} for id in new_gene_monomer_ids]

				new_gene_mRNA_counts_late_gens = [{} for id in new_gene_mRNA_ids]
				new_gene_monomer_counts_late_gens = [{} for id in new_gene_monomer_ids]

			print("Variant: ",variant)
			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_counts_var = np.ones(len(all_cells))
				new_gene_monomer_counts_var = np.ones(len(all_cells))
				for gen in range(len(all_cells)):
					sim_dir = all_cells[gen]
					simOutDir = os.path.join(sim_dir, 'simOut')
					mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
					monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))

					mRNA_counts_gen = mRNA_counts_reader.readColumn('mRNA_counts')
					new_gene_mRNA_counts_gen = mRNA_counts_gen[:, new_gene_mRNA_indexes[i]]
					new_gene_mRNA_counts_var[gen] = np.mean(new_gene_mRNA_counts_gen)

					monomer_counts_gen = monomer_counts_reader.readColumn('monomerCounts')
					new_gene_monomer_counts_gen = monomer_counts_gen[:, new_gene_monomer_indexes[i]]
					new_gene_monomer_counts_var[gen] = np.mean(new_gene_monomer_counts_gen)

				new_gene_mRNA_counts[i][variant] = np.log10(new_gene_mRNA_counts_var + 1)
				new_gene_monomer_counts[i][variant] = np.log10(new_gene_monomer_counts_var + 1)

				if len(all_cells) >= MIN_LATE_CELL_INDEX:
					new_gene_mRNA_counts_early_gens[i][variant] = new_gene_mRNA_counts[i][variant][early_cell_index]
					new_gene_monomer_counts_early_gens[i][variant] = new_gene_monomer_counts[i][variant][early_cell_index]

					new_gene_mRNA_counts_late_gens[i][variant] = new_gene_mRNA_counts[i][variant][late_cell_index]
					new_gene_monomer_counts_late_gens[i][variant] = new_gene_monomer_counts[i][variant][late_cell_index]

		# Plotting
		print("---Plotting---")
		# ALL GENS
		for i in range(len(new_gene_mRNA_ids)):
			_, axes = plt.subplots(2, 1, figsize=(10, 10))
			self.hist(axes[0], new_gene_monomer_counts[i], 'Log10(' + new_gene_monomer_ids[i][:-3] +' Counts + 1)', bin_width=0.25, sf=2,xlim=[-1,8])
			self.hist(axes[1], new_gene_mRNA_counts[i], 'Log10(' + new_gene_mRNA_ids[i][:-3] +' Counts + 1)', bin_width=0.25, sf=2,xlim=[-1,8])
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName+'_all_gens_'+new_gene_monomer_ids[i][:-3], metadata)

		if len(all_cells) >= MIN_LATE_CELL_INDEX:
			# EARLY GENS
			for i in range(len(new_gene_mRNA_ids)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.hist(axes[0], new_gene_monomer_counts_early_gens[i], 'Log10(' + new_gene_monomer_ids[i][:-3] + ' Counts + 1)',
						  bin_width=0.25, sf=2, xlim=[-1, 8])
				self.hist(axes[1], new_gene_mRNA_counts_early_gens[i], 'Log10(' + new_gene_mRNA_ids[i][:-3] + ' Counts + 1)',
						  bin_width=0.25, sf=2, xlim=[-1, 8])
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_early_gens_' + new_gene_monomer_ids[i][:-3], metadata)

			# LATE GENS
			for i in range(len(new_gene_mRNA_ids)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.hist(axes[0], new_gene_monomer_counts_late_gens[i],
						  'Log10(' + new_gene_monomer_ids[i][:-3] + ' Counts + 1)',
						  bin_width=0.25, sf=2, xlim=[-1, 8])
				self.hist(axes[1], new_gene_mRNA_counts_late_gens[i],
						  'Log10(' + new_gene_mRNA_ids[i][:-3] + ' Counts + 1)',
						  bin_width=0.25, sf=2, xlim=[-1, 8])
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_late_gens_' + new_gene_monomer_ids[i][:-3], metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
