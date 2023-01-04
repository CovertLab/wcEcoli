"""
Plot mRNA and protein counts for new genes, colored by variant, for all generations, early generations, and late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, stacked_cell_identification
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS


FONT_SIZE=9
MAX_CELL_LENGTH = 180
#MAX_CELL_LENGTH += 1 # comment out this line to filter sims that reach the max time of 180 min
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized
MAX_VARIANT = 10 # do not include any variant >= this index

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def scatter(self, ax, new_gene_data,doubling_time_data, xlabel,ylabel, xlim=None,ylim=None, sf=1):
		for variant, variant_data in new_gene_data.items():
			color = COLORS[variant % len(COLORS)]
			mean = doubling_time_data[variant].mean()
			std = doubling_time_data[variant].std()
			mean_new_gene = variant_data.mean()
			std_new_gene = variant_data.std()
			ax.scatter(variant_data, doubling_time_data[variant], color=color, alpha=0.5,
				label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
			ax.scatter(mean_new_gene,mean,color=color,alpha=0.5,marker='x')

		if xlim:
			ax.set_xlim(xlim)
		if ylim:
			ax.set_ylim(ylim)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		doubling_times = {}
		doubling_times_early_gens = {}
		doubling_times_late_gens = {}

		# Data extraction
		print("---Data Extraction---")
		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:
			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt < MAX_CELL_LENGTH]

			if len(all_cells) >= MIN_LATE_CELL_INDEX:
				all_cells_gens = [int(c.split("/")[-2][-6:]) for c in all_cells]
				early_cell_index = [i for i,v in enumerate(all_cells_gens) if v < MIN_LATE_CELL_INDEX]
				late_cell_index = [i for i,v in enumerate(all_cells_gens) if v >= MIN_LATE_CELL_INDEX]

				dt_early_cells = dt[early_cell_index]
				dt_late_cells = dt[late_cell_index]

				doubling_times_early_gens[variant] = dt_early_cells[dt_early_cells < MAX_CELL_LENGTH ]
				doubling_times_late_gens[variant] = dt_late_cells[dt_late_cells < MAX_CELL_LENGTH ]

				### TODO: need some assert statement here for MIN_LATE_CELL_INDEX???

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

			avg_new_gene_mRNA_counts = read_stacked_columns(all_cells, 'mRNACounts', 'mRNA_counts',fun=lambda x: np.mean(x[:,new_gene_mRNA_indexes],axis=0))
			avg_new_gene_monomer_counts = read_stacked_columns(all_cells, 'MonomerCounts', 'monomerCounts',fun=lambda x: np.mean(x[:,new_gene_monomer_indexes],axis=0))

			for i in range(len(new_gene_mRNA_ids)): ### TODO: MORE EFFICIENT WAY TO DO THIS?
				new_gene_mRNA_counts[i][variant] = np.log10(avg_new_gene_mRNA_counts[:,i] + 1)
				new_gene_monomer_counts[i][variant] = np.log10(avg_new_gene_monomer_counts[:,i] + 1)

				new_gene_mRNA_counts_early_gens[i][variant] = new_gene_mRNA_counts[i][variant][early_cell_index]
				new_gene_monomer_counts_early_gens[i][variant] = new_gene_monomer_counts[i][variant][early_cell_index]

				new_gene_mRNA_counts_late_gens[i][variant] = new_gene_mRNA_counts[i][variant][late_cell_index]
				new_gene_monomer_counts_late_gens[i][variant] = new_gene_monomer_counts[i][variant][late_cell_index]

		# Plotting
		print("---Plotting---")
		for i in range(len(new_gene_mRNA_ids)):
			_, axes = plt.subplots(2, 1, figsize=(10, 10))
			self.scatter(axes[0], new_gene_monomer_counts[i], doubling_times, 'Log10(' + new_gene_monomer_ids[i][:-3] +' Counts + 1)', 'Doubling Time (min)', sf=2, xlim=[-1,8], ylim=[30, 185])
			self.scatter(axes[1], new_gene_mRNA_counts[i], doubling_times, 'Log10(' + new_gene_mRNA_ids[i][:-3] +' Counts + 1)', 'Doubling Time (min)', sf=2, xlim=[-1,8], ylim=[30, 185])
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName+'_all_gens_'+new_gene_monomer_ids[i][:-3], metadata)

		if len(all_cells) >= MIN_LATE_CELL_INDEX:
			# EARLY GENS
			for i in range(len(new_gene_mRNA_ids)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.scatter(axes[0], new_gene_monomer_counts_early_gens[i], doubling_times_early_gens,
							 'Log10(' + new_gene_monomer_ids[i][:-3] + ' Counts + 1)', 'Doubling Time (min)', sf=2,
							 xlim=[-1, 8], ylim=[30, 185])
				self.scatter(axes[1], new_gene_mRNA_counts_early_gens[i], doubling_times_early_gens,
							 'Log10(' + new_gene_mRNA_ids[i][:-3] + ' Counts + 1)', 'Doubling Time (min)', sf=2,
							 xlim=[-1, 8], ylim=[30, 185])
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_early_gens_' + new_gene_monomer_ids[i][:-3], metadata)

			# LATE GENS
			for i in range(len(new_gene_mRNA_ids)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.scatter(axes[0], new_gene_monomer_counts_late_gens[i], doubling_times_late_gens,
							 'Log10(' + new_gene_monomer_ids[i][:-3] + ' Counts + 1)', 'Doubling Time (min)', sf=2,
							 xlim=[-1, 8], ylim=[30, 185])
				self.scatter(axes[1], new_gene_mRNA_counts_late_gens[i], doubling_times_late_gens,
							 'Log10(' + new_gene_mRNA_ids[i][:-3] + ' Counts + 1)', 'Doubling Time (min)', sf=2,
							 xlim=[-1, 8], ylim=[30, 185])
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + '_late_gens_' + new_gene_monomer_ids[i][:-3],
							 metadata)

		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
