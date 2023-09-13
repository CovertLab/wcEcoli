"""
Plot (histogram and scatterplot) mRNA and protein counts for each new gene,
colored by variant, for all generations, early generations, and/or late
(i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure,\
	read_stacked_columns, stacked_cell_threshold_mask
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS\
	as COLORS, labeled_indexable_hist, labeled_indexable_scatter

# 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_timeout_cells = 1

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index
MAX_CELL_INDEX = 8 # do not include any generation >= this index

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

MAX_CELL_LENGTH = 180
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

MAX_YLIM_PLOT = MAX_CELL_LENGTH


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# print("Running analysis script with exclude_timeout_cells=",
		# 	  exclude_timeout_cells,
		# 	  " and exclude_early_gens=", exclude_early_gens)

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in new_gene_mRNA_ids]
		if len(new_gene_mRNA_ids) == 0:
			print("This plot is intended to be run on simulations where the"
				  " new gene option was enabled, but no new gene mRNAs were "
				  "found.")
			return
		if len(new_gene_monomer_ids) == 0:
			print("This plot is intended to be run on simulations where the "
				  "new gene option was enabled, but no new gene proteins "
				  "were "
				  "found.")
			return
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids),\
			'number of new gene monomers and mRNAs should be equal'

		# Data extraction
		doubling_times = {}
		generations = {}
		new_gene_mRNA_counts = [{} for id_ in new_gene_mRNA_ids]
		new_gene_monomer_counts = [{} for id_ in new_gene_monomer_ids]

		variants = self.ap.get_variants()
		min_variant = min(variants)

		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			all_cells = self.ap.get_cells(variant=[variant],
										  only_successful=True)
			if len(all_cells) == 0:
				continue

			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', MAX_CELL_LENGTH,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			all_cells_gens = np.array([int(os.path.basename(os.path.dirname(
				cell_path))[-6:]) for cell_path in all_cells])[
				exclude_timeout_cell_mask]
			generations[variant] = all_cells_gens

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
									  fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[exclude_timeout_cell_mask]

			# New gene mRNA and monomer counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract mRNA indexes for each new gene
				mRNA_counts_reader = TableReader(os.path.join(simOutDir,
															  'RNACounts'))
				mRNA_idx_dict = {rna[:-3]: i for i, rna in
					enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
				new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id)
										 for mRNA_id in new_gene_mRNA_ids]

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(
					simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
					enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id)
					for monomer_id in new_gene_monomer_ids]

			avg_new_gene_mRNA_counts = read_stacked_columns(all_cells,
				'RNACounts', 'mRNA_counts',
				fun=lambda x: np.mean(x[:,new_gene_mRNA_indexes],axis=0))
			avg_new_gene_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				fun=lambda x: np.mean(x[:,new_gene_monomer_indexes],axis=0))

			avg_new_gene_mRNA_counts = \
				avg_new_gene_mRNA_counts[exclude_timeout_cell_mask,]
			avg_new_gene_monomer_counts = \
				avg_new_gene_monomer_counts[exclude_timeout_cell_mask,]

			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_counts[i][variant] = \
					np.log10(avg_new_gene_mRNA_counts[:,i] + 1)
				new_gene_monomer_counts[i][variant] = \
					np.log10(avg_new_gene_monomer_counts[:,i] + 1)

		# Plotting
		std_bin_width = 0.25
		std_sf = 2
		std_xlim = [-1,8]
		std_ylim = [30,MAX_YLIM_PLOT+10]
		std_ylab = 'Doubling Time (min)'

		data_start = [0]
		data_end = [MAX_CELL_INDEX]
		plot_label = ['_all_gens_']
		if exclude_early_gens:
			data_start += [0,MIN_LATE_CELL_INDEX]
			data_end += [MIN_LATE_CELL_INDEX,MAX_CELL_INDEX]
			plot_label += ['_early_gens_', '_late_gens_']

		for i in range(len(new_gene_mRNA_ids)):
			std_mRNA_xlab = 'Log10(' + new_gene_mRNA_ids[i] + ' Counts + 1)'
			std_monomer_xlab = 'Log10(' + new_gene_monomer_ids[i][:-3] \
							   + ' Counts + 1)'
			for j in range(len(data_start)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				labeled_indexable_scatter(self, axes[0],
							new_gene_monomer_counts[i],
							doubling_times, generations, data_start[j],
							data_end[j], COLORS, std_monomer_xlab, std_ylab,
							sf=std_sf, xlim=std_xlim, ylim=std_ylim)
				labeled_indexable_scatter(self, axes[1],
							new_gene_mRNA_counts[i],
							doubling_times, generations, data_start[j],
							data_end[j], COLORS, std_mRNA_xlab, std_ylab,
							sf=std_sf, xlim=std_xlim, ylim=std_ylim)
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName+plot_label[j]+
							 'scatterplot_' + new_gene_monomer_ids[i][:-3],
							 metadata)

				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				labeled_indexable_hist(
					self, axes[0], new_gene_monomer_counts[i], generations,
					data_start[j], data_end[j], COLORS, std_monomer_xlab,
					bin_width=std_bin_width, sf=std_sf,
					xlim=std_xlim)
				labeled_indexable_hist(
					self, axes[1], new_gene_mRNA_counts[i], generations,
					data_start[j], data_end[j], COLORS, std_mRNA_xlab,
					bin_width=std_bin_width, sf=std_sf, xlim=std_xlim)
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName + plot_label[j]
							 + 'histogram_' + new_gene_monomer_ids[i][:-3],
							 metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
