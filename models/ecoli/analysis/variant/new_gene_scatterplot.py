"""
Plot mRNA and protein counts for each new gene, colored by variant, for all generations, early generations, and/or late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, index_of_first
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS

exclude_timeout_cells = 1 # 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_early_gens = 1 # 1 to plot early (before MIN_LATE_CELL_INDEX), and late generationss in addition to all generations

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index
MAX_CELL_INDEX = 8 # do not include any generation >= this index
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized
MAX_CELL_LENGTH = 180
if exclude_timeout_cells:
	MAX_CELL_LENGTH += 1000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def scatter(self, ax, new_gene_data,doubling_time_data, data_start, data_end, xlabel,ylabel, xlim=None,ylim=None, sf=1):
		for variant, variant_data in new_gene_data.items():
			variant_data = variant_data[data_start:data_end]
			dt_data = doubling_time_data[variant][data_start:data_end]

			color = COLORS[variant % len(COLORS)]
			mean = dt_data.mean()
			std = dt_data.std()
			mean_new_gene = variant_data.mean()
			ax.scatter(variant_data, dt_data, color=color, alpha=0.5,
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
		print("Running analysis script with exclude_timeout_cells=", exclude_timeout_cells, " and exclude_early_gens=", exclude_early_gens)

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'], monomer_sim_data['id']))
		new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
		assert len(new_gene_mRNA_ids) != 0, 'no new gene mRNAs found'
		assert len(new_gene_monomer_ids) != 0, 'no new gene proteins found'
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), 'number of new gene monomers and mRNAs should be equal'

		# Data extraction
		print("---Data Extraction---")
		doubling_times = {}
		new_gene_mRNA_counts = [{} for id in new_gene_mRNA_ids]
		new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]

		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:
			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()

			exclude_timeout_index = index_of_first(dt,MAX_CELL_LENGTH)
			doubling_times[variant] = dt[:exclude_timeout_index]

			# New gene mRNA and monomer counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract mRNA indexes for each new gene
				mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
				mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
				new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids] # need to add the [c] to it?

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
							   enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]

			avg_new_gene_mRNA_counts = read_stacked_columns(all_cells, 'mRNACounts', 'mRNA_counts',fun=lambda x: np.mean(x[:,new_gene_mRNA_indexes],axis=0))
			avg_new_gene_monomer_counts = read_stacked_columns(all_cells, 'MonomerCounts', 'monomerCounts',fun=lambda x: np.mean(x[:,new_gene_monomer_indexes],axis=0))

			avg_new_gene_mRNA_counts = avg_new_gene_mRNA_counts[dt < MAX_CELL_LENGTH]
			avg_new_gene_monomer_counts = avg_new_gene_monomer_counts[dt < MAX_CELL_LENGTH]

			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_counts[i][variant] = np.log10(avg_new_gene_mRNA_counts[:,i] + 1)
				new_gene_monomer_counts[i][variant] = np.log10(avg_new_gene_monomer_counts[:,i] + 1)

		# Plotting
		print("---Plotting---")
		std_sf = 2
		std_xlim = [-1,8]
		std_ylim = [30,185]
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
			std_monomer_xlab = 'Log10(' + new_gene_monomer_ids[i][:-3] + ' Counts + 1)'
			for j in range(len(data_start)):
				_, axes = plt.subplots(2, 1, figsize=(10, 10))
				self.scatter(axes[0], new_gene_monomer_counts[i], doubling_times, data_start[j], data_end[j], std_monomer_xlab, std_ylab, sf=std_sf,
							 xlim=std_xlim, ylim=std_ylim)
				self.scatter(axes[1], new_gene_mRNA_counts[i], doubling_times, data_start[j], data_end[j], std_mRNA_xlab, std_ylab, sf=std_sf,
							 xlim=std_xlim, ylim=std_ylim)
				plt.tight_layout()
				exportFigure(plt, plotOutDir, plotOutFileName+plot_label[j]+new_gene_monomer_ids[i][:-3], metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
