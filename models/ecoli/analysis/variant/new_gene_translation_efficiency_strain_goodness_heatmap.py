### TODO: Merge into new_gene_translation_efficiency_heatmaps.py (kept
# separate now for convenience of running on large simulation set)

"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Plots:
- Average "strain goodness": avg new gene protein counts from the final
generation, multiplied by the growth rate (log_e(2)/doubling time)

TODO:
- Accomodate more than one new gene
"""

import numpy as np
from matplotlib import pyplot as plt
import math

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, \
	read_stacked_columns, stacked_cell_threshold_mask
from wholecell.analysis.plotting_tools import heatmap

import os.path
import pickle

# 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_timeout_cells = 1

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE=9
MAX_VARIANT = 55 # do not include any variant >= this index
MAX_CELL_INDEX = 16 # do not include any generation >= this index

"""
Count number of sims that reach this generation (remember index 7 
corresponds to generation 8)
"""
COUNT_INDEX = 15

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

MAX_CELL_LENGTH = 36000
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		print("Running analysis script with exclude_timeout_cells=",
			  exclude_timeout_cells, " and exclude_early_gens=",
			  exclude_early_gens)

		# Map variant indices to expression factors and translation efficiency
		# values
		if 'new_gene_expression_factors' not in metadata or \
				'new_gene_translation_efficiency_values' not in metadata:
			print("This plot is intended to be run on simulations where the"
				  " new gene expression-translation efficiency variant was "
				  "enabled, but no parameters for this variant were found.")

		new_gene_expression_factors = metadata[
			'new_gene_expression_factors']
		new_gene_translation_efficiency_values = metadata[
			'new_gene_translation_efficiency_values']

		separator = len(new_gene_translation_efficiency_values)

		variants = self.ap.get_variants()
		variant_index_to_values = {}
		variant_index_to_list_indices = {}
		variant_mask = np.zeros(( # Track whether we ran this sim
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)), dtype=bool)

		for index in variants:
			if index >= MAX_VARIANT:
				continue

			if index == 0:
				expression_list_index = 0
				trl_eff_list_index = len(
					new_gene_translation_efficiency_values) - 1
				expression_variant_index = 0
				# Note: this value should not matter since gene is knocked out
				trl_eff_value = 0
			else:
				trl_eff_list_index = index % separator
				if trl_eff_list_index == 0:
					expression_list_index = index // separator
				else:
					expression_list_index = index // separator + 1

				expression_variant_index = new_gene_expression_factors[
											   expression_list_index]
				trl_eff_value = new_gene_translation_efficiency_values[
					trl_eff_list_index]
			variant_index_to_values[index] = np.array([
				expression_variant_index, trl_eff_value])
			variant_index_to_list_indices[index] = np.array([
				expression_list_index, trl_eff_list_index])
			variant_mask[trl_eff_list_index, expression_list_index] = True

		# Create data structures that we will use for the heatmaps
		completed_gens_heatmap = np.zeros(( 1,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)))
		# TODO: Expand to Accomodate Multiple New Genes
		avg_new_gene_strain_goodness_heatmap = np.zeros((1,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1

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
		print("---Data Extraction---")
		doubling_times = {}
		reached_count_gen = {}
		generations = {}
		new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]
		new_gene_strain_goodness = [{} for id in new_gene_monomer_ids]

		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant],
										  only_successful=True)
			if len(all_cells) == 0:
				continue

			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', MAX_CELL_LENGTH,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			all_cells_gens = np.array([int(os.path.basename(os.path.dirname(
				cell_path))[-6:]) for cell_path in all_cells])[exclude_timeout_cell_mask]
			generations[variant] = all_cells_gens

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
									  fun=lambda x: (x[-1] - x[0]) / 60.)
			doubling_times[variant] = dt[exclude_timeout_cell_mask]

			# Count the number of simulations that reach gen COUNT_INDEX + 1
			num_count_gen = len(self.ap.get_cells(variant=[variant],
							  generation = [COUNT_INDEX],
												  only_successful=True))
			num_zero_gen = len(self.ap.get_cells(variant=[variant],
							  generation = [0], only_successful=True))
			reached_count_gen[variant] = num_count_gen / num_zero_gen

			# New gene mRNA and monomer counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(
					simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
									enumerate(
										monomer_counts_reader.readAttribute(
											'monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id)
											for monomer_id in
											new_gene_monomer_ids]

			avg_new_gene_monomer_counts = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts', fun=lambda x:
				np.mean( x[:,new_gene_monomer_indexes], axis=0))

			avg_new_gene_monomer_counts = \
				avg_new_gene_monomer_counts[exclude_timeout_cell_mask,]

			# Strain goodness calculation
			count_gen_mask = (generations[variant] == COUNT_INDEX)[exclude_timeout_cell_mask]
			for i in range(len(new_gene_monomer_ids)):
				new_gene_monomer_counts[i][variant] = \
					avg_new_gene_monomer_counts[:, i]
				# TODO: resolve this edge case in a nicer way
				if len(new_gene_monomer_counts[i][variant].shape) == 2:
					count_gen_mask = count_gen_mask[0]
				avg_new_gene_monomer_counts_for_strain_goodness = \
					avg_new_gene_monomer_counts[:, i][count_gen_mask]
				dt_for_strain_goodness = doubling_times[variant][
					count_gen_mask].squeeze()
				new_gene_strain_goodness[i][variant] = \
					avg_new_gene_monomer_counts_for_strain_goodness\
					/dt_for_strain_goodness * math.log(2)
				print(avg_new_gene_monomer_counts_for_strain_goodness)
				print(dt_for_strain_goodness)
				print(new_gene_strain_goodness[i][variant])

			# Add values to heatmap data structures
			exp_index, trl_eff_index = variant_index_to_list_indices[variant]
			completed_gens_heatmap[0, trl_eff_index, exp_index] = \
				round(reached_count_gen[variant],2)
			i = 0 ### TODO: accomodate multiple new genes
			if len(new_gene_strain_goodness[i][variant]) != 0:
				avg_new_gene_strain_goodness_heatmap[0, trl_eff_index, exp_index] \
					= round(np.mean(new_gene_strain_goodness[i][variant]))
		
		# Plotting
		print("---Plotting---")
		plot_descr = ["_all_gens"]

		# New Gene Monomer Counts
		fig, ax = plt.subplots(1, 1, figsize=(10, 5))
		heatmap(self, ax, variant_mask,
					 avg_new_gene_strain_goodness_heatmap[0, :, :],
					 completed_gens_heatmap[0, :, :],
					 "Expression Variant",
					 "Translation Efficiency Value (Normalized)",
					 new_gene_expression_factors,
					 new_gene_translation_efficiency_values,
					 "Strain Goodness Metric", "x-small")
		fig.tight_layout()
		plt.show()
		exportFigure(plt, plotOutDir, 'new_gene_strain_goodness' + plot_descr[0])

		plt.close('all')

if __name__ == "__main__":
	Plot().cli()