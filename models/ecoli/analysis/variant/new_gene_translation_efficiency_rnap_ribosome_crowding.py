"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Plots:
- Number of TUs overcrowded by RNA polymerases
- Number of Monomers overcrowded by ribosomes
- TODO: Add heatmap for number of time steps new genes are overcrowded on avg (generation, seed) by RNA polymerases and ribsosomes

Here, overcrowded is defined as the actual probability being less than the
target probability transcript on average for at least one generation in at
least one seed for that variant index.
"""

### TODO: filter early vs late gens, sims that timed out, etc

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
from models.ecoli.sim.variants.new_gene_expression_and_translation_efficiency \
	import NEW_GENE_EXPRESSION_FACTORS, NEW_GENE_TRANSLATION_EFFICIENCY_VALUES

# 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_timeout_cells = 1

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE=9
MAX_VARIANT = 43 # do not include any variant >= this index
MAX_CELL_INDEX = 16 # do not include any generation >= this index

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

MAX_CELL_LENGTH = 18000
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

MAX_YLIM_PLOT = MAX_CELL_LENGTH

"""
Count number of sims that reach this generation (remember index 7 
corresponds to generation 8)
"""
COUNT_INDEX = 15

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def heatmap(self, ax, mask, data, completion_data, xlabel, ylabel, xlabels,
				ylabels, title):
		im = ax.imshow(data, cmap="GnBu")
		ax.set_xticks(np.arange(len(xlabels)))
		ax.set_xticklabels(xlabels)
		ax.set_yticks(np.arange(len(
			ylabels)))
		ax.set_yticklabels(ylabels)
		plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
				 rotation_mode="anchor")
		for i in range(len(ylabels)):
			for j in range(len(xlabels)):
				if mask[i,j]:
					col = "k"
					if completion_data[i,j] < 0.9:
						col = "r"
					text = ax.text(j, i, data[i, j],
								   ha="center", va="center", color=col)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
		ax.set_title(title)

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# TODO: READ IN EXP AND TRL EFF LISTS FROM SIM METADATA INSTEAD OF IMPORT
		SEPARATOR = len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES)

		reached_count_gen = {}
		generations = {}

		variants = self.ap.get_variants()
		variant_index_to_values = {}
		variant_index_to_list_indices = {}
		variant_mask = np.zeros((  # Track whether we ran this sim
			len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES),
			len(NEW_GENE_EXPRESSION_FACTORS)), dtype=bool)

		# Match variant index to new gene expression and translation
		# efficicency values
		for index in variants:
			if index >= MAX_VARIANT:
				continue

			if index == 0:
				expression_list_index = 0
				trl_eff_list_index = len(
					NEW_GENE_TRANSLATION_EFFICIENCY_VALUES) - 1
				expression_variant_index = 0
				# Note: this value should not matter since gene is knocked out
				trl_eff_value = 0
			else:
				trl_eff_list_index = index % SEPARATOR
				if trl_eff_list_index == 0:
					expression_list_index = index // SEPARATOR
				else:
					expression_list_index = index // SEPARATOR + 1

				expression_variant_index = NEW_GENE_EXPRESSION_FACTORS[
											   expression_list_index]
				trl_eff_value = NEW_GENE_TRANSLATION_EFFICIENCY_VALUES[
					trl_eff_list_index]
			variant_index_to_values[index] = np.array([
				expression_variant_index, trl_eff_value])
			variant_index_to_list_indices[index] = np.array([
				expression_list_index, trl_eff_list_index])
			variant_mask[trl_eff_list_index, expression_list_index] = True

		# Create data structures that we will use for the heatmaps
		completed_gens_heatmap = np.zeros(
			(1, len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES),
			 len(NEW_GENE_EXPRESSION_FACTORS)))
		rnap_crowding_heatmap = np.zeros((3,
			len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES),
			len(NEW_GENE_EXPRESSION_FACTORS))) - 1
		ribosome_crowding_heatmap = np.zeros((3,
			len(NEW_GENE_TRANSLATION_EFFICIENCY_VALUES),
			len(NEW_GENE_EXPRESSION_FACTORS))) - 1

		# Data extraction
		print("---Data Extraction---")
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
				cell_path))[-6:]) for cell_path in all_cells])[
				exclude_timeout_cell_mask]
			generations[variant] = all_cells_gens

			# Count the number of simulations that reach gen COUNT_INDEX + 1
			num_count_gen = len(self.ap.get_cells(variant=[variant],
												  generation=[COUNT_INDEX],
												  only_successful=True))
			num_zero_gen = len(self.ap.get_cells(variant=[variant],
												 generation=[0],
												 only_successful=True))
			reached_count_gen[variant] = num_count_gen / num_zero_gen

			# RNA polymerase overcrowding
			avg_actual_rna_synth_prob = read_stacked_columns(all_cells,
				'RnaSynthProb', 'actual_rna_synth_prob',
				fun=lambda x: np.mean(x, axis=0))
			avg_target_rna_synth_prob = read_stacked_columns(all_cells,
				'RnaSynthProb', 'target_rna_synth_prob',
				fun=lambda x: np.mean(x, axis=0))

			# Get indexes of transcription units that on
			# average were overcrowded in any generation for any seed
			avg_overcrowded_tu_indexes = np.where(
				sum(avg_actual_rna_synth_prob <
					avg_target_rna_synth_prob) > 0)[0]
			n_overcrowded_tus = len(avg_overcrowded_tu_indexes)

			# Ribosome overcrowding
			avg_actual_prob_translation_per_transcript = read_stacked_columns(all_cells,
				'RibosomeData', 'actual_prob_translation_per_transcript',
				fun=lambda x: np.mean(x, axis = 0))
			avg_target_prob_translation_per_transcript = read_stacked_columns(
				all_cells,
				'RibosomeData', 'target_prob_translation_per_transcript',
				fun=lambda x: np.mean(x, axis = 0))

			# Get indexes of proteins corresponding to mRNAs that on
			# average were overcrowded in any generation for any seed
			avg_overcrowded_monomer_indexes = np.where(
				sum(avg_actual_prob_translation_per_transcript <
				avg_target_prob_translation_per_transcript) > 0)[0]
			n_overcrowded_monomers = len(avg_overcrowded_monomer_indexes)

			# Add values to heatmap data structures
			exp_index, trl_eff_index = variant_index_to_list_indices[
				variant]
			completed_gens_heatmap[0, trl_eff_index, exp_index] = \
				round(reached_count_gen[variant], 2)
			rnap_crowding_heatmap[0, trl_eff_index, exp_index] = \
				n_overcrowded_tus
			ribosome_crowding_heatmap[0, trl_eff_index, exp_index] = \
				n_overcrowded_monomers

			if exclude_early_gens == 1:
				# Add early gen values to the heatmap structure
				early_cell_mask = generations[variant] < MIN_LATE_CELL_INDEX
				if len(early_cell_mask) == 1:
					early_cell_mask = early_cell_mask[0]

				rnap_crowding_heatmap[1, trl_eff_index, exp_index] = \
					len(np.where(sum((avg_actual_rna_synth_prob <
					avg_target_rna_synth_prob)[early_cell_mask, :]) > 0)[0])
				ribosome_crowding_heatmap[1, trl_eff_index, exp_index] = \
					len(np.where(sum((avg_actual_prob_translation_per_transcript <
					avg_target_prob_translation_per_transcript)[
					early_cell_mask,:]) > 0)[0])

				# Add late gen values to the heatmap structure
				late_cell_mask = np.logical_and((generations[variant] >=
								  MIN_LATE_CELL_INDEX),
								 (generations[variant] < MAX_CELL_INDEX))
				if len(late_cell_mask) == 1:
					late_cell_mask = late_cell_mask[0]
				if sum(late_cell_mask) != 0:
					rnap_crowding_heatmap[2, trl_eff_index, exp_index] = \
						len(np.where(sum((avg_actual_rna_synth_prob <
						avg_target_rna_synth_prob)[late_cell_mask, :]) > 0)[0])
					ribosome_crowding_heatmap[2, trl_eff_index, exp_index] = \
						len(np.where(sum((avg_actual_prob_translation_per_transcript <
						avg_target_prob_translation_per_transcript)[
						late_cell_mask,:]) > 0)[0])

		# Plotting
		print("---Plotting---")
		plot_descr = ["_all_gens"]
		if exclude_early_gens == 1:
			plot_descr += ["_early_gens", "_late_gens"]

		for j in range(len(plot_descr)):

			# RNA Polymerase Crowding
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask,
						 rnap_crowding_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 NEW_GENE_EXPRESSION_FACTORS,
						 NEW_GENE_TRANSLATION_EFFICIENCY_VALUES,
						 "RNA Polymerase Crowding: # of TUs")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir,
						 'rnap_crowding_heatmap' +
						 plot_descr[j])

			# Ribosome Crowding
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask,
						 ribosome_crowding_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 NEW_GENE_EXPRESSION_FACTORS,
						 NEW_GENE_TRANSLATION_EFFICIENCY_VALUES,
						 "Ribosome Crowding: # of Monomers")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir,
						 'ribosome_crowding_heatmap' +
						 plot_descr[j])

		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
