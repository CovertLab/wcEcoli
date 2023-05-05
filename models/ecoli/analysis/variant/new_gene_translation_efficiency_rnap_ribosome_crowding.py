"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Plots:
- Number of TUs overcrowded by RNA polymerases
- Number of Monomers overcrowded by ribosomes
- Number of time steps new genes are overcrowded on avg (generation,
seed) by RNA polymerases and ribsosomes

Here, overcrowded is defined as the actual probability being less than the
target probability transcript on average for at least one generation in at
least one seed for that variant index.
"""

### TODO: filter sims that timed out, accomodate multiple new genes

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure,\
	read_stacked_columns, stacked_cell_threshold_mask
from wholecell.analysis.plotting_tools import heatmap

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

		reached_count_gen = {}
		generations = {}

		variants = self.ap.get_variants()
		min_variant = min(variants)
		variant_index_to_values = {}
		variant_index_to_list_indices = {}
		variant_mask = np.zeros((  # Track whether we ran this sim
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)), dtype=bool)

		# Match variant index to new gene expression and translation
		# efficicency values
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
		completed_gens_heatmap = np.zeros(
			(1, len(new_gene_translation_efficiency_values),
			 len(new_gene_expression_factors)))
		rnap_crowding_heatmap = np.zeros((3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		ribosome_crowding_heatmap = np.zeros((3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		# TODO: Expand to Accomodate Multiple New Genes
		avg_time_new_gene_rnap_overcrowded_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		avg_time_new_gene_ribosome_overcrowded_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']][
			'id'].tolist()
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
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
			'number of new gene monomers and mRNAs should be equal'

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

			# New gene mRNA and monomer counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract mRNA indexes for each new gene
				mRNA_counts_reader = TableReader(os.path.join(simOutDir,
															  'RNACounts'))
				mRNA_idx_dict = {rna[:-3]: i for i, rna in
								 enumerate(mRNA_counts_reader.readAttribute(
									 'mRNA_ids'))}
				new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id)
										 for mRNA_id in new_gene_mRNA_ids]

				# Extract RNA indexes for each new gene
				rnap_reader = TableReader(os.path.join(simOutDir,
															  'RnaSynthProb'))
				RNA_idx_dict = {rna[:-3]: i for i, rna in
								 enumerate(rnap_reader.readAttribute(
									 'rnaIds'))}
				new_gene_RNA_indexes = [RNA_idx_dict.get(mRNA_id)
										 for mRNA_id in new_gene_mRNA_ids]

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

			# Average fraction of time steps that RNA polymerase
			# overcrowding occurs for new genes per generation
			new_gene_num_time_steps_rnap_overcrowded = read_stacked_columns(
				all_cells, 'RnaSynthProb', 'tu_is_overcrowded',
				fun=lambda x: np.sum(x[:,new_gene_RNA_indexes],
				axis=0) / (x[:,new_gene_RNA_indexes].shape[0]))

			# Average fraction of time steps that ribosome
			# overcrowding occurs for new genes per generation
			new_gene_num_time_steps_ribosome_overcrowded = read_stacked_columns(
				all_cells, 'RibosomeData', 'mRNA_is_overcrowded',
				fun=lambda x: np.sum(x[:, new_gene_monomer_indexes],
				axis=0)/(x[:,new_gene_monomer_indexes].shape[0]))

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

			i = 0 ### TODO: accomodate multiple new genes
			avg_time_new_gene_rnap_overcrowded_heatmap[0, trl_eff_index, exp_index] =\
				round(np.mean(new_gene_num_time_steps_rnap_overcrowded[i,:]), 2)
			avg_time_new_gene_ribosome_overcrowded_heatmap[0, trl_eff_index, exp_index] =\
				round(np.mean(new_gene_num_time_steps_ribosome_overcrowded[i, :]), 2)

			rnap_crowding_heatmap[0, trl_eff_index, exp_index] = \
				n_overcrowded_tus
			ribosome_crowding_heatmap[0, trl_eff_index, exp_index] = \
				n_overcrowded_monomers

			if exclude_early_gens == 1:
				# Add early gen values to the heatmap structure
				early_cell_mask = generations[variant] < MIN_LATE_CELL_INDEX
				if len(early_cell_mask) == 1:
					early_cell_mask = early_cell_mask[0]

				i = 0  ### TODO: accomodate multiple new genes
				avg_time_new_gene_rnap_overcrowded_heatmap[
					1, trl_eff_index, exp_index] = \
					round(np.mean(new_gene_num_time_steps_rnap_overcrowded[:,
						i][early_cell_mask]), 2)
				avg_time_new_gene_ribosome_overcrowded_heatmap[
					1, trl_eff_index, exp_index] = \
					round(np.mean(new_gene_num_time_steps_ribosome_overcrowded[:,
						i][early_cell_mask]),2)

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
					i = 0  ### TODO: accomodate multiple new genes
					avg_time_new_gene_rnap_overcrowded_heatmap[
						2, trl_eff_index, exp_index] = \
						round(np.mean(
							new_gene_num_time_steps_rnap_overcrowded[:,
							i][late_cell_mask]), 2)
					avg_time_new_gene_ribosome_overcrowded_heatmap[
						2, trl_eff_index, exp_index] = \
						round(np.mean(
							new_gene_num_time_steps_ribosome_overcrowded[:,
							i][late_cell_mask]), 2)
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
			# New Gene RNA Polymerase Crowding - Fraction of Time Steps
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			heatmap(self, ax, variant_mask,
						 avg_time_new_gene_rnap_overcrowded_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Fraction of Time RNA Polymerase Crowded New Gene")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'new_gene_rnap_crowding_heatmap' +
						 plot_descr[j])

			# New Gene Ribosome Crowding - Fraction of Time Steps
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			heatmap(self, ax, variant_mask,
						 avg_time_new_gene_ribosome_overcrowded_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Fraction of Time Ribosome Crowded New Gene")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir,
						 'new_gene_ribosome_crowding_heatmap' +
						 plot_descr[j])

			# RNA Polymerase Crowding
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			heatmap(self, ax, variant_mask,
						 rnap_crowding_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "RNA Polymerase Crowding: # of TUs")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir,
						 'rnap_crowding_heatmap' +
						 plot_descr[j])

			# Ribosome Crowding
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			heatmap(self, ax, variant_mask,
						 ribosome_crowding_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Ribosome Crowding: # of Monomers")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir,
						 'ribosome_crowding_heatmap' +
						 plot_descr[j])

		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
