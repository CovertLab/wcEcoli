"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Plots:
- Average doubling time
- Percent of sims that successfully reached a given generation number
- Average new gene mRNA count
- Average new gene protein count
- Average new gene protein mass fraction
- Average number of ribosomes
- Average number of RNA polymerases
- Average ppGpp concentration

TODO:
- Accomodate more than one new gene
"""

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, \
	read_stacked_columns, stacked_cell_threshold_mask, \
	read_stacked_bulk_molecules, stacked_cell_identification
from unum.units import g, mol

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
MAX_VARIANT = 43 # do not include any variant >= this index
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
	### TODO: move to analysis_tools
	def heatmap(self, ax, mask, data, completion_data, xlabel, ylabel, xlabels,
				ylabels, title, textsize = "medium"):
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
								   ha="center", va="center", color=col,
								   fontsize=textsize)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
		ax.set_title(title)

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
		doubling_times_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		completed_gens_heatmap = np.zeros(( 1,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors)))
		avg_rnap_counts_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		avg_ribosome_counts_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		avg_ppgpp_counts_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		# TODO: Expand to Accomodate Multiple New Genes
		avg_new_gene_mRNA_counts_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		avg_new_gene_monomer_counts_heatmap = np.zeros(( 3,
			len(new_gene_translation_efficiency_values),
			len(new_gene_expression_factors))) - 1
		avg_new_gene_monomer_mass_fraction_heatmap = np.zeros(( 3,
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
		new_gene_mRNA_counts = [{} for id in new_gene_mRNA_ids]
		new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]
		new_gene_monomer_mass_fraction = [{} for id in new_gene_monomer_ids]
		new_gene_monomer_masses = [1 for id in new_gene_monomer_ids]
		for i in range(len(new_gene_monomer_ids)):
			new_gene_monomer_masses[i] = float(
				(sim_data.getter.get_mass(new_gene_monomer_ids[i]) / 1000
				 * 0.0000016605402) / (1 * g / mol))  # convert from g/mol to fg
		rnap_counts = {}
		ribosome_counts = {}
		ppgpp_counts = {}

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

				# Extract mRNA indexes for each new gene
				mRNA_counts_reader = TableReader(os.path.join(simOutDir,
															  'RNACounts'))
				mRNA_idx_dict = {rna[:-3]: i for i, rna in
								 enumerate(mRNA_counts_reader.readAttribute(
									 'mRNA_ids'))}
				new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id)
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

			avg_new_gene_mRNA_counts = read_stacked_columns(
				all_cells, 'RNACounts', 'mRNA_counts', fun=lambda
					x: np.mean( x[:, new_gene_mRNA_indexes], axis=0))
			avg_new_gene_monomer_counts = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts', fun=lambda x:
				np.mean( x[:,new_gene_monomer_indexes], axis=0))

			avg_new_gene_mRNA_counts = \
				avg_new_gene_mRNA_counts[exclude_timeout_cell_mask,]
			avg_new_gene_monomer_counts = \
				avg_new_gene_monomer_counts[exclude_timeout_cell_mask,]

			# Total protein mass
			avg_protein_mass = read_stacked_columns(
				all_cells, 'Mass', 'proteinMass', fun=lambda x: np.mean(x))
			avg_protein_mass = avg_protein_mass[exclude_timeout_cell_mask]

			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_counts[i][variant] = \
					np.log10(avg_new_gene_mRNA_counts[:, i] + 1)
				new_gene_monomer_counts[i][variant] = \
					np.log10(avg_new_gene_monomer_counts[:, i] + 1)
				new_gene_monomer_mass_fraction[i][variant] = \
					(avg_new_gene_monomer_counts.flatten() *
					 new_gene_monomer_masses[i]) / avg_protein_mass.flatten()

			# Ribosome counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')
				uniqueMoleculeCounts = TableReader(os.path.join(
					simOutDir, "UniqueMoleculeCounts"))
				ribosomeIndex = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')

			avg_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				fun=lambda x: np.mean(x[:, ribosomeIndex], axis=0))
			ribosome_counts[variant] = avg_ribosome_counts[
				exclude_timeout_cell_mask]

			# RNA polymerase counts
			rnapId = ["APORNAP-CPLX[c]"]
			(rnapCountsBulk,) = read_stacked_bulk_molecules(all_cells,
															(rnapId,))
			cell_id_vector = stacked_cell_identification(all_cells,
														 'Main',
														 'time')
			cell_ids, idx, cell_total_timesteps = np.unique(
				cell_id_vector, return_inverse=True, return_counts=True)
			sum_rnap_counts = np.bincount(idx, weights=rnapCountsBulk)
			avg_rnap_counts = sum_rnap_counts / cell_total_timesteps
			rnap_counts[variant] = avg_rnap_counts[exclude_timeout_cell_mask]

			# ppGpp counts
			avg_ppgpp_counts = read_stacked_columns(
				all_cells, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True, fun=lambda x: np.mean(x)).squeeze()
			ppgpp_counts[variant] = avg_ppgpp_counts[exclude_timeout_cell_mask]

			# Add values to heatmap data structures
			exp_index, trl_eff_index = variant_index_to_list_indices[variant]
			doubling_times_heatmap[0, trl_eff_index, exp_index] = round(
				np.mean(
				doubling_times[variant]))
			completed_gens_heatmap[0, trl_eff_index, exp_index] = \
				round(reached_count_gen[variant],2)
			i = 0 ### TODO: accomodate multiple new genes
			avg_new_gene_mRNA_counts_heatmap[0, trl_eff_index, exp_index] = \
				round(np.mean(new_gene_mRNA_counts[i][variant]), 2)
			avg_new_gene_monomer_counts_heatmap[0, trl_eff_index, exp_index]\
				= \
				round(np.mean(new_gene_monomer_counts[i][variant]), 2)
			avg_new_gene_monomer_mass_fraction_heatmap[0, trl_eff_index,\
														exp_index]	= \
				round(np.mean(new_gene_monomer_mass_fraction[i][variant]), 2)
			avg_ribosome_counts_heatmap[0, trl_eff_index, exp_index] = int(
				round(
				np.mean(ribosome_counts[variant])))
			avg_rnap_counts_heatmap[0, trl_eff_index, exp_index] = round(
				np.mean(rnap_counts[variant]))
			avg_ppgpp_counts_heatmap[0, trl_eff_index, exp_index] = round(
				np.mean(ppgpp_counts[variant]), 1)

			if exclude_early_gens == 1:
				# Add early gen values to the heatmap structure
				early_cell_mask = generations[variant] < MIN_LATE_CELL_INDEX
				if len(early_cell_mask) == 1:
					early_cell_mask = early_cell_mask[0]

				doubling_times_heatmap[1, trl_eff_index, exp_index] = round(
					np.mean(doubling_times[variant][early_cell_mask]))
				i = 0  ### TODO: accomodate multiple new genes
				avg_new_gene_mRNA_counts_heatmap[1, trl_eff_index, exp_index]\
					= \
					round(np.mean(new_gene_mRNA_counts[i][variant][early_cell_mask]), 2)
				avg_new_gene_monomer_counts_heatmap[
					1, trl_eff_index, exp_index] \
					= \
					round(np.mean(new_gene_monomer_counts[i][variant][early_cell_mask]), 2)

				avg_new_gene_monomer_mass_fraction_heatmap[1, trl_eff_index, \
														   exp_index] = \
					round(np.mean(new_gene_monomer_mass_fraction[i][variant][early_cell_mask]),
						  2)
				avg_ribosome_counts_heatmap[1, trl_eff_index, exp_index] = \
					int(round(np.mean(ribosome_counts[variant][early_cell_mask])))
				avg_rnap_counts_heatmap[1, trl_eff_index, exp_index] = round(
					np.mean(rnap_counts[variant][early_cell_mask]))
				avg_ppgpp_counts_heatmap[1, trl_eff_index, exp_index] = round(
					np.mean(ppgpp_counts[variant][early_cell_mask]), 1)

				# Add late gen values to the heatmap structure
				late_cell_mask = np.logical_and((generations[variant] >=
								  MIN_LATE_CELL_INDEX), \
								 (generations[variant] < MAX_CELL_INDEX))
				if len(late_cell_mask) == 1:
					late_cell_mask = late_cell_mask[0]
				if sum(late_cell_mask) != 0:
					doubling_times_heatmap[2, trl_eff_index, exp_index] = round(
						np.mean(doubling_times[variant][late_cell_mask]))
					i = 0  ### TODO: accomodate multiple new genes
					avg_new_gene_mRNA_counts_heatmap[2, trl_eff_index,exp_index] = \
						round(np.mean(new_gene_mRNA_counts[i][variant][late_cell_mask]), 2)
					avg_new_gene_monomer_counts_heatmap[2, trl_eff_index, exp_index] = \
						round(np.mean(new_gene_monomer_counts[i][variant][late_cell_mask]), 2)
					avg_new_gene_monomer_mass_fraction_heatmap[2, trl_eff_index, exp_index] = \
						round(np.mean(new_gene_monomer_mass_fraction[i][variant][late_cell_mask]),
							  2)
					avg_ribosome_counts_heatmap[2, trl_eff_index, exp_index] = \
						int(round(np.mean(ribosome_counts[variant][late_cell_mask])))
					avg_rnap_counts_heatmap[2, trl_eff_index, exp_index] = \
						round(np.mean(rnap_counts[variant][late_cell_mask]))
					avg_ppgpp_counts_heatmap[2, trl_eff_index, exp_index] = \
						round(np.mean(ppgpp_counts[variant][late_cell_mask]), 1)
		
		# Plotting
		print("---Plotting---")
		plot_descr = ["_all_gens"]
		if exclude_early_gens == 1:
			plot_descr += ["_early_gens", "_late_gens"]

		# Percent Completion
		fig, ax = plt.subplots(1, 1, figsize=(10, 5))
		self.heatmap(ax, variant_mask, completed_gens_heatmap[0, :, :],
					 completed_gens_heatmap[0, :, :],
					 "Expression Variant",
					 "Translation Efficiency Value (Normalized)",
					 new_gene_expression_factors,
					 new_gene_translation_efficiency_values,
					 "Percentage of Sims that Reached Generation " \
					 + str(COUNT_INDEX + 1))
		fig.tight_layout()
		plt.show()
		exportFigure(plt, plotOutDir, 'completed_gens_heatmap')

		for j in range(len(plot_descr)):

			# Doubling Time
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask, doubling_times_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Doubling Times")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'doubling_time_heatmap' +
						 plot_descr[j])

			# New Gene mRNA Counts
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask,
						 avg_new_gene_mRNA_counts_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Log(New Gene mRNA Counts+1)")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'new_gene_mRNA_heatmap' + plot_descr[j])

			# New Gene Monomer Counts
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask,
						 avg_new_gene_monomer_counts_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Log(New Gene Protein Counts+1)")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'new_gene_monomer_heatmap' + plot_descr[j])

			# New Gene Monomer Mass Fraction
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask,
						 avg_new_gene_monomer_mass_fraction_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "New Gene Monomer Mass Fraction")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir,
						 'new_gene_monomer_mass_fraction_heatmap' + plot_descr[j])

			# Ribosome Counts
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask, avg_ribosome_counts_heatmap[j, :,:],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "Ribosome Counts", "x-small")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'ribosome_heatmap' + plot_descr[j])

			# RNA Polymerase Counts
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask, avg_rnap_counts_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "RNA Polymerase Counts", "x-small")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'rnap_heatmap' + plot_descr[j])

			# ppGpp Concentration
			fig, ax = plt.subplots(1, 1, figsize=(10, 5))
			self.heatmap(ax, variant_mask, avg_ppgpp_counts_heatmap[j, :, :],
						 completed_gens_heatmap[0, :, :],
						 "Expression Variant",
						 "Translation Efficiency Value (Normalized)",
						 new_gene_expression_factors,
						 new_gene_translation_efficiency_values,
						 "ppGpp Concentration (uM)")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'ppgpp_heatmap' + plot_descr[j])

			plt.close('all')

if __name__ == "__main__":
	Plot().cli()
