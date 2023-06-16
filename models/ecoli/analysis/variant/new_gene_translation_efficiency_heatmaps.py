"""
Plot one value per index via heatmap for
new_gene_expression_and_translation_efficiency variant.

Possible Plots:
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
from wholecell.analysis.plotting_tools import heatmap
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
COUNT_INDEX = 5 #15 # TODO: CHANGE THIS BACK TO 16

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

MAX_CELL_LENGTH = 36000
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def save_new_gene_heatmap(self, heatmap_data, h, i, trl_eff_index,
							  exp_index, curr_heatmap_data, heatmap_details,
							  early_cell_mask, late_cell_mask):
		heatmap_data[h][i][0, trl_eff_index, exp_index] = round(
			np.mean(curr_heatmap_data), heatmap_details[h][
				'num_digits_rounding'])

		if exclude_early_gens:
			heatmap_data[h][i][1, trl_eff_index, exp_index] = \
				round(
					np.mean(curr_heatmap_data[early_cell_mask]),
					heatmap_details[h]['num_digits_rounding'])
			if sum(late_cell_mask) != 0:
				heatmap_data[h][i][2, trl_eff_index, exp_index] = \
					round(
						np.mean(curr_heatmap_data[late_cell_mask]),
						heatmap_details[h]['num_digits_rounding'])

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		print("Running analysis script with exclude_timeout_cells=",
			  exclude_timeout_cells, " and exclude_early_gens=",
			  exclude_early_gens)

		# Specify which subset of heatmaps should be made
		# Completed_gens heatmap is always made, because it is used to
		# create the other heatmaps, and should not be included here
		heatmaps_to_make = {"doubling_times_heatmap", "cell_volume_heatmap",
							"cell_mass_heatmap", "cell_dry_mass_heatmap",
							"cell_mRNA_mass_heatmap",
							"cell_protein_mass_heatmap",
							"ppgpp_concentration_heatmap",
							"rnap_counts_heatmap",
							"ribosome_counts_heatmap",
							#"rnap_crowding_heatmap", # TODO This one is definitely slower
							#"ribosome_crowding_heatmap" # TODO This one is definitely slower
							"new_gene_mRNA_counts_heatmap",
							"new_gene_mRNA_NTP_fraction_heatmap",
							"new_gene_monomer_counts_heatmap",
							"new_gene_mRNA_mass_fraction_heatmap",
							"new_gene_monomer_mass_fraction_heatmap",
							"new_gene_rnap_init_rate_heatmap",
							"new_gene_ribosome_init_rate_heatmap",
							"new_gene_rnap_time_overcrowded_heatmap",
							"new_gene_ribosome_time_overcrowded_heatmap"
							}

		# Placeholders for lambda functions
		ribosome_index = -1
		rnap_id = ["APORNAP-CPLX[c]"]
		new_gene_mRNA_indexes = []
		new_gene_monomer_indexes = []

		# Details needed to create all possible heatmaps
		# TODO: make sure all the heatmaps have units
		heatmap_details = {
			"doubling_times_heatmap" :
				{'data_table': 'Main',
				 'data_column': 'time',
				 'remove_first': False,
				 'function_to_apply': lambda x: (x[-1] - x[0]) / 60.,
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'medium',
				 'plot_title': 'Doubling Time (minutes)',
				 },
			"cell_volume_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'cellVolume',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x),
				 'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'Cell Volume (fL)',
				 },
			"cell_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'cellMass',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x),
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Cell Mass (fg)',
				 },
			"cell_dry_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'dryMass',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x),
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Dry Cell Mass (fg)',
				 },
			"cell_mRNA_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'mRnaMass',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x),
				 'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'Total mRNA Mass (fg)',
				 },
			"cell_protein_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'proteinMass',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x),
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Total Protein Mass (fg)',
				 },
			"ppgpp_concentration_heatmap":
				{'data_table': 'GrowthLimits',
				 'data_column': 'ppgpp_conc',
				 'remove_first': True,
				 'function_to_apply': lambda x: np.mean(x),
				 'default_value': -1,
				 'num_digits_rounding': 1,
				 'box_text_size': 'medium',
				 'plot_title': 'ppGpp Concentration (uM)',
				 },
			"rnap_counts_heatmap":
				{'data_table': 'N/A',
				 'data_column': 'N/A',
				 'remove_first': False,
				 'function_to_apply': lambda x: 0,
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNA Polymerase Counts',
				 },
			"ribosome_counts_heatmap":
				{'data_table': 'UniqueMoleculeCounts',
				 'data_column': 'uniqueMoleculeCounts',
				 'remove_first': False,
				 'function_to_apply':
					 lambda x: np.mean(x[:, ribosome_index], axis=0),
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Ribosome Counts',
				 },
			"rnap_crowding_heatmap":
				{'data_table': 'RnaSynthProb',
				 'data_column': 'N/A',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x, axis = 0),
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'medium',
				 'plot_title': 'RNA Polymerase Crowding: # of TUs',
				 },
			"ribosome_crowding_heatmap":
				{'data_table': 'RibosomeData',
				 'data_column': 'N/A',
				 'remove_first': False,
				 'function_to_apply': lambda x: np.mean(x, axis = 0),
				 'default_value': -1,
				 'num_digits_rounding': 0,
				 'box_text_size': 'medium',
				 'plot_title': 'Ribosome Crowding: # of Monomers',
				 },
			"new_gene_mRNA_counts_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'Log(New Gene mRNA Counts+1)',
				 },
			"new_gene_monomer_counts_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'Log(New Gene Protein Counts+1)',
				 },
			"new_gene_mRNA_mass_fraction_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'New Gene mRNA Mass Fraction',
				 },
			"new_gene_mRNA_NTP_fraction_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 4,
				 'box_text_size': 'small',
				 'plot_title': 'New Gene mRNA Mass Fraction',
				 },
			"new_gene_monomer_mass_fraction_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'New Gene Protein Mass Fraction',
				 },
			"new_gene_rnap_init_rate_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'New Gene RNAP Initialization Rate',
				 },
			"new_gene_ribosome_init_rate_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'New Gene Ribosome Initalization Rate',
				 },
			"new_gene_rnap_time_overcrowded_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'Fraction of Time RNA Polymerase Overcrowded '
							   'New Gene',
				 },
			"new_gene_ribosome_time_overcrowded_heatmap":
				{'default_value': -1,
				 'num_digits_rounding': 2,
				 'box_text_size': 'medium',
				 'plot_title': 'Fraction of Time Ribosome Overcrowded New '
							   'Gene',
				 },
		}

		assert "completed_gens_heatmap" not in heatmaps_to_make, \
			"the completed_gens_heatmap is run by default, do not include in heatmaps_to_make"
		for h in heatmaps_to_make:
			assert h in heatmap_details, "Heatmap " + h + " is not an option"

		new_gene_heatmaps = {h for h in heatmaps_to_make
								 if h.startswith("new_gene_")}
		special_new_gene_heatmaps = {"new_gene_mRNA_NTP_fraction_heatmap"}
		standard_new_gene_heatmaps = new_gene_heatmaps - special_new_gene_heatmaps
		non_new_gene_heatmaps = heatmaps_to_make - new_gene_heatmaps
		special_non_new_gene_heatmaps = {"rnap_counts_heatmap",
										 "rnap_crowding_heatmap",
										 "ribosome_crowding_heatmap"}
		standard_non_new_gene_heatmaps = non_new_gene_heatmaps - \
										 special_non_new_gene_heatmaps
		# Dictionary for storing the data (numpy arrays) for each heatmap
		heatmap_data = {}


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


		if "new_gene_mRNA_NTP_fraction_heatmap" in heatmaps_to_make:
			# Determine NTP ids
			with open(simDataFile, 'rb') as f:
				sim_data = pickle.load(f)
			ntp_ids = list(sim_data.ntp_code_to_id_ordered.values())

			# Determine number of NTPs per new gene mRNA
			new_gene_mRNA_ntp_counts = [{} for id in new_gene_mRNA_ids]
			all_rna_counts_ACGU = sim_data.process.transcription.rna_data[
				'counts_ACGU'].asNumber()
			rna_ids = sim_data.process.transcription.rna_data['id']
			rna_id_to_index_mapping = {rna[:-3]: i for i, rna in
									   enumerate(rna_ids)}
			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_index = rna_id_to_index_mapping[
					new_gene_mRNA_ids[i]]
				for ntp_index in range(len(ntp_ids)):
					new_gene_mRNA_ntp_counts[i][ntp_ids[ntp_index]] = \
						all_rna_counts_ACGU[new_gene_mRNA_index, ntp_index]

			all_mRNA_counts_ACGU = \
				all_rna_counts_ACGU[sim_data.process.transcription.rna_data[
					"is_mRNA"]]


		# Create data structures that to use for the heatmaps
		heatmap_data["completed_gens_heatmap"] = np.zeros((
					1, len(new_gene_translation_efficiency_values),
					len(new_gene_expression_factors)))
		for h in non_new_gene_heatmaps:
			heatmap_data[h] = np.zeros((
				3, len(new_gene_translation_efficiency_values),
				len(new_gene_expression_factors))) + heatmap_details[h][
				'default_value']
		for h in standard_new_gene_heatmaps:
			heatmap_data[h] = [np.zeros((
				3, len(new_gene_translation_efficiency_values),
				len(new_gene_expression_factors))) + heatmap_details[h][
				'default_value'] for x in range(len(new_gene_mRNA_ids))]
		if "new_gene_mRNA_NTP_fraction_heatmap" in heatmaps_to_make:
			heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"] = [{} for x
				in range(len(new_gene_mRNA_ids))]
			for i in range(len(new_gene_mRNA_ids)):
				for ntp_id in ntp_ids:
					heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"][i][
						ntp_id] = np.zeros((3,
						len(new_gene_translation_efficiency_values),
						len(new_gene_expression_factors))) + heatmap_details[
						"new_gene_mRNA_NTP_fraction_heatmap"]['default_value']


		# Data extraction
		print("---Data Extraction---")
		reached_count_gen = {}
		generations = {}
		new_gene_mRNA_counts = [{} for id in new_gene_mRNA_ids]
		new_gene_monomer_counts = [{} for id in new_gene_monomer_ids]
		new_gene_monomer_masses = [1 for id in new_gene_monomer_ids]
		for i in range(len(new_gene_monomer_ids)):
			new_gene_monomer_masses[i] = float(
				(sim_data.getter.get_mass(new_gene_monomer_ids[i]) / 1000
				 * 0.0000016605402) / (1 * g / mol))  # convert from g/mol to fg
		new_gene_mRNA_masses = [1 for id in new_gene_mRNA_ids]
		for i in range(len(new_gene_mRNA_ids)):
			new_gene_mRNA_masses[i] = float(
				(sim_data.getter.get_mass(new_gene_mRNA_ids[i]) / 1000
				 * 0.0000016605402) / (1 * g / mol))  # convert from g/mol to fg
		all_mRNA_ntp_totals = {}  # {variant: {NTP id: values}}

		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant],
										  only_successful=True)
			exp_index, trl_eff_index = variant_index_to_list_indices[variant]

			if len(all_cells) == 0:
				continue

			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', MAX_CELL_LENGTH,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			all_cells_gens = np.array([int(os.path.basename(os.path.dirname(
				cell_path))[-6:]) for cell_path in all_cells])[exclude_timeout_cell_mask]
			generations[variant] = all_cells_gens

			if exclude_early_gens == 1:
				# Add early gen values to the heatmap structure
				early_cell_mask = generations[variant] < MIN_LATE_CELL_INDEX
				if len(early_cell_mask) == 1:
					early_cell_mask = early_cell_mask[0]
				# Add late gen values to the heatmap structure
				late_cell_mask = np.logical_and((generations[variant] >=
									MIN_LATE_CELL_INDEX), \
									(generations[variant] < MAX_CELL_INDEX))
				if len(late_cell_mask) == 1:
					late_cell_mask = late_cell_mask[0]

			# Completed Gens Heatmap: Count the number of simulations that
			# reach gen COUNT_INDEX + 1
			num_count_gen = len(self.ap.get_cells(variant=[variant],
												  generation=[COUNT_INDEX],
												  only_successful=True))
			num_zero_gen = len(self.ap.get_cells(variant=[variant],
												 generation=[0],
												 only_successful=True))
			reached_count_gen[variant] = num_count_gen / num_zero_gen
			heatmap_data["completed_gens_heatmap"][0, trl_eff_index,
				exp_index]	= round(reached_count_gen[variant], 2)

			# Get ribosome index for ribosome counts heatmap (if needed)
			if variant == min_variant and "ribosome_counts_heatmap" in heatmaps_to_make:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')
				uniqueMoleculeCounts = TableReader(os.path.join(
					simOutDir, "UniqueMoleculeCounts"))
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')

			# Save the data for heatmaps that are standard and don't involve
			# new genes
			for h in standard_non_new_gene_heatmaps:
				curr_heatmap_data = read_stacked_columns(all_cells,
					heatmap_details[h]['data_table'],
					heatmap_details[h]['data_column'],
					remove_first=heatmap_details[h]['remove_first'],
					fun=heatmap_details[h]['function_to_apply'])[
					exclude_timeout_cell_mask]
				heatmap_data[h][0, trl_eff_index, exp_index] = round(
					np.mean(curr_heatmap_data), heatmap_details[h][
						'num_digits_rounding'])
				if exclude_early_gens:
					heatmap_data[h][1, trl_eff_index, exp_index] = round(
						np.mean(curr_heatmap_data[early_cell_mask]),
						heatmap_details[h]['num_digits_rounding'])
					if sum(late_cell_mask) != 0:
						heatmap_data[h][2, trl_eff_index, exp_index] = round(
							np.mean(curr_heatmap_data[late_cell_mask]),
							heatmap_details[h]['num_digits_rounding'])

			# Save the data for heatmaps that don't involve new genes but
			# involve nonstandard data loading/manipulation

			# RNAP Counts
			if 'rnap_counts_heatmap' in heatmaps_to_make:
				(rnapCountsBulk,) = read_stacked_bulk_molecules(
					all_cells,(rnap_id,))
				cell_id_vector = stacked_cell_identification(all_cells, 'Main',
					'time')
				cell_ids, idx, cell_total_timesteps = np.unique(
					cell_id_vector, return_inverse=True, return_counts=True)
				sum_rnap_counts = np.bincount(idx, weights=rnapCountsBulk)
				avg_rnap_counts = (sum_rnap_counts / cell_total_timesteps)[
					exclude_timeout_cell_mask]

				heatmap_data["rnap_counts_heatmap"][0, trl_eff_index, exp_index] = round(
					np.mean(avg_rnap_counts), heatmap_details[
						"rnap_counts_heatmap"][
						'num_digits_rounding'])
				if exclude_early_gens:
					heatmap_data["rnap_counts_heatmap"][1, trl_eff_index, exp_index] = round(
						np.mean(avg_rnap_counts[early_cell_mask]),
						heatmap_details["rnap_counts_heatmap"]['num_digits_rounding'])
					if sum(late_cell_mask) != 0:
						heatmap_data["rnap_counts_heatmap"][2, trl_eff_index, exp_index] = round(
							np.mean(avg_rnap_counts[late_cell_mask]),
							heatmap_details["rnap_counts_heatmap"]['num_digits_rounding'])

			# RNA Polymerase Overcrowding
			if "rnap_crowding_heatmap" in heatmaps_to_make:
				avg_actual_rna_synth_prob = read_stacked_columns(
					all_cells,
					'RnaSynthProb', 'actual_rna_synth_prob',
					fun=lambda x: np.mean(x, axis=0))
				avg_target_rna_synth_prob = read_stacked_columns(
					all_cells,
					'RnaSynthProb', 'target_rna_synth_prob',
					fun=lambda x: np.mean(x, axis=0))

				# Get indexes of transcription units that on
				# average were overcrowded in any generation for any seed
				avg_overcrowded_tu_indexes = np.where(
					sum(avg_actual_rna_synth_prob <
						avg_target_rna_synth_prob) > 0)[0]
				n_overcrowded_tus = len(avg_overcrowded_tu_indexes)

				heatmap_data["rnap_crowding_heatmap"][0, trl_eff_index,
					exp_index] = n_overcrowded_tus
				if exclude_early_gens == 1:
					heatmap_data["rnap_crowding_heatmap"][1, trl_eff_index, exp_index] = \
						len(np.where(sum((avg_actual_rna_synth_prob <
						avg_target_rna_synth_prob)[early_cell_mask, :]) > 0)[0])
					if sum(late_cell_mask) != 0:
						heatmap_data["rnap_crowding_heatmap"][2, trl_eff_index, exp_index] = \
							len(np.where(sum((avg_actual_rna_synth_prob <
							avg_target_rna_synth_prob)[late_cell_mask, :]) > 0)[0])

			# Ribosome Overcrowding
			if "ribosome_crowding_heatmap" in heatmaps_to_make:
				avg_actual_prob_translation_per_transcript = read_stacked_columns(
					all_cells,
					'RibosomeData', 'actual_prob_translation_per_transcript',
					fun=lambda x: np.mean(x, axis=0))
				avg_target_prob_translation_per_transcript = read_stacked_columns(
					all_cells,
					'RibosomeData', 'target_prob_translation_per_transcript',
					fun=lambda x: np.mean(x, axis=0))

				# Get indexes of proteins corresponding to mRNAs that on
				# average were overcrowded in any generation for any seed
				avg_overcrowded_monomer_indexes = np.where(
					sum(avg_actual_prob_translation_per_transcript <
						avg_target_prob_translation_per_transcript) > 0)[0]
				n_overcrowded_monomers = len(avg_overcrowded_monomer_indexes)

				heatmap_data["ribosome_crowding_heatmap"][0, trl_eff_index,
					exp_index]	= n_overcrowded_monomers
				if exclude_early_gens == 1:
					heatmap_data["ribosome_crowding_heatmap"][1, trl_eff_index, exp_index] = \
						len(np.where(sum((avg_actual_prob_translation_per_transcript <
						avg_target_prob_translation_per_transcript)[
						early_cell_mask, :]) > 0)[0])
					if sum(late_cell_mask) != 0:
						heatmap_data["ribosome_crowding_heatmap"][2, trl_eff_index, exp_index] = \
							len(np.where(sum((avg_actual_prob_translation_per_transcript <
							avg_target_prob_translation_per_transcript)[
							late_cell_mask, :]) > 0)[0])


			# New Gene Based Heatmaps
			if not new_gene_heatmaps:
				continue

			# New Gene mRNA and Monomer Counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract cistron indexes for each new gene
				rnap_reader = TableReader(os.path.join(simOutDir,
													   'RnapData'))
				cistron_idx_dict = {cis: i for i, cis in
									enumerate(rnap_reader.readAttribute(
										'cistron_ids'))}
				new_gene_cistron_indexes = [cistron_idx_dict.get(mRNA_id)
											for mRNA_id in new_gene_mRNA_ids]

				# Extract RNA indexes for each new gene
				rnap_reader = TableReader(os.path.join(simOutDir,
													   'RnaSynthProb'))
				RNA_idx_dict = {rna[:-3]: i for i, rna in
								enumerate(rnap_reader.readAttribute(
									'rnaIds'))}
				new_gene_RNA_indexes = [RNA_idx_dict.get(mRNA_id)
										for mRNA_id in new_gene_mRNA_ids]

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

			avg_new_gene_mRNA_counts = (read_stacked_columns(
				all_cells, 'RNACounts', 'mRNA_counts', fun=lambda
					x: np.mean( x[:, new_gene_mRNA_indexes], axis=0)))[exclude_timeout_cell_mask,]
			avg_new_gene_monomer_counts = (read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts', fun=lambda x:
				np.mean( x[:,new_gene_monomer_indexes], axis=0)))[exclude_timeout_cell_mask,]

			# New Gene mRNA Counts
			if "new_gene_mRNA_counts_heatmap" in heatmaps_to_make:
				for i in range(len(new_gene_mRNA_ids)):
					new_gene_mRNA_counts[i][variant] = \
						np.log10(avg_new_gene_mRNA_counts[:, i] + 1)
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_mRNA_counts_heatmap", i, trl_eff_index,
						exp_index, new_gene_mRNA_counts[i][variant],
						heatmap_details, early_cell_mask, late_cell_mask);

			# New Gene Protein Counts
			if "new_gene_monomer_counts_heatmap" in heatmaps_to_make:
				for i in range(len(new_gene_mRNA_ids)):
					new_gene_monomer_counts[i][variant] = \
						np.log10(avg_new_gene_monomer_counts[:, i] + 1)
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_monomer_counts_heatmap", i,
						trl_eff_index, exp_index, new_gene_monomer_counts[i][variant],
						heatmap_details, early_cell_mask, late_cell_mask)

			# New Gene mRNA Mass Fraction
			if "new_gene_mRNA_mass_fraction_heatmap":
				# Total mRNA mass
				avg_mRNA_mass = read_stacked_columns(
					all_cells, 'Mass', 'mRnaMass', fun=lambda x: np.mean(x))
				avg_mRNA_mass = avg_mRNA_mass[exclude_timeout_cell_mask]
				for i in range(len(new_gene_mRNA_ids)):
					new_gene_mRNA_mass_fraction = \
						(avg_new_gene_mRNA_counts[:,i] *
						new_gene_mRNA_masses[i]) / avg_mRNA_mass[:,i]
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_mRNA_mass_fraction_heatmap",
						i, trl_eff_index, exp_index,
						new_gene_mRNA_mass_fraction, heatmap_details,
						early_cell_mask, late_cell_mask)

			# New Gene mRNA NTP Fraction
			if "new_gene_mRNA_NTP_fraction_heatmap" in heatmaps_to_make:
				# Total NTPs in all mRNAs
				avg_mRNA_counts = read_stacked_columns(
					all_cells, 'RNACounts', 'mRNA_counts', fun=lambda
						x: np.mean(x, axis=0))
				all_mRNA_ntp_totals[variant] = {}

				for ntp_index in range(len(ntp_ids)):
					ntp_id = ntp_ids[ntp_index]
					all_mRNA_ntp_totals[variant][ntp_id] = \
						(avg_mRNA_counts @ all_mRNA_counts_ACGU[:, ntp_index])[
							exclude_timeout_cell_mask]
					for i in range(len(new_gene_mRNA_ids)):
						heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"][i][
							ntp_id][0,trl_eff_index, exp_index] = round(np.mean(
							(avg_new_gene_mRNA_counts[:,i] *
							 new_gene_mRNA_ntp_counts[i][ntp_id]) /
							all_mRNA_ntp_totals[variant][ntp_id]),
							heatmap_details["new_gene_mRNA_NTP_fraction_heatmap"][
								'num_digits_rounding'])
						if exclude_early_gens == 1:
							heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"][
								i][ntp_id][1, trl_eff_index, exp_index] = \
								round(np.mean(
								(avg_new_gene_mRNA_counts[:,i][early_cell_mask] *
								new_gene_mRNA_ntp_counts[i][ntp_id]) /
								all_mRNA_ntp_totals[variant][ntp_id][early_cell_mask]),
								heatmap_details["new_gene_mRNA_NTP_fraction_heatmap"][
								'num_digits_rounding'])

							if sum(late_cell_mask) != 0:
								heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"][
									i][ntp_id][2, trl_eff_index, exp_index] = \
									round(np.mean(
									(avg_new_gene_mRNA_counts[:,i][late_cell_mask] *
									new_gene_mRNA_ntp_counts[i][ntp_id]) /
									all_mRNA_ntp_totals[variant][ntp_id][late_cell_mask]),
									heatmap_details["new_gene_mRNA_NTP_fraction_heatmap"][
									'num_digits_rounding'])

			# New Gene Protein Mass Fraction
			if "new_gene_monomer_mass_fraction_heatmap" in heatmaps_to_make:
				avg_protein_mass = (read_stacked_columns(
					all_cells, 'Mass', 'proteinMass',
					fun=lambda x: np.mean(x)))[exclude_timeout_cell_mask,]
				for i in range(len(new_gene_mRNA_ids)):
					new_gene_monomer_mass_fraction = \
						(avg_new_gene_monomer_counts[:,i] *
						 new_gene_monomer_masses[i]) / avg_protein_mass[:,i]
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_monomer_mass_fraction_heatmap", i,
						trl_eff_index, exp_index,
						new_gene_monomer_mass_fraction,
						heatmap_details, early_cell_mask, late_cell_mask)

			# New Gene RNAP Initialization Rate
			if "new_gene_rnap_init_rate_heatmap" in heatmaps_to_make:
				avg_new_gene_copy_number = (read_stacked_columns(
					all_cells, 'RnaSynthProb', 'gene_copy_number',
					fun=lambda x: np.mean(x[:, new_gene_cistron_indexes],
					axis=0)))[exclude_timeout_cell_mask,]
				avg_new_gene_rnap_init_rates = (read_stacked_columns(
					all_cells, 'RnapData', 'rna_init_event_per_cistron', fun=lambda
						x: np.mean(x[:, new_gene_cistron_indexes],
						axis=0)))[exclude_timeout_cell_mask,] / avg_new_gene_copy_number
				for i in range(len(new_gene_mRNA_ids)):
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_rnap_init_rate_heatmap", i, trl_eff_index,
						exp_index, avg_new_gene_rnap_init_rates[:,i],
						heatmap_details, early_cell_mask, late_cell_mask)

			# New Gene Ribosome Initialization Rate
			if "new_gene_ribosome_init_rate_heatmap" in heatmaps_to_make:
				avg_new_gene_ribosome_init_rates = (read_stacked_columns(
					all_cells, 'RibosomeData', 'ribosome_init_event_per_monomer',
					fun=lambda x: np.mean(x[:, new_gene_monomer_indexes],
					axis=0)))[exclude_timeout_cell_mask,] / \
					avg_new_gene_mRNA_counts
				for i in range(len(new_gene_mRNA_ids)):
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_ribosome_init_rate_heatmap",i, trl_eff_index,
						exp_index, avg_new_gene_ribosome_init_rates[:,i],
						heatmap_details, early_cell_mask, late_cell_mask)

			# Fraction of Time RNAP Overcrowds New Gene
			if "new_gene_rnap_time_overcrowded_heatmap" in heatmaps_to_make:
				# Average fraction of time steps that RNA polymerase
				# overcrowding occurs for new genes per generation
				new_gene_num_time_steps_rnap_overcrowded = (
					read_stacked_columns(all_cells, 'RnaSynthProb',
					'tu_is_overcrowded',
					fun=lambda x: np.sum(x[:, new_gene_RNA_indexes], axis=0)
					/ (x[:, new_gene_RNA_indexes].shape[0]))
					)[exclude_timeout_cell_mask,]
				for i in range(len(new_gene_mRNA_ids)):
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_rnap_time_overcrowded_heatmap", i,
						trl_eff_index, exp_index,
						new_gene_num_time_steps_rnap_overcrowded[:,i],
						heatmap_details, early_cell_mask, late_cell_mask)

			# Fraction of Time Ribosome Overcrowds New Gene
			if "new_gene_ribosome_time_overcrowded_heatmap" in \
					heatmaps_to_make:
				# Average fraction of time steps that ribosome
				# overcrowding occurs for new genes per generation
				new_gene_num_time_steps_ribosome_overcrowded = (
					read_stacked_columns(all_cells, 'RibosomeData', 'mRNA_is_overcrowded',
					fun=lambda x: np.sum(x[:, new_gene_monomer_indexes], axis=0) / (
					x[:, new_gene_monomer_indexes].shape[0]))
					)[exclude_timeout_cell_mask,]
				for i in range(len(new_gene_mRNA_ids)):
					self.save_new_gene_heatmap(heatmap_data,
						"new_gene_ribosome_time_overcrowded_heatmap", i,
						trl_eff_index, exp_index,
						new_gene_num_time_steps_ribosome_overcrowded[:,i],
						heatmap_details, early_cell_mask, late_cell_mask)

		# Plotting
		print("---Plotting---")
		plot_descr = ["_all_gens"]
		if exclude_early_gens == 1:
			plot_descr += ["_early_gens", "_late_gens"]
		heatmap_x_label = "Expression Variant"
		heatmap_y_label = "Translation Efficiency Value"

		# TODO: Change these? Make them dynamic?
		figsize_x = 10
		figsize_y = 5

		# Plot percent completion heatmap
		fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
		heatmap(self, ax, variant_mask,
				heatmap_data["completed_gens_heatmap"][0,:,:],
				heatmap_data["completed_gens_heatmap"][0,:,:],
				heatmap_x_label,
				heatmap_y_label,
				new_gene_expression_factors,
				new_gene_translation_efficiency_values,
				"Percentage of Sims that Reached Generation " \
				+ str(COUNT_INDEX + 1))
		fig.tight_layout()
		plt.show()
		exportFigure(plt, plotOutDir, 'completed_gens_heatmap')

		# Plot heatmaps that do not involve new genes
		for h in non_new_gene_heatmaps:
			for j in range(len(plot_descr)):
				fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
				heatmap(self, ax, variant_mask,
						heatmap_data[h][j,:,:],
						heatmap_data["completed_gens_heatmap"][0,:,:],
						heatmap_x_label,
						heatmap_y_label,
						new_gene_expression_factors,
						new_gene_translation_efficiency_values,
						heatmap_details[h]['plot_title'],
						heatmap_details[h]['box_text_size'])
				fig.tight_layout()
				plt.show()
				exportFigure(plt, plotOutDir, h + plot_descr[j])

			plt.close("all")

		# Plot standard new gene heatmaps
		for h in standard_new_gene_heatmaps:
			for i in range(len(new_gene_mRNA_ids)):
				for j in range(len(plot_descr)):
					fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
					heatmap(self, ax, variant_mask,
							heatmap_data[h][i][j,:,:],
							heatmap_data["completed_gens_heatmap"][0,:,:],
							heatmap_x_label,
							heatmap_y_label,
							new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							heatmap_details[h]['plot_title'] + ": " + new_gene_mRNA_ids[i][:-4],
							heatmap_details[h]['box_text_size'])
					fig.tight_layout()
					plt.show()
					exportFigure(plt, plotOutDir, h + "_" +
						new_gene_mRNA_ids[i][:-4] + plot_descr[j])

			plt.close('all')

		# New Gene mRNA NTP Fraction
		if "new_gene_mRNA_NTP_fraction_heatmap" in heatmaps_to_make:
			for i in range(len(new_gene_mRNA_ids)):
				for j in range(len(plot_descr)):
					for ntp_id in ntp_ids:
						fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
						heatmap(self, ax, variant_mask,
							heatmap_data["new_gene_mRNA_NTP_fraction_heatmap"][
								i][ntp_id][j, :, :],
							heatmap_data["completed_gens_heatmap"][0, :, :],
							heatmap_x_label,
							heatmap_y_label,
							new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							heatmap_details["new_gene_mRNA_NTP_fraction_heatmap"][
								'plot_title'] + ": " + new_gene_mRNA_ids[i][:-4],
							heatmap_details["new_gene_mRNA_NTP_fraction_heatmap"][
								'box_text_size'])
						fig.tight_layout()
						plt.show()
						exportFigure(plt, plotOutDir,
							'new_gene_mRNA_' + ntp_id[:-3] + '_fraction_heatmap'
							+ "_" + new_gene_mRNA_ids[i][:-4] + plot_descr[j])

		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
