"""
Heatmap dashboard for the new_gene_trl_eff_sweep variant.

Adapted from new_gene_kcat_translation_efficiency_heatmaps.py. Heatmaps are
2 rows (no-kcat, kcat) x (N_TRL_EFF + 1) columns (control + trl_eff values).
Only 1 expression factor (8.5), so the expression axis collapses into the
kcat row dimension.

Uses the same HEATMAPS_TO_MAKE_LIST and data extraction methods from the
parent class via inheritance.
"""

import numpy as np
from matplotlib import pyplot as plt
import math
import os.path
import pickle

from models.ecoli.analysis.variant.new_gene_kcat_translation_efficiency_heatmaps import (
	Plot as KcatHeatmapsPlot,
	HEATMAPS_TO_MAKE_LIST,
	DASHBOARD_FLAG,
	STD_DEV_FLAG,
	COUNT_INDEX,
	MIN_CELL_INDEX,
	MAX_CELL_INDEX,
	FONT_SIZE,
	capacity_gene_monomer_id,
	capacity_gene_common_name,
)
from models.ecoli.sim.variants.new_gene_internal_shift import (
	determine_new_gene_ids_and_indices,
)
from models.ecoli.sim.variants.new_gene_trl_eff_sweep import (
	TRL_EFF_VALUES,
	EXPRESSION_FACTOR,
	N_TRL_EFF,
	KCAT_HALF_START,
)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.plotting_tools import heatmap


class Plot(KcatHeatmapsPlot):
	"""
	Subclass that overrides do_plot to handle the new_gene_trl_eff_sweep
	variant index layout (2 rows: no-kcat/kcat, N_TRL_EFF+1 columns:
	control + trl_eff values).

	All data extraction methods are inherited from the parent class.
	"""

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			self.sim_data = pickle.load(f)

		# Determine new gene mRNA and monomer ids
		(self.new_gene_mRNA_ids, self.new_gene_indices,
			self.new_gene_monomer_ids,
			self.new_gene_monomer_indices) = determine_new_gene_ids_and_indices(
			self.sim_data)

		# Axis labels for heatmap:
		# - "columns" (exp axis) = 2 values: no-kcat, kcat
		# - "rows" (trl_eff axis) = N_TRL_EFF + 1 values: control + trl_eff values
		n_cols = 2  # [no-kcat, kcat]
		n_rows = N_TRL_EFF + 1  # control + 21 trl_eff values

		# Column labels (expression axis in parent code)
		col_labels = ['No kcat', '+kcat']
		# Row labels (trl_eff axis in parent code)
		row_labels = ['Control (KO)'] + [str(v) for v in TRL_EFF_VALUES]

		# Map variant indexes to (trl_eff_index, exp_index) for heatmap data.
		# trl_eff_index = local_index (0=control, 1..21=trl_eff values)
		# exp_index = 0 (no-kcat) or 1 (kcat)
		variants = self.ap.get_variants()
		variant_index_to_list_indices = {}
		variant_mask = np.zeros((n_rows, n_cols), dtype=bool)

		for index in variants:
			is_kcat = (index >= KCAT_HALF_START)
			local_index = index - KCAT_HALF_START if is_kcat else index
			exp_index = 1 if is_kcat else 0
			trl_eff_index = local_index  # 0=control, 1..21=trl_eff

			variant_index_to_list_indices[index] = np.array([
				exp_index, trl_eff_index])
			variant_mask[trl_eff_index, exp_index] = True

		# Set up heatmap details (inherited from parent)
		# Re-execute the heatmap_details setup from parent
		heatmaps_to_make = set(HEATMAPS_TO_MAKE_LIST)
		default_is_nonstandard_data_retrieval = False
		default_is_nonstandard_plot = False
		default_value = -1
		default_remove_first = False
		default_function_to_apply = lambda x: np.mean(x)
		default_num_digits_rounding = 2
		default_box_text_size = 'medium'

		self.heatmap_details = {
			"doubling_times_heatmap":
				{'data_table': 'Main',
				 'data_column': 'time',
				 'function_to_apply': lambda x: (x[-1] - x[0]) / 60.,
				 'num_digits_rounding': 0,
				 'plot_title': 'Doubling Time (minutes)',
				},
			"cell_volume_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'cellVolume',
				 'plot_title': 'Cell Volume (fL)',
				},
			"cell_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'cellMass',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Cell Mass (fg)',
				},
			"cell_dry_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'dryMass',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Dry Cell Mass (fg)',
				},
			"cell_mRNA_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'mRnaMass',
				 'plot_title': 'Total mRNA Mass (fg)',
				},
			"cell_protein_mass_heatmap":
				{'data_table': 'Mass',
				 'data_column': 'proteinMass',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Total Protein Mass (fg)',
				},
			"ppgpp_concentration_heatmap":
				{'data_table': 'GrowthLimits',
				 'data_column': 'ppgpp_conc',
				 'remove_first': True,
				 'num_digits_rounding': 1,
				 'plot_title': 'ppGpp Concentration (uM)',
				},
			"ribosome_init_events_heatmap":
				{'data_table': 'RibosomeData',
				 'data_column': 'didInitialize',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Ribosome Init Events Per Time Step',
				},
			"rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNA Polymerase (RNAP) Counts',
				},
			"ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'UniqueMoleculeCounts',
				 'data_column': 'uniqueMoleculeCounts',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Ribosome Counts',
				},
			"active_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Active RNA Polymerase (RNAP) Counts',
				},
			"active_ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'UniqueMoleculeCounts',
				 'data_column': 'uniqueMoleculeCounts',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Active Ribosome Counts',
				},
			"free_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Free RNA Polymerase (RNAP) Counts',
				},
			"free_ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'UniqueMoleculeCounts',
				 'data_column': 'uniqueMoleculeCounts',
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Free Ribosome Counts',
				},
			"rnap_ribosome_counts_ratio_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 4,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNAP Counts / Ribosome Counts',
				},
			"rnap_crowding_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'RnaSynthProb',
				 'function_to_apply': lambda x: np.mean(x, axis=0),
				 'num_digits_rounding': 0,
				 'plot_title': 'RNAP Crowding: # of TUs',
				},
			"ribosome_crowding_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'data_table': 'RibosomeData',
				 'function_to_apply': lambda x: np.mean(x, axis=0),
				 'num_digits_rounding': 0,
				 'plot_title': 'Ribosome Crowding: # of Monomers',
				},
			"weighted_avg_translation_efficiency_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Weighted Avg Translation Efficiency',
				 'num_digits_rounding': 2,
				},
			"rrna_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'rRNA RNAP Counts',
				},
			"rrna_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'rRNA RNAP Portion',
				 'num_digits_rounding': 4,
				},
			"rnap_subunit_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNAP Subunit RNAP Counts',
				},
			"rnap_subunit_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'RNAP Subunit RNAP Portion',
				 'num_digits_rounding': 4,
				},
			"rnap_subunit_ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'RNAP Subunit Ribosome Counts',
				},
			"rnap_subunit_ribosome_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'RNAP Subunit Ribosome Portion',
				 'num_digits_rounding': 4,
				},
			"ribosomal_protein_rnap_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Ribosomal Protein RNAP Counts',
				},
			"ribosomal_protein_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Ribosomal Protein RNAP Portion',
				 'num_digits_rounding': 4,
				},
			"ribosomal_protein_ribosome_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				 'plot_title': 'Ribosomal Protein Ribosome Counts',
				},
			"ribosomal_protein_ribosome_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Ribosomal Protein Ribosome Portion',
				 'num_digits_rounding': 4,
				},
			"rna_synth_prob_max_p_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'RNA Synth Prob Max P',
				 'num_digits_rounding': 4,
				},
			"protein_init_prob_max_p_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Protein Init Prob Max P',
				 'num_digits_rounding': 4,
				},
			"capacity_gene_mRNA_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene mRNA Counts: ' + capacity_gene_common_name,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				},
			"capacity_gene_monomer_counts_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene Monomer Counts: ' + capacity_gene_common_name,
				 'num_digits_rounding': 0,
				 'box_text_size': 'x-small',
				},
			"capacity_gene_mRNA_mass_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene mRNA Mass Fraction: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				},
			"capacity_gene_monomer_mass_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene Monomer Mass Fraction: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				},
			"capacity_gene_mRNA_counts_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene mRNA Counts Fraction: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				},
			"capacity_gene_monomer_counts_fraction_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene Monomer Counts Fraction: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				},
			"capacity_gene_rnap_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene RNAP Portion: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				},
			"capacity_gene_ribosome_portion_heatmap":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title':
				 'Capacity Gene Ribosome Portion: ' + capacity_gene_common_name,
				 'num_digits_rounding': 4,
				},
			"new_gene_yield_per_glucose":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'New Gene fg Protein Yield per fg Glucose',
				 'num_digits_rounding': 3,
				},
			"new_gene_yield_per_hour":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'New Gene fg Protein Yield per Hour',
				 'num_digits_rounding': 2,
				},
			"glucose_consumption_rate":
				{'is_nonstandard_data_retrieval': True,
				 'plot_title': 'Average Glucose Consumption Rate (fg/hr)',
				 'num_digits_rounding': 1,
				},
		}

		# Check validity and fill in defaults
		self.total_heatmaps_to_make = 0
		for h in heatmaps_to_make:
			assert h in self.heatmap_details, f"Heatmap {h} is not an option"
			self.heatmap_details[h]['is_new_gene_heatmap'] = h.startswith("new_gene_")
			self.heatmap_details[h].setdefault(
				'is_nonstandard_data_retrieval',
				default_is_nonstandard_data_retrieval)
			self.heatmap_details[h].setdefault(
				'is_nonstandard_plot', default_is_nonstandard_plot)
			self.heatmap_details[h].setdefault('default_value', default_value)
			self.heatmap_details[h].setdefault(
				'remove_first', default_remove_first)
			self.heatmap_details[h].setdefault(
				'function_to_apply', default_function_to_apply)
			self.heatmap_details[h].setdefault(
				'num_digits_rounding', default_num_digits_rounding)
			self.heatmap_details[h].setdefault(
				'box_text_size', default_box_text_size)
			if not h.startswith("new_gene_"):
				self.total_heatmaps_to_make += 1
			elif h == "new_gene_mRNA_NTP_fraction_heatmap":
				self.ntp_ids = list(
					self.sim_data.ntp_code_to_id_ordered.values())
				self.total_heatmaps_to_make += len(self.ntp_ids)
			else:
				self.total_heatmaps_to_make += len(self.new_gene_mRNA_ids)

		# Create data structures for heatmaps
		# Dimensions: (initial_index, trl_eff_index, exp_index)
		# where trl_eff_index is local_index (0=control, 1..N_TRL_EFF)
		# and exp_index is 0=no-kcat, 1=kcat
		self.heatmap_data = {}
		self.heatmap_data["completed_gens_heatmap"] = np.zeros((
			1, n_rows, n_cols))
		for h in heatmaps_to_make:
			if not self.heatmap_details[h]['is_new_gene_heatmap']:
				self.heatmap_data[h] = {}
				self.heatmap_data[h]["mean"] = np.zeros((
					1, n_rows, n_cols)
					) + self.heatmap_details[h]['default_value']
				self.heatmap_data[h]["std_dev"] = np.zeros((
					1, n_rows, n_cols)
					) + self.heatmap_details[h]['default_value']
			else:
				if h == "new_gene_mRNA_NTP_fraction_heatmap":
					self.heatmap_data[h] = {}
					self.heatmap_data[h]["mean"] = {}
					self.heatmap_data[h]["std_dev"] = {}
					for ntp_id in self.ntp_ids:
						self.heatmap_data[h]["mean"][ntp_id] = np.zeros(
							(len(self.new_gene_mRNA_ids), n_rows, n_cols)
							) + self.heatmap_details[h]['default_value']
						self.heatmap_data[h]["std_dev"][ntp_id] = np.zeros(
							(len(self.new_gene_mRNA_ids), n_rows, n_cols)
							) + self.heatmap_details[h]['default_value']
				else:
					self.heatmap_data[h] = {}
					self.heatmap_data[h]["mean"] = np.zeros((
						len(self.new_gene_mRNA_ids), n_rows, n_cols)
						) + self.heatmap_details[h]['default_value']
					self.heatmap_data[h]["std_dev"] = np.zeros((
						len(self.new_gene_mRNA_ids), n_rows, n_cols)
						) + self.heatmap_details[h]['default_value']

		# Data extraction
		print("---Data Extraction---")
		reached_count_gen = {}
		for variant in variants:
			print("Variant: ", variant)
			all_cells = self.ap.get_cells(
				variant=[variant],
				generation=np.arange(MIN_CELL_INDEX, MAX_CELL_INDEX),
				only_successful=True)
			exp_index, trl_eff_index = variant_index_to_list_indices[variant]
			if len(all_cells) == 0:
				continue

			# Completed Gens Heatmap
			num_count_gen = len(self.ap.get_cells(
				variant=[variant], generation=[COUNT_INDEX],
				only_successful=True))
			num_zero_gen = len(self.ap.get_cells(
				variant=[variant], generation=[0], only_successful=True))
			reached_count_gen[variant] = (
				num_count_gen / num_zero_gen if num_zero_gen > 0 else 0)
			self.heatmap_data["completed_gens_heatmap"][
				0, trl_eff_index, exp_index] = round(
				reached_count_gen[variant], 2)

			# Extract data for each heatmap
			for h in heatmaps_to_make:
				self.extract_heatmap_data(
					all_cells, h, trl_eff_index, exp_index)

		# Plotting
		print("---Plotting---")
		plot_suffix = "_gens_" + str(MIN_CELL_INDEX) + "_through_" + str(
			MAX_CELL_INDEX)
		heatmap_x_label = "Kcat Constraint"
		heatmap_y_label = "Translation Efficiency"
		figsize_x = 2 + 2 * n_cols / 3
		figsize_y = 2 * n_rows / 2

		if DASHBOARD_FLAG == 0 or DASHBOARD_FLAG == 2:
			self.plot_heatmaps(
				False, variant_mask, heatmap_x_label, heatmap_y_label,
				col_labels, row_labels, 'mean', figsize_x,
				figsize_y, plotOutDir, plot_suffix)

			if STD_DEV_FLAG:
				self.plot_heatmaps(
					False, variant_mask, heatmap_x_label, heatmap_y_label,
					col_labels, row_labels, 'std_dev',
					figsize_x, figsize_y, plotOutDir, plot_suffix)

		if DASHBOARD_FLAG == 1 or DASHBOARD_FLAG == 2:
			self.plot_heatmaps(
				True, variant_mask, heatmap_x_label, heatmap_y_label,
				col_labels, row_labels, 'mean', figsize_x,
				figsize_y, plotOutDir, plot_suffix)

			if STD_DEV_FLAG:
				self.plot_heatmaps(
					True, variant_mask, heatmap_x_label, heatmap_y_label,
					col_labels, row_labels, 'std_dev',
					figsize_x, figsize_y, plotOutDir, plot_suffix)

	def plot_heatmaps(
			self, is_dashboard, variant_mask, heatmap_x_label, heatmap_y_label,
			new_gene_expression_factors, new_gene_translation_efficiency_values,
			summary_statistic, figsize_x, figsize_y, plotOutDir, plot_suffix):
		"""
		Plots all heatmaps. Overrides parent to use the trl_eff_sweep axis
		labels (col_labels as expression_factors, row_labels as trl_eff_values).
		"""
		if summary_statistic == 'std_dev':
			plot_suffix = plot_suffix + "_std_dev"
		elif summary_statistic != 'mean':
			raise Exception(
				"'mean' and 'std_dev' are the only currently supported"
				" summary statistics")

		heatmaps_to_make = set(HEATMAPS_TO_MAKE_LIST)

		if is_dashboard:
			if self.total_heatmaps_to_make > 3:
				dashboard_ncols = 4
				dashboard_nrows = math.ceil(
					(self.total_heatmaps_to_make + 1) / dashboard_ncols)
			else:
				dashboard_ncols = self.total_heatmaps_to_make + 1
				dashboard_nrows = 1
			fig, axs = plt.subplots(
				nrows=dashboard_nrows,
				ncols=dashboard_ncols,
				figsize=(figsize_y * dashboard_ncols,
						 figsize_x * dashboard_nrows),
				layout='constrained')
			if dashboard_nrows == 1:
				axs = np.reshape(axs, (1, dashboard_ncols))

			# Percent Completion Heatmap
			heatmap(
				self, axs[0, 0], variant_mask,
				self.heatmap_data["completed_gens_heatmap"][0, :, :],
				self.heatmap_data["completed_gens_heatmap"][0, :, :],
				new_gene_expression_factors,
				new_gene_translation_efficiency_values,
				heatmap_x_label,
				heatmap_y_label,
				f"Percentage of Sims That Reached Generation {COUNT_INDEX + 1}")
			row_ax = 0
			col_ax = 1

			for h in HEATMAPS_TO_MAKE_LIST:
				if not self.heatmap_details[h]["is_nonstandard_plot"]:
					stop_index = 1
					title_addition = ""
					if self.heatmap_details[h]["is_new_gene_heatmap"]:
						stop_index = len(self.new_gene_mRNA_ids)
					for i in range(stop_index):
						if self.heatmap_details[h]["is_new_gene_heatmap"]:
							title_addition = f": {self.new_gene_mRNA_ids[i][:-4]}"
						self.make_single_heatmap(
							h, axs[row_ax, col_ax], variant_mask,
							heatmap_x_label, heatmap_y_label, i,
							new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							summary_statistic, title_addition)
						col_ax += 1
						if col_ax == dashboard_ncols:
							col_ax = 0
							row_ax += 1
				elif h == "new_gene_mRNA_NTP_fraction_heatmap":
					for i in range(len(self.new_gene_mRNA_ids)):
						for ntp_id in self.ntp_ids:
							self.make_new_gene_mRNA_NTP_fraction_heatmap(
								h, axs[row_ax, col_ax], variant_mask,
								heatmap_x_label, heatmap_y_label, i,
								new_gene_expression_factors,
								new_gene_translation_efficiency_values,
								summary_statistic, ntp_id)
							fig.tight_layout()
							plt.show()
							exportFigure(
								plt, plotOutDir,
								f'new_gene_mRNA_{ntp_id[:-3]}_fraction_heatmap'
								f'_{self.new_gene_mRNA_ids[i][:-4]}'
								f'{plot_suffix}')
							col_ax += 1
							if col_ax == dashboard_ncols:
								col_ax = 0
								row_ax += 1
				else:
					raise Exception(
						f"Heatmap {h} is neither a standard plot nor a"
						f" nonstandard plot that has specific instructions for"
						f" plotting.")
			fig.tight_layout()
			exportFigure(plt, plotOutDir,
				f"00_trl_eff_sweep_dashboard{plot_suffix}")
			plt.close("all")

		else:
			# Plot percent completion heatmap
			fig, ax = plt.subplots(1, 1, figsize=(figsize_x, figsize_y))
			heatmap(
				self, ax, variant_mask,
				self.heatmap_data["completed_gens_heatmap"][0, :, :],
				self.heatmap_data["completed_gens_heatmap"][0, :, :],
				new_gene_expression_factors,
				new_gene_translation_efficiency_values,
				heatmap_x_label,
				heatmap_y_label,
				f"Percentage of Sims that Reached Generation {COUNT_INDEX + 1}")
			fig.tight_layout()
			plt.show()
			exportFigure(plt, plotOutDir, 'completed_gens_heatmap')

			for h in HEATMAPS_TO_MAKE_LIST:
				if not self.heatmap_details[h]["is_nonstandard_plot"]:
					stop_index = 1
					title_addition = ""
					filename_addition = ""
					if self.heatmap_details[h]["is_new_gene_heatmap"]:
						stop_index = len(self.new_gene_mRNA_ids)
					for i in range(stop_index):
						if self.heatmap_details[h]["is_new_gene_heatmap"]:
							title_addition = f": {self.new_gene_mRNA_ids[i][:-4]}"
							filename_addition = f"_{self.new_gene_mRNA_ids[i][:-4]}"
						fig, ax = plt.subplots(1, 1,
							figsize=(figsize_x, figsize_y))
						self.make_single_heatmap(
							h, ax, variant_mask, heatmap_x_label,
							heatmap_y_label, i,
							new_gene_expression_factors,
							new_gene_translation_efficiency_values,
							summary_statistic, title_addition)
						fig.tight_layout()
						plt.show()
						exportFigure(plt, plotOutDir,
							h + filename_addition + plot_suffix)
						plt.close()
				elif h == "new_gene_mRNA_NTP_fraction_heatmap":
					for i in range(len(self.new_gene_mRNA_ids)):
						for ntp_id in self.ntp_ids:
							fig, ax = plt.subplots(1, 1,
								figsize=(figsize_x, figsize_y))
							self.make_new_gene_mRNA_NTP_fraction_heatmap(
								h, ax, variant_mask, heatmap_x_label,
								heatmap_y_label, i,
								new_gene_expression_factors,
								new_gene_translation_efficiency_values,
								summary_statistic, ntp_id)
							fig.tight_layout()
							plt.show()
							exportFigure(
								plt, plotOutDir,
								f'new_gene_mRNA_{ntp_id[:-3]}_fraction_heatmap'
								f'_{self.new_gene_mRNA_ids[i][:-4]}'
								f'{plot_suffix}')
				else:
					raise Exception(
						f"Heatmap {h} is neither a standard plot nor a"
						f" nonstandard plot that has specific instructions for"
						f" plotting.")


if __name__ == "__main__":
	Plot().cli()
