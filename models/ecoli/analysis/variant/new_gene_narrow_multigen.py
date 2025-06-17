"""
Analyze how production machinery counts change over time and across variants
"""

import pickle
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, read_stacked_bulk_molecules)
from wholecell.io.tablereader import TableReader

poster_colors = {
    "light_gray": (0.75, 0.75, 0.75),
    "poster_green": (66/255, 170/255, 154/255),
    "poster_blue": (27/255, 132/255, 198/255),
    "poster_purple": (188/255, 140/255, 191/255),
    "poster_gold": (221/255, 203/255, 119/255),
    "poster_light_blue": (136/255, 205/255, 240/255),
    "poster_red": (202/255, 0/255, 32/255),
}

START_GEN_INDEX = 7
END_GEN_INDEX = 13

SELECTED_SEED_INDEXES = [
			0, 1, 2, 3, 4, 5, 6, 7]

# For SHEROCK
SELECTED_VARIANT_INDEXES = [
	0,  # Expression factor: 0, Translation efficiency: 0, Doubling time: 53
	23,  # Expression factor: 8.5, Translation efficiency: 0, Doubling time: 50
	8,  # Expression factor: 8.5, Translation Efficiency: 2.20, Doubling time: 66
	1,  # Expression factor: 8.5, Translation Efficiency: 3.25, Doubling time: 78
	]
VAR_TO_COLOR = {
	0: poster_colors["poster_gold"],
	1: poster_colors["poster_green"],
	8: poster_colors["poster_purple"],
	23: poster_colors["poster_blue"],
	}

# FOR LOCAL TESTING
SELECTED_VARIANT_INDEXES = [0, 20, 21]
VAR_TO_COLOR = {
	0: poster_colors["poster_gold"],
	20: poster_colors["poster_green"],
	21: poster_colors["poster_purple"],
	}



class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		plotOutDir = os.path.join(plotOutDir, plotOutFileName)

		all_cells_dict = {}
		all_time_dict = {}
		all_time_no_first_dict = {}
		all_num_time_steps_dict = {}
		all_gen_start_index_dict = {}
		all_gen_end_index_dict = {}
		all_total_rnap_counts_ratio_dict = {}
		all_total_ribosome_counts_ratio_dict = {}
		for seed_index in SELECTED_SEED_INDEXES:
			all_cells_dict[seed_index] = {}
			all_time_dict[seed_index] = {}
			all_time_no_first_dict[seed_index] = {}
			all_num_time_steps_dict[seed_index] = {}
			all_gen_start_index_dict[seed_index] = {}
			all_gen_end_index_dict[seed_index] = {}
			all_total_rnap_counts_ratio_dict[seed_index] = {}
			all_total_ribosome_counts_ratio_dict[seed_index] = {}

		# Loop through variant indexes
		for i, variant_index in enumerate(SELECTED_VARIANT_INDEXES):
			for seed_index in SELECTED_SEED_INDEXES:
				all_cells = self.ap.get_cells(
					variant=[variant_index],
					seed=[seed_index],
					generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
					only_successful=True)
				if len(all_cells) == 0:
					continue

				time = read_stacked_columns(
					all_cells, 'Main', 'time',
					ignore_exception=True).squeeze()
				time_no_first = read_stacked_columns(
					all_cells, 'Main', 'time',
					remove_first=True, ignore_exception=True).squeeze()
				num_time_steps = read_stacked_columns(
					all_cells, 'Main', 'time',
					fun=lambda x: len(x)).squeeze()
				dt = read_stacked_columns(
					all_cells, 'Main', 'time',
					fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
				gen_labels = np.repeat(np.arange(len(dt)), num_time_steps)
				unique_gen_labels = np.unique(gen_labels)
				gen_start_index = np.array(
					[gen_labels.tolist().index(i) for i in unique_gen_labels])
				gen_end_index = np.concatenate((
					np.array(gen_start_index[1:] - 1), np.array([len(gen_labels) - 1])))

				all_cells_dict[seed_index][variant_index] = all_cells
				all_time_dict[seed_index][variant_index] = time
				all_time_no_first_dict[seed_index][variant_index] = time_no_first
				all_num_time_steps_dict[seed_index][variant_index] = num_time_steps
				all_gen_start_index_dict[seed_index][variant_index] = gen_start_index
				all_gen_end_index_dict[seed_index][variant_index] = gen_end_index

		plot_types = ["counts", "ratios"]
		for seed_index in SELECTED_SEED_INDEXES:
			variants_set = all_cells_dict[seed_index].keys()
			all_cells_set = all_cells_dict[seed_index]
			time_set = all_time_dict[seed_index]
			time_no_first_set = all_time_no_first_dict[seed_index]
			num_time_steps_set = all_num_time_steps_dict[seed_index]

			# if all_cells_set is empty, skip
			if len(all_cells_set) == 0:
				continue

			for plot_type in plot_types:
				# Create figure
				total_plots = 30 # TODO: Modularize and get rid of this magic number
				mpl.rcParams['axes.spines.right'] = False
				mpl.rcParams['axes.spines.top'] = False
				plt.figure(figsize = (14, total_plots*3))

				plot_num = 1
				ax1 = plt.subplot(total_plots, 1, plot_num)

				# Doubling Time
				for variant_index in variants_set:
					cell_paths = all_cells_set[variant_index]
					time = time_set[variant_index]
					num_time_steps = num_time_steps_set[variant_index]
					color = VAR_TO_COLOR[variant_index]
					dt = read_stacked_columns(
						cell_paths, 'Main', 'time',
						fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
					dt_to_plot = np.repeat(dt, num_time_steps)
					plt.plot(
						time / 60.0, dt_to_plot, color=color, linewidth=1,
						label=f"Variant {variant_index}")
				plt.ylabel("Doubling Time (min)", fontsize=8)
				plt.xlabel("Time (min)", fontsize=8)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				plot_num += 1

				# Cell mass
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				for variant_index in variants_set:
					cell_paths = all_cells_set[variant_index]
					time = time_set[variant_index]
					color = VAR_TO_COLOR[variant_index]
					cell_mass = read_stacked_columns(
						cell_paths, 'Mass', 'cellMass',
						ignore_exception=True).squeeze()
					if plot_type == "counts":
						plt.plot(time / 60.0, cell_mass, color=color, linewidth=1,
								 label=f"Variant {variant_index}")
						plt.ylabel("Cell Mass (fg)", fontsize=8)
					elif plot_type == "ratios":
						num_time_steps = num_time_steps_set[variant_index]
						gen_starts = all_gen_start_index_dict[seed_index][variant_index]
						initial_cell_mass = cell_mass[gen_starts]
						initial_cell_mass_vec = np.repeat(initial_cell_mass, num_time_steps)
						cell_mass_ratio = cell_mass / initial_cell_mass_vec
						plt.plot(time / 60.0, cell_mass_ratio, color=color, linewidth=1,
								 label=f"Variant {variant_index}")
						plt.axhline(y=2, color='k', linestyle='--', linewidth=0.5)
						plt.ylabel("Initial Cell Mass Ratio", fontsize=8)
						plt.ylim(0.9, 2.5)
				plt.xlabel("Time (min)", fontsize=8)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				plot_num += 1

				# Total RNAP counts
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				for variant_index in variants_set:
					cell_paths = all_cells_set[variant_index]
					sim_dir = cell_paths[0]
					simOutDir = os.path.join(sim_dir, 'simOut')
					time = time_set[variant_index]
					color = VAR_TO_COLOR[variant_index]
					# Inactive
					rnap_id = [sim_data.molecule_ids.full_RNAP]
					(inactive_rnap_counts,) = read_stacked_bulk_molecules(
						cell_paths, (rnap_id,), ignore_exception=True)
					# Active
					uniqueMoleculeCounts = TableReader(
						os.path.join(simOutDir, "UniqueMoleculeCounts"))
					active_rnap_index = uniqueMoleculeCounts.readAttribute(
						"uniqueMoleculeIds").index('active_RNAP')
					active_rnap_counts = read_stacked_columns(
						cell_paths, 'UniqueMoleculeCounts',
						'uniqueMoleculeCounts',
						ignore_exception=True)[:, active_rnap_index]
					total_rnap_counts = inactive_rnap_counts + active_rnap_counts
					if plot_type == "counts":
						plt.plot(time / 60.0, total_rnap_counts, color=color, linewidth=1,
								 label=f"Variant {variant_index}")
						plt.ylabel("Total RNAP Counts", fontsize=8)
					elif plot_type == "ratios":
						num_time_steps = num_time_steps_set[variant_index]
						gen_starts = all_gen_start_index_dict[seed_index][variant_index]
						initial_inactive_rnap_counts = inactive_rnap_counts[
							gen_starts]
						initial_active_rnap_counts = active_rnap_counts[
							gen_starts]
						initial_total_rnap_counts = (
							initial_inactive_rnap_counts + initial_active_rnap_counts)
						initial_total_rnap_counts_vec = np.repeat(initial_total_rnap_counts, num_time_steps)
						total_rnap_counts_ratio = total_rnap_counts / initial_total_rnap_counts_vec
						plt.plot(time / 60.0, total_rnap_counts_ratio, color=color, linewidth=1,
								 label=f"Variant {variant_index}")
						plt.axhline(y=2, color='k', linestyle='--', linewidth=0.5)
						plt.ylabel("Initial Total RNAP Counts Ratio", fontsize=8)
						all_total_rnap_counts_ratio_dict[seed_index][variant_index] = total_rnap_counts_ratio
						plt.ylim(0.9, 2.5)
				plt.xlabel("Time (min)", fontsize=8)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				plot_num += 1

				# Total Ribosome counts
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				for variant_index in variants_set:
					cell_paths = all_cells_set[variant_index]
					sim_dir = cell_paths[0]
					simOutDir = os.path.join(sim_dir, 'simOut')
					time = time_set[variant_index]
					color = VAR_TO_COLOR[variant_index]
					# Inactive
					complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
					complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
					(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
						cell_paths, (complex_id_30s, complex_id_50s), ignore_exception=True)
					inactive_ribosome_counts = np.minimum(
						complex_counts_30s, complex_counts_50s)
					# Active
					unique_molecule_counts_table = TableReader(
						os.path.join(simOutDir, "UniqueMoleculeCounts"))
					ribosome_index = unique_molecule_counts_table.readAttribute(
						"uniqueMoleculeIds").index('active_ribosome')
					active_ribosome_counts = read_stacked_columns(
						cell_paths, 'UniqueMoleculeCounts',
						'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
					# Total
					total_ribosome_counts = inactive_ribosome_counts + active_ribosome_counts
					if plot_type == "counts":
						plt.plot(time / 60.0, total_ribosome_counts, color=color, linewidth=1,
								 label=f"Variant {variant_index}")
						plt.ylabel("Total Ribosome Counts", fontsize=8)
					elif plot_type == "ratios":
						num_time_steps = num_time_steps_set[variant_index]
						gen_starts = all_gen_start_index_dict[seed_index][variant_index]
						initial_inactive_ribosome_counts = inactive_ribosome_counts[
							gen_starts]
						initial_active_ribosome_counts = active_ribosome_counts[
							gen_starts]
						initial_total_ribosome_counts = (
							initial_inactive_ribosome_counts + initial_active_ribosome_counts)
						initial_total_ribosome_counts_vec = np.repeat(initial_total_ribosome_counts, num_time_steps)
						total_ribosome_counts_ratio = total_ribosome_counts / initial_total_ribosome_counts_vec
						plt.plot(time / 60.0, total_ribosome_counts_ratio, color=color, linewidth=1,
								 label=f"Variant {variant_index}")
						plt.axhline(y=2, color='k', linestyle='--', linewidth=0.5)
						plt.ylabel("Initial Total Ribosome Counts Ratio", fontsize=8)
						all_total_ribosome_counts_ratio_dict[seed_index][variant_index] = total_ribosome_counts_ratio
						plt.ylim(0.9, 2.5)
				plt.xlabel("Time (min)", fontsize=8)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				plot_num += 1

				# Plot total RNAP - total Ribosome counts ratios
				if plot_type == "ratios":
					plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					for variant_index in variants_set:
						time = time_set[variant_index]
						color = VAR_TO_COLOR[variant_index]
						total_rnap_counts_ratio = all_total_rnap_counts_ratio_dict[seed_index][variant_index]
						total_ribosome_counts_ratio = all_total_ribosome_counts_ratio_dict[seed_index][variant_index]
						plt.plot(
							time / 60.0, total_rnap_counts_ratio - total_ribosome_counts_ratio,
							color=color, linewidth=1, label=f"Variant {variant_index}")
						plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
						plt.ylabel("Initial Total RNAP - Total Ribosome Counts Ratio", fontsize=8)
					plt.xlabel("Time (min)", fontsize=8)
					plt.legend(fontsize=6, loc='lower left', ncol=2)
					plot_num += 1

				# Plot difference between total RNAP ratio and the total RNAP ratio for variant 0
				if plot_type == "ratios":
					plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					for variant_index in variants_set:
						time = time_set[variant_index]
						color = VAR_TO_COLOR[variant_index]
						total_rnap_counts_ratio = all_total_rnap_counts_ratio_dict[seed_index][variant_index]
						total_rnap_counts_ratio_0 = all_total_rnap_counts_ratio_dict[seed_index][0]
						min_length = min(len(total_rnap_counts_ratio), len(total_rnap_counts_ratio_0))
						total_rnap_counts_ratio = total_rnap_counts_ratio[:min_length]
						total_rnap_counts_ratio_0 = total_rnap_counts_ratio_0[:min_length]
						if len(time) > min_length:
							time = time[:min_length]
						plt.plot(
							time / 60.0, total_rnap_counts_ratio - total_rnap_counts_ratio_0,
							color=color, linewidth=1, label=f"Variant {variant_index}")
						plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
						plt.ylabel("Initial Total RNAP Counts Ratio Diff", fontsize=8)
					plt.xlabel("Time (min)", fontsize=8)
					plt.legend(fontsize=6, loc='lower left', ncol=2)
					plt.ylim(-1.0, 1.0)
					plot_num += 1

				# Plot difference between total Ribosome ratio and the total Ribosome ratio for variant 0
				if plot_type == "ratios":
					plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					for variant_index in variants_set:
						time = time_set[variant_index]
						color = VAR_TO_COLOR[variant_index]
						total_ribosome_counts_ratio = all_total_ribosome_counts_ratio_dict[seed_index][variant_index]
						total_ribosome_counts_ratio_0 = all_total_ribosome_counts_ratio_dict[seed_index][0]
						min_length = min(len(total_ribosome_counts_ratio), len(total_ribosome_counts_ratio_0))
						total_ribosome_counts_ratio = total_ribosome_counts_ratio[:min_length]
						total_ribosome_counts_ratio_0 = total_ribosome_counts_ratio_0[:min_length]
						# Ensure time vectors are the same length
						if len(time) > min_length:
							time = time[:min_length]
						plt.plot(
							time / 60.0, total_ribosome_counts_ratio - total_ribosome_counts_ratio_0,
							color=color, linewidth=1, label=f"Variant {variant_index}")
						plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
						plt.ylabel("Initial Total Ribosome Counts Ratio Diff", fontsize=8)
					plt.xlabel("Time (min)", fontsize=8)
					plt.legend(fontsize=6, loc='lower left', ncol=2)
					plt.ylim(-1.0, 1.0)
					plot_num += 1

					if plot_type == "ratios":
						plt.subplot(total_plots, 1, plot_num, sharex=ax1)
						for variant_index in variants_set:
							time = time_set[variant_index]
							color = VAR_TO_COLOR[variant_index]
							total_rnap_counts_ratio = all_total_rnap_counts_ratio_dict[seed_index][
								variant_index]
							total_rnap_counts_ratio_0 = all_total_rnap_counts_ratio_dict[seed_index][0]
							min_length = min(len(total_rnap_counts_ratio),
											 len(total_rnap_counts_ratio_0))
							total_rnap_counts_ratio = total_rnap_counts_ratio[:min_length]
							total_rnap_counts_ratio_0 = total_rnap_counts_ratio_0[:min_length]
							total_ribosome_counts_ratio = all_total_ribosome_counts_ratio_dict[seed_index][variant_index]
							total_ribosome_counts_ratio_0 = all_total_ribosome_counts_ratio_dict[seed_index][0]
							total_ribosome_counts_ratio = total_ribosome_counts_ratio[:min_length]
							total_ribosome_counts_ratio_0 = total_ribosome_counts_ratio_0[:min_length]
							if len(time) > min_length:
								time = time[:min_length]
							plt.plot(
								time / 60.0,
								(total_rnap_counts_ratio - total_rnap_counts_ratio_0) - (total_ribosome_counts_ratio - total_ribosome_counts_ratio_0),
								color=color, linewidth=1, label=f"Variant {variant_index}")
							plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
							plt.ylabel("RNAP - Ribo Counts Ratio Diff to Gen Initial 0", fontsize=8)
						plt.xlabel("Time (min)", fontsize=8)
						plt.legend(fontsize=6, loc='lower left', ncol=2)
						plot_num += 1

						# Doubling Time
						plt.subplot(total_plots, 1, plot_num, sharex=ax1)
						for variant_index in variants_set:
							cell_paths = all_cells_set[variant_index]
							time = time_set[variant_index]
							num_time_steps = num_time_steps_set[variant_index]
							color = VAR_TO_COLOR[variant_index]
							dt = read_stacked_columns(
								cell_paths, 'Main', 'time',
								fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
							dt_to_plot = np.repeat(dt, num_time_steps)
							plt.plot(
								time / 60.0, dt_to_plot, color=color, linewidth=1,
								label=f"Variant {variant_index}")
						plt.ylabel("Doubling Time (min)", fontsize=8)
						plt.xlabel("Time (min)", fontsize=8)
						plt.legend(fontsize=6, loc='upper left', ncol=2)
						plot_num += 1

					# RNAP Subunit Protein Counts
					RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
					monomer_counts_reader = TableReader(
						os.path.join(simOutDir, "MonomerCounts"))
					monomer_idx_dict = {
						monomer: i for i, monomer in enumerate(
						monomer_counts_reader.readAttribute('monomerIds'))}
					for i, RNAP_subunit_monomer_id in enumerate(RNAP_subunit_monomer_ids):
						RNAP_subunit_monomer_index = monomer_idx_dict.get(
							RNAP_subunit_monomer_id)
						if RNAP_subunit_monomer_index is None:
							print(
								f"Monomer {RNAP_subunit_monomer_id} not found in MonomerCounts.")
							continue
						plt.subplot(total_plots, 1, plot_num, sharex=ax1)
						for variant_index in variants_set:
							cell_paths = all_cells_dict[seed_index][variant_index]
							time = all_time_dict[seed_index][variant_index]
							num_time_steps = all_num_time_steps_dict[seed_index][variant_index]
							color = VAR_TO_COLOR[variant_index]
							RNAP_subunit_monomer_counts = read_stacked_columns(
								cell_paths, 'MonomerCounts', 'monomerCounts',
								ignore_exception=True)[:, RNAP_subunit_monomer_index].squeeze()
							if RNAP_subunit_monomer_id == "EG10893-MONOMER[c]":
								# Divide by 2 to account for stochiometry of RNAP subunits
								RNAP_subunit_monomer_counts = RNAP_subunit_monomer_counts / 2.0
							if plot_type == "counts":
								# Plot the RNAP subunit monomer counts over time
								plt.plot(
									time / 60.0, RNAP_subunit_monomer_counts, color=color,
									linewidth=1,
									label=f"Variant {variant_index}")
								plt.ylabel(f"{RNAP_subunit_monomer_id} Counts", fontsize=8)
							elif plot_type == "ratios":
								# Plot the ratio of RNAP subunit monomer counts to the initial counts
								# at the start of each generation
								gen_starts = all_gen_start_index_dict[seed_index][variant_index]
								initial_monomer_counts = RNAP_subunit_monomer_counts[gen_starts]
								initial_monomer_counts_vec = np.repeat(initial_monomer_counts,
																	   num_time_steps)
								monomer_counts_ratio = RNAP_subunit_monomer_counts / initial_monomer_counts_vec
								plt.plot(
									time / 60.0, monomer_counts_ratio, color=color, linewidth=1,
									label=f"Variant {variant_index}")
								plt.axhline(y=2, color='k', linestyle='--', linewidth=0.5)
								plt.ylabel(f"Initial {RNAP_subunit_monomer_id} Counts Ratio",
										   fontsize=8)
								plt.ylim(0.9, 2.5)
						plt.xlabel("Time (min)", fontsize=8)
						plt.legend(fontsize=6, loc='upper left', ncol=2)
						plot_num += 1

					# RNAP Subunit Protein Counts Ratio Difference
					if plot_type == "ratios":
						for i, RNAP_subunit_monomer_id in enumerate(RNAP_subunit_monomer_ids):
							RNAP_subunit_monomer_index = monomer_idx_dict.get(
								RNAP_subunit_monomer_id)
							if RNAP_subunit_monomer_index is None:
								print(
									f"Monomer {RNAP_subunit_monomer_id} not found in MonomerCounts.")
								continue
							plt.subplot(total_plots, 1, plot_num, sharex=ax1)
							for variant_index in variants_set:
								time = all_time_dict[seed_index][variant_index]
								color = VAR_TO_COLOR[variant_index]
								RNAP_subunit_monomer_counts = read_stacked_columns(
									all_cells_dict[seed_index][variant_index], 'MonomerCounts',
									'monomerCounts', ignore_exception=True)[:, RNAP_subunit_monomer_index].squeeze()
								RNAP_subunit_monomer_counts_0 = read_stacked_columns(
									all_cells_dict[seed_index][0], 'MonomerCounts',
									'monomerCounts', ignore_exception=True)[:, RNAP_subunit_monomer_index].squeeze()
								if RNAP_subunit_monomer_id == "EG10893-MONOMER[c]":
									# Divide by 2 to account for stochiometry of RNAP subunits
									RNAP_subunit_monomer_counts = RNAP_subunit_monomer_counts / 2.0
									RNAP_subunit_monomer_counts_0 = RNAP_subunit_monomer_counts_0 / 2.0
								gen_starts = all_gen_start_index_dict[seed_index][variant_index]
								initial_RNAP_subunit_monomer_counts = RNAP_subunit_monomer_counts[
									gen_starts]
								initial_RNAP_subunit_monomer_counts_vec = np.repeat(
									initial_RNAP_subunit_monomer_counts, num_time_steps_set[variant_index])
								RNAP_subunit_monomer_counts_ratio = (
										RNAP_subunit_monomer_counts / initial_RNAP_subunit_monomer_counts_vec)
								initial_RNAP_subunit_monomer_counts_0 = RNAP_subunit_monomer_counts_0[
									gen_starts]
								initial_RNAP_subunit_monomer_counts_vec_0 = np.repeat(
									initial_RNAP_subunit_monomer_counts_0, num_time_steps_set[0])
								RNAP_subunit_monomer_counts_ratio_0 = (
									RNAP_subunit_monomer_counts_0 / initial_RNAP_subunit_monomer_counts_vec_0)
								min_length = min(len(RNAP_subunit_monomer_counts_ratio),
												 len(RNAP_subunit_monomer_counts_ratio_0))
								RNAP_subunit_monomer_counts_ratio = RNAP_subunit_monomer_counts_ratio[:min_length]
								RNAP_subunit_monomer_counts_ratio_0 = RNAP_subunit_monomer_counts_ratio_0[:min_length]
								if len(time) > min_length:
									time = time[:min_length]
								plt.plot(
									time / 60.0,
									RNAP_subunit_monomer_counts_ratio - RNAP_subunit_monomer_counts_ratio_0,
									color=color, linewidth=1, label=f"Variant {variant_index}")
							plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
							plt.ylabel(f"{RNAP_subunit_monomer_id} Counts Ratio Diff",
									   fontsize=8)
							plt.xlabel("Time (min)", fontsize=8)
							plt.legend(fontsize=6, loc='lower left', ncol=2)
							plot_num += 1

				# Save figure
				print(f"\nSeed: {seed_index}, Plot Type: {plot_type}")
				print("Total number of plots made: ", plot_num - 1)
				plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
				variants_str = "_".join([str(v) for v in variants_set])
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_{plot_type}_seed_{seed_index}_var_{variants_str}",
					metadata)
				plt.close("all")

			# RNAP and Ribosome Initial Ratio direct comparison - for one variant and one seed
			ribosome_color = poster_colors["poster_blue"]
			rnap_color = poster_colors["poster_red"]
			for variant_index in variants_set:
				cell_paths = all_cells_dict[seed_index][variant_index]
				time = all_time_dict[seed_index][variant_index]
				time_no_first = all_time_no_first_dict[seed_index][variant_index]
				num_time_steps = all_num_time_steps_dict[seed_index][variant_index]

				# Inactive RNAP
				rnap_id = [sim_data.molecule_ids.full_RNAP]
				(inactive_rnap_counts,) = read_stacked_bulk_molecules(
					cell_paths, (rnap_id,), ignore_exception=True)
				# Active RNAP
				uniqueMoleculeCounts = TableReader(
					os.path.join(cell_paths[0], "simOut", "UniqueMoleculeCounts"))
				active_rnap_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_RNAP')
				active_rnap_counts = read_stacked_columns(
					cell_paths, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts', ignore_exception=True)[:, active_rnap_index]
				total_rnap_counts = inactive_rnap_counts + active_rnap_counts
				gen_starts = all_gen_start_index_dict[seed_index][variant_index]
				initial_inactive_rnap_counts = inactive_rnap_counts[
					gen_starts]
				initial_active_rnap_counts = active_rnap_counts[
					gen_starts]
				initial_total_rnap_counts = (
						initial_inactive_rnap_counts + initial_active_rnap_counts)
				initial_total_rnap_counts_vec = np.repeat(initial_total_rnap_counts, num_time_steps)
				total_rnap_counts_ratio = total_rnap_counts / initial_total_rnap_counts_vec

				# Inactive Ribosome
				complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
				complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
				(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
					cell_paths, (complex_id_30s, complex_id_50s), ignore_exception=True)
				inactive_ribosome_counts = np.minimum(complex_counts_30s, complex_counts_50s)
				# Active Ribosome
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')
				active_ribosome_counts = read_stacked_columns(
					cell_paths, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
				total_ribosome_counts = inactive_ribosome_counts + active_ribosome_counts
				gen_starts = all_gen_start_index_dict[seed_index][variant_index]
				initial_inactive_ribosome_counts = inactive_ribosome_counts[
					gen_starts]
				initial_active_ribosome_counts = active_ribosome_counts[
					gen_starts]
				initial_total_ribosome_counts = (
						initial_inactive_ribosome_counts + initial_active_ribosome_counts)
				initial_total_ribosome_counts_vec = np.repeat(initial_total_ribosome_counts,
															  num_time_steps)
				total_ribosome_counts_ratio = total_ribosome_counts / initial_total_ribosome_counts_vec

				plt.figure(figsize=(8, 6))
				plt.plot(
					time / 60.0, total_rnap_counts_ratio, color=rnap_color,
					linewidth=1, label = "Total RNAP Counts")
				plt.plot(
					time / 60.0, total_ribosome_counts_ratio, color=ribosome_color,
					linewidth=1, label = "Total Ribosome Counts")
				plt.xlabel("Time (min)", fontsize=8)
				plt.ylabel("Counts Ratio to Gen Initial", fontsize=8)
				plt.axhline(y=2, color='k', linestyle='--', linewidth=0.5)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				plt.title(f"Variant {variant_index}, Seed {seed_index}", fontsize=8)
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_rnap_ribosome_ratio_seed_{seed_index}_var_{variant_index}",
					metadata)
				plt.close("all")

				plt.figure(figsize=(8, 6))
				plt.plot(
					time / 60.0, total_rnap_counts_ratio - total_ribosome_counts_ratio,
					color=poster_colors["poster_purple"], linewidth=1)
				plt.xlabel("Time (min)", fontsize=8)
				plt.ylabel("RNAP - Ribo Counts Ratio to Gen Initial", fontsize=8)
				plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
				plt.title(f"Variant {variant_index}, Seed {seed_index}", fontsize=8)
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_rnap_ribosome_ratio_diff_seed_{seed_index}_var_{variant_index}",
					metadata)
				plt.close("all")

				plt.figure(figsize=(8, 6))
				total_rnap_counts_ratio_0 = all_total_rnap_counts_ratio_dict[seed_index][0]
				total_ribosome_counts_ratio_0 = all_total_ribosome_counts_ratio_dict[seed_index][0]
				min_length = min(len(total_rnap_counts_ratio), len(total_rnap_counts_ratio_0))
				total_rnap_counts_ratio = total_rnap_counts_ratio[:min_length]
				total_rnap_counts_ratio_0 = total_rnap_counts_ratio_0[:min_length]
				total_ribosome_counts_ratio = total_ribosome_counts_ratio[:min_length]
				total_ribosome_counts_ratio_0 = total_ribosome_counts_ratio_0[:min_length]
				if len(time) > min_length:
					time = time[:min_length]
				plt.plot(
					time / 60.0, total_rnap_counts_ratio - total_rnap_counts_ratio_0,
					color=rnap_color, linewidth=1, label = "Total RNAP Counts Ratio Diff")
				plt.plot(
					time / 60.0, total_ribosome_counts_ratio - total_ribosome_counts_ratio_0,
					color=ribosome_color, linewidth=1, label = "Total Ribosome Counts Ratio Diff")
				plt.xlabel("Time (min)", fontsize=8)
				plt.ylabel("Counts Ratio Diff to Gen Initial", fontsize=8)
				plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				plt.title(f"Variant {variant_index}, Seed {seed_index}", fontsize=8)
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_rnap_ribosome_ratio_diff0_seed_{seed_index}_var_{variant_index}",
					metadata)
				plt.close("all")

				plt.figure(figsize=(8, 6))
				plt.plot(
					time / 60.0,
					(total_rnap_counts_ratio - total_rnap_counts_ratio_0) - (total_ribosome_counts_ratio - total_ribosome_counts_ratio_0),
					color=poster_colors["poster_purple"], linewidth=1)
				plt.xlabel("Time (min)", fontsize=8)
				plt.ylabel("RNAP - Ribo Counts Ratio Diff to Gen Initial 0", fontsize=8)
				plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
				plt.title(f"Variant {variant_index}, Seed {seed_index}", fontsize=8)
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_rnap_ribosome_ratio_diff0_diff_seed_{seed_index}_var_{variant_index}",
					metadata)
				plt.close("all")

				# Compare the protein counts of the RNAP subunits
				RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
				monomer_counts_reader = TableReader(
					os.path.join(cell_paths[0], "simOut", "MonomerCounts"))
				monomer_idx_dict = {
					monomer: i for i, monomer in enumerate(
						monomer_counts_reader.readAttribute('monomerIds'))}
				plt.figure(figsize=(8, 6))
				for i, RNAP_subunit_monomer_id in enumerate(RNAP_subunit_monomer_ids):
					RNAP_subunit_monomer_index = monomer_idx_dict.get(
						RNAP_subunit_monomer_id)
					if RNAP_subunit_monomer_index is None:
						print(
							f"Monomer {RNAP_subunit_monomer_id} not found in MonomerCounts.")
						continue
					RNAP_subunit_monomer_counts = read_stacked_columns(
						cell_paths, 'MonomerCounts', 'monomerCounts',
						ignore_exception=True)[:, RNAP_subunit_monomer_index].squeeze()
					RNAP_subunit_monomer_counts_0 = read_stacked_columns(
						all_cells_dict[seed_index][0], 'MonomerCounts',
						'monomerCounts', ignore_exception=True)[:, RNAP_subunit_monomer_index].squeeze()
					if RNAP_subunit_monomer_id == "EG10893-MONOMER[c]":
						# Divide by 2 to account for stochiometry of RNAP subunits
						RNAP_subunit_monomer_counts = RNAP_subunit_monomer_counts / 2.0
						RNAP_subunit_monomer_counts_0 = RNAP_subunit_monomer_counts_0 / 2.0
					gen_starts = all_gen_start_index_dict[seed_index][variant_index]
					initial_RNAP_subunit_monomer_counts = RNAP_subunit_monomer_counts[
						gen_starts]
					initial_RNAP_subunit_monomer_counts_vec = np.repeat(
						initial_RNAP_subunit_monomer_counts, num_time_steps)
					RNAP_subunit_monomer_counts_ratio = (
							RNAP_subunit_monomer_counts / initial_RNAP_subunit_monomer_counts_vec)
					initial_RNAP_subunit_monomer_counts_0 = RNAP_subunit_monomer_counts_0[
						gen_starts]
					initial_RNAP_subunit_monomer_counts_vec_0 = np.repeat(
						initial_RNAP_subunit_monomer_counts_0, num_time_steps_set[0])
					RNAP_subunit_monomer_counts_ratio_0 = (
						RNAP_subunit_monomer_counts_0 / initial_RNAP_subunit_monomer_counts_vec_0)
					min_length = min(len(RNAP_subunit_monomer_counts_ratio),
									 len(RNAP_subunit_monomer_counts_ratio_0))
					RNAP_subunit_monomer_counts_ratio = RNAP_subunit_monomer_counts_ratio[:min_length]
					RNAP_subunit_monomer_counts_ratio_0 = RNAP_subunit_monomer_counts_ratio_0[:min_length]
					if len(time) > min_length:
						time_plot = time[:min_length]
					else:
						time_plot = time
					plt.plot(
						time_plot / 60.,
						RNAP_subunit_monomer_counts_ratio - RNAP_subunit_monomer_counts_ratio_0,
						linewidth=1,
						label=f"{RNAP_subunit_monomer_id} Counts Ratio Diff")
				plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
				plt.xlabel("Time (min)", fontsize=8)
				plt.ylabel(f"RNAP Subunit Protein Counts Ratio Diff", fontsize=8)
				plt.title(f"Variant {variant_index}, Seed {seed_index}", fontsize=8)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_RNAP_subunit_ratio_diff_seed_{seed_index}_var_{variant_index}",
					metadata)
				plt.close("all")

				time = all_time_dict[seed_index][variant_index]
				# Plot the RNAP subunit monomer counts over time
				plt.figure(figsize=(8, 6))
				for i, RNAP_subunit_monomer_id in enumerate(RNAP_subunit_monomer_ids):
					RNAP_subunit_monomer_index = monomer_idx_dict.get(
						RNAP_subunit_monomer_id)
					if RNAP_subunit_monomer_index is None:
						print(
							f"Monomer {RNAP_subunit_monomer_id} not found in MonomerCounts.")
						continue
					RNAP_subunit_monomer_counts = read_stacked_columns(
						cell_paths, 'MonomerCounts', 'monomerCounts',
						ignore_exception=True)[:, RNAP_subunit_monomer_index].squeeze()
					if RNAP_subunit_monomer_id == "EG10893-MONOMER[c]":
						# Divide by 2 to account for stochiometry of RNAP subunits
						RNAP_subunit_monomer_counts = RNAP_subunit_monomer_counts / 2.0
					plt.plot(
						time / 60.0, RNAP_subunit_monomer_counts, linewidth=1,
						label=f"{RNAP_subunit_monomer_id} Counts")
				plt.plot(
					time / 60.0, total_rnap_counts, color="black",
					linewidth=1, label="Total RNAP Counts",
					linestyle='--', alpha=0.5)
				plt.xlabel("Time (min)", fontsize=8)
				plt.ylabel("Counts", fontsize=8)
				plt.title(f"Variant {variant_index}, Seed {seed_index}", fontsize=8)
				plt.legend(fontsize=6, loc='upper left', ncol=2)
				exportFigure(
					plt, plotOutDir,
					plotOutFileName + f"_RNAP_subunit_counts_seed_{seed_index}_var_{variant_index}",
					metadata)
				plt.close("all")


if __name__ == "__main__":
	Plot().cli()
