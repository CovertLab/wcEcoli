"""
Analyze how production machinery counts change over time and across at most two
variants.
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

START_GEN_INDEX = 8
END_GEN_INDEX = 9 # Not inclusive
GEN_RANGE = np.arange(START_GEN_INDEX, END_GEN_INDEX)

DEFAULT_COLOR = poster_colors["light_gray"]

VARIANT_CONTROL_COLOR = poster_colors["light_gray"]
VARIANT_A_COLOR = poster_colors["poster_blue"]
VARIANT_B_COLOR = poster_colors["poster_light_blue"]

# FOR SHERLOCK
VARIANT_INDEX_CONTROL = 0
VARIANT_INDEX_A = 1
VARIANT_INDEX_B = 2 # Optional, set to -1 if only one variant is desired

# # FOR LOCAL DEVELOPMENT
# VARIANT_INDEX_CONTROL = 0
# VARIANT_INDEX_A = 21
# VARIANT_INDEX_B = 20 # Optional, set to -1 if only one variant is desired

N_SEEDS = 50

std_linewidth = 3


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		plotOutDir = os.path.join(plotOutDir, plotOutFileName)

		avg_counterfactual_diff_A = {}
		avg_counterfactual_diff_B = {}

		for seed_index in range(N_SEEDS):

			# Check to make sure the seed data exists for all variants
			all_cells_control = self.ap.get_cells(
				variant=[VARIANT_INDEX_CONTROL],
				seed=[seed_index],
				generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
				only_successful=True)
			if len(all_cells_control) == 0:
				continue
			if VARIANT_INDEX_A != -1:
				all_cells_A = self.ap.get_cells(
					variant=[VARIANT_INDEX_A],
					seed=[seed_index],
					generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
					only_successful=True)
				if len(all_cells_A) == 0:
					continue
			if VARIANT_INDEX_B != -1:
				all_cells_B = self.ap.get_cells(
					variant=[VARIANT_INDEX_B],
					seed=[seed_index],
					generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
					only_successful=True)
				if len(all_cells_B) == 0:
					continue

			# Get data for control variant
			# Load generation labels
			# dt_control = read_stacked_columns(
			# 	all_cells_control, 'Main', 'time',
			# 	fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			# num_time_steps_control = read_stacked_columns(
			# 	all_cells_control, 'Main', 'time',
			# 	fun=lambda x: len(x)).squeeze()
			# gen_labels_control = np.repeat(np.arange(len(dt_control)), num_time_steps_control)
			# unique_gen_labels_control = np.unique(gen_labels_control)
			# gen_start_indexes_control = np.array(
			# 	[gen_labels_control.tolist().index(i) for i in unique_gen_labels_control])
			# gen_end_indexes_control = np.concatenate((
			# 	np.array(gen_start_indexes_control[1:] - 1), np.array([len(gen_labels_control) - 1])))
			time_control = read_stacked_columns(
				all_cells_control, 'Main', 'time', ignore_exception=True)
			initial_time_absolute_control = time_control[0]
			time_control = time_control - initial_time_absolute_control
			sim_dir_control = all_cells_control[0]
			simOutDir_control = os.path.join(sim_dir_control, 'simOut')

			# RNAP Counts
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts_control,) = read_stacked_bulk_molecules(
				all_cells_control, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts_control = TableReader(
				os.path.join(simOutDir_control, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts_control.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts_control = read_stacked_columns(
				all_cells_control, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			total_rnap_counts_control = inactive_rnap_counts_control + active_rnap_counts_control
			rel_rnap_counts_control = total_rnap_counts_control / total_rnap_counts_control[0]

			# Ribosome Counts
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s_control, complex_counts_50s_control) = read_stacked_bulk_molecules(
				all_cells_control, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts_control = np.minimum(
				complex_counts_30s_control, complex_counts_50s_control)
			# Active
			unique_molecule_counts_table_control = TableReader(
				os.path.join(simOutDir_control, "UniqueMoleculeCounts"))
			ribosome_index_control = unique_molecule_counts_table_control.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			unique_molecule_counts_table_control.close()
			active_ribosome_counts_control = read_stacked_columns(
				all_cells_control, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index_control]
			# Total
			total_ribosome_counts_control = inactive_ribosome_counts_control + active_ribosome_counts_control
			rel_ribosome_counts_control = total_ribosome_counts_control / total_ribosome_counts_control[0]


			# Get data for variant A
			# Load generation labels
			# dt_A = read_stacked_columns(
			# 	all_cells_A, 'Main', 'time',
			# 	fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			# num_time_steps_A = read_stacked_columns(
			# 	all_cells_A, 'Main', 'time',
			# 	fun=lambda x: len(x)).squeeze()
			# gen_labels_A = np.repeat(np.arange(len(dt_A)), num_time_steps_A)
			# unique_gen_labels_A = np.unique(gen_labels_A)
			# gen_start_indexes_A = np.array(
			# 	[gen_labels_A.tolist().index(i) for i in unique_gen_labels_A])
			# gen_end_indexes_A = np.concatenate((
			# 	np.array(gen_start_indexes_A[1:] - 1), np.array([len(gen_labels_A) - 1])))
			time_A = read_stacked_columns(
				all_cells_A, 'Main', 'time', ignore_exception=True)
			initial_time_absolute_A = time_A[0]
			time_A = time_A - initial_time_absolute_A
			sim_dir = all_cells_A[0]
			simOutDir = os.path.join(sim_dir, 'simOut')

			# RNAP Counts
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts_A,) = read_stacked_bulk_molecules(
				all_cells_A, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts_A = read_stacked_columns(
				all_cells_A, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			total_rnap_counts_A = inactive_rnap_counts_A + active_rnap_counts_A
			rel_rnap_counts_A = total_rnap_counts_A / total_rnap_counts_A[0]

			# Ribosome Counts
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s_A, complex_counts_50s_A) = read_stacked_bulk_molecules(
				all_cells_A, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts_A = np.minimum(
				complex_counts_30s_A, complex_counts_50s_A)
			# Active
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			unique_molecule_counts_table.close()
			active_ribosome_counts_A = read_stacked_columns(
				all_cells_A, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# Total
			total_ribosome_counts_A = inactive_ribosome_counts_A + active_ribosome_counts_A
			rel_ribosome_counts_A = total_ribosome_counts_A / total_ribosome_counts_A[0]


			# Get data for variant B
			if VARIANT_INDEX_B != -1:
				# Load generation labels
				# dt_B = read_stacked_columns(
				# 	all_cells_B, 'Main', 'time',
				# 	fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
				# num_time_steps_B = read_stacked_columns(
				# 	all_cells_B, 'Main', 'time',
				# 	fun=lambda x: len(x)).squeeze()
				# gen_labels_B = np.repeat(np.arange(len(dt_B)), num_time_steps_B)
				# unique_gen_labels_B = np.unique(gen_labels_B)
				# gen_start_indexes_B = np.array(
				# 	[gen_labels_B.tolist().index(i) for i in unique_gen_labels_B])
				# gen_end_indexes_B = np.concatenate((
				# 	np.array(gen_start_indexes_B[1:] - 1), np.array([len(gen_labels_B) - 1])))
				time_B = read_stacked_columns(
					all_cells_B, 'Main', 'time', ignore_exception=True)
				initial_time_absolute_B = time_B[0]
				time_B = time_B - initial_time_absolute_B
				sim_dir = all_cells_B[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# RNAP Counts
				# Inactive
				rnap_id = [sim_data.molecule_ids.full_RNAP]
				(inactive_rnap_counts_B,) = read_stacked_bulk_molecules(
					all_cells_B, (rnap_id,), ignore_exception=True)
				# Active
				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				active_rnap_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_RNAP')
				active_rnap_counts_B = read_stacked_columns(
					all_cells_B, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts',
					ignore_exception=True)[:, active_rnap_index]
				total_rnap_counts_B = inactive_rnap_counts_B + active_rnap_counts_B
				rel_rnap_counts_B = total_rnap_counts_B / total_rnap_counts_B[0]

				# Ribosome Counts
				# Inactive
				complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
				complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
				(complex_counts_30s_B, complex_counts_50s_B) = read_stacked_bulk_molecules(
					all_cells_B, (complex_id_30s, complex_id_50s), ignore_exception=True)
				inactive_ribosome_counts_B = np.minimum(
					complex_counts_30s_B, complex_counts_50s_B)
				# Active
				unique_molecule_counts_table = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosome_index = unique_molecule_counts_table.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')
				unique_molecule_counts_table.close()
				active_ribosome_counts_B = read_stacked_columns(
					all_cells_B, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
				# Total
				total_ribosome_counts_B = inactive_ribosome_counts_B + active_ribosome_counts_B
				rel_ribosome_counts_B = total_ribosome_counts_B / total_ribosome_counts_B[0]

			# Truncate data to the minimum length across all variants
			if VARIANT_INDEX_B != -1:
				min_length = min(
					len(time_control), len(time_A), len(time_B))
				time_B = time_B[:min_length]
				rel_rnap_counts_B = rel_rnap_counts_B[:min_length]
				rel_ribosome_counts_B = rel_ribosome_counts_B[:min_length]
			else:
				min_length = min(len(time_control), len(time_A))
			time_control = time_control[:min_length]
			rel_rnap_counts_control = rel_rnap_counts_control[:min_length]
			rel_ribosome_counts_control = rel_ribosome_counts_control[:min_length]
			time_A = time_A[:min_length]
			rel_rnap_counts_A = rel_rnap_counts_A[:min_length]
			rel_ribosome_counts_A = rel_ribosome_counts_A[:min_length]

			counterfactual_rel_rnap_counts_A = rel_rnap_counts_control - rel_rnap_counts_A
			counterfactual_rel_ribosome_counts_A = rel_ribosome_counts_control - rel_ribosome_counts_A
			counterfactual_diff_A = counterfactual_rel_rnap_counts_A - counterfactual_rel_ribosome_counts_A
			avg_counterfactual_diff_A[seed_index] = np.mean(counterfactual_diff_A)

			if VARIANT_INDEX_B != -1:
				counterfactual_rel_rnap_counts_B = rel_rnap_counts_control - rel_rnap_counts_B
				counterfactual_rel_ribosome_counts_B = rel_ribosome_counts_control - rel_ribosome_counts_B
				counterfactual_diff_B = counterfactual_rel_rnap_counts_B - counterfactual_rel_ribosome_counts_B
				avg_counterfactual_diff_B[seed_index] = np.mean(counterfactual_diff_B)

			if seed_index == 0:
				total_plots = N_SEEDS # TODO: Modularize and get rid of this magic number
				mpl.rcParams['axes.spines.right'] = False
				mpl.rcParams['axes.spines.top'] = False
				plt.rcParams['font.size'] = 20
				fig = plt.figure(figsize=(18, total_plots * 3))

				plot_num = 1
				ax1 = plt.subplot(total_plots, 1, plot_num)
				ax1.axhline(0, color="black", linewidth=std_linewidth)
				# Plot counterfactual_diff
				ax1.plot(
					time_A / 60.0, counterfactual_diff_A,
					color=VARIANT_A_COLOR, linewidth=std_linewidth,
					label=f"Var. {VARIANT_INDEX_A}")
				if VARIANT_INDEX_B != -1:
					ax1.plot(
						time_B / 60.0, counterfactual_diff_B,
						color=VARIANT_B_COLOR, linewidth=std_linewidth,
						label=f"Var. {VARIANT_INDEX_B}")
				ax1.set_ylabel("C.F. Rel. RNAP - Ribo")
				ax1.set_title(f"Seed {seed_index + 1}")
				ax1.legend()
				ax1.set_xlabel("Time (min)")

			else:
				plot_num += 1
				ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				ax.axhline(0, color="black", linewidth=std_linewidth)
				# Plot counterfactual_diff
				ax.plot(
					time_A / 60.0, counterfactual_diff_A,
					color=VARIANT_A_COLOR, linewidth=std_linewidth,
					label=f"Var. {VARIANT_INDEX_A}")
				if VARIANT_INDEX_B != -1:
					ax.plot(
						time_B / 60.0, counterfactual_diff_B,
						color=VARIANT_B_COLOR, linewidth=std_linewidth,
						label=f"Var. {VARIANT_INDEX_B}")
				ax.set_ylabel("C.F. Rel. RNAP - Ribo")
				ax.set_title(f"Seed {seed_index + 1}")
				ax.set_xlabel("Time (min)")
				ax.legend()

		# Save figure
		fig.subplots_adjust(hspace=0.5)
		if VARIANT_INDEX_B != -1:
			variants_str = f"{VARIANT_INDEX_A}_vs_{VARIANT_INDEX_B}"
		else:
			variants_str = f"{VARIANT_INDEX_A}"
		exportFigure(
			plt, plotOutDir,
			plotOutFileName + f"_var_{variants_str}",
			metadata)
		plt.close("all")

		# Make a histogram of the counterfactual differences
		fig2 = plt.figure(figsize=(12, 6))
		ax2 = fig2.add_subplot(111)
		counterfactual_diffs_A = np.array(list(avg_counterfactual_diff_A.values()))
		counterfactual_diffs_B = np.array(list(avg_counterfactual_diff_B.values()))
		ax2.hist(
			counterfactual_diffs_A.flatten(), color=VARIANT_A_COLOR,
			alpha=0.5, label=f"Var. {VARIANT_INDEX_A}")
		ax2.axvline(
			counterfactual_diffs_A.mean(), color=VARIANT_A_COLOR,
			linestyle='--', linewidth=std_linewidth, label=f"Avg. Var. {VARIANT_INDEX_A}")
		if VARIANT_INDEX_B != -1:
			ax2.hist(
				counterfactual_diffs_B.flatten(), color=VARIANT_B_COLOR,
				alpha=0.5, label=f"Var. {VARIANT_INDEX_B}")
			ax2.axvline(
				counterfactual_diffs_B.mean(), color=VARIANT_B_COLOR,
				linestyle='--', linewidth=std_linewidth, label=f"Avg. Var. {VARIANT_INDEX_B}")
		ax2.axvline(0, color="black", linewidth=std_linewidth)
		ax2.set_xlabel("C.F. Rel. RNAP - Ribo")
		ax2.set_ylabel("Count")
		ax2.set_title("Counterfactual Differences Histogram")
		ax2.legend()
		exportFigure(
			plt, plotOutDir,
			plotOutFileName + f"_var_{variants_str}_histogram",
			metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
