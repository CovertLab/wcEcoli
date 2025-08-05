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
from models.ecoli.analysis.variant.new_gene_narrow_multigen_schematic_poster import \
	SELECTED_SEED_INDEXES
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
END_GEN_INDEX = 8 # Not inclusive
GEN_RANGE = np.arange(START_GEN_INDEX, END_GEN_INDEX)

VARIANT_INDEX_A = 0

SELECTED_SEED_INDEXES = range(50)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		plotOutDir = os.path.join(plotOutDir, plotOutFileName)
		if not os.path.exists(plotOutDir):
			os.makedirs(plotOutDir)

		num_times_resized = 0

		for i, seed_index in enumerate(SELECTED_SEED_INDEXES):
			all_cells_A = self.ap.get_cells(
				variant=[VARIANT_INDEX_A],
				seed=[seed_index],
				generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
				only_successful=True)
			if len(all_cells_A) == 0:
				continue
			sim_dir = all_cells_A[0]
			simOutDir = os.path.join(sim_dir, 'simOut')

			time_A = read_stacked_columns(
				all_cells_A, 'Main', 'time', ignore_exception=True)
			initial_time_absolute_A = time_A[0]
			time_A = time_A - initial_time_absolute_A

			if seed_index == SELECTED_SEED_INDEXES[0]:
				# Set up data tables
				max_num_time_steps = len(time_A)

				num_time_steps_matrix = np.full((
					len(SELECTED_SEED_INDEXES), 2),
					-1.0, dtype=float)

				time_matrix = np.full((
					len(SELECTED_SEED_INDEXES), 1 + max_num_time_steps),
					-1.0, dtype=float)

				rnap_matrix = np.full((
					len(SELECTED_SEED_INDEXES), 1 + max_num_time_steps),
					-1.0, dtype=float)

				ribosome_matrix = np.full((
					len(SELECTED_SEED_INDEXES), 1 + max_num_time_steps),
					-1.0, dtype=float)

				mass_matrix = np.full((
					len(SELECTED_SEED_INDEXES), 1 + max_num_time_steps),
					-1.0, dtype=float)

			elif len(time_A) > max_num_time_steps:
				# Update data table sizes
				num_times_resized += 1
				num_extra_columns = len(time_A) - max_num_time_steps
				max_num_time_steps = len(time_A)
				extra_columns = np.full((
					len(SELECTED_SEED_INDEXES), num_extra_columns),
					-1.0, dtype=float)

				time_matrix = np.hstack((time_matrix, extra_columns))

				rnap_matrix = np.hstack((rnap_matrix, extra_columns))

				ribosome_matrix = np.hstack((ribosome_matrix, extra_columns))

				mass_matrix = np.hstack((mass_matrix, extra_columns))

			num_time_steps_matrix[i, 0] = seed_index
			num_time_steps_matrix[i, 1] = len(time_A)

			time_matrix[i, 0] = seed_index
			time_A = time_A.flatten()
			time_matrix[i, 1:(1 + len(time_A))] = time_A

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

			rnap_matrix[i, 0] = seed_index
			total_rnap_counts_A = total_rnap_counts_A.flatten()
			rnap_matrix[i, 1:(1 + len(time_A))] = total_rnap_counts_A

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

			ribosome_matrix[i, 0] = seed_index
			total_ribosome_counts_A = total_ribosome_counts_A.flatten()
			ribosome_matrix[i, 1:(1 + len(time_A))] = total_ribosome_counts_A

			# Cell Mass
			cell_mass_A = read_stacked_columns(
				all_cells_A, 'Mass', 'cellMass', ignore_exception=True).squeeze()

			mass_matrix[i, 0] = seed_index
			cell_mass_A = cell_mass_A.flatten()
			mass_matrix[i, 1:(1 + len(time_A))] = cell_mass_A


			# # Instantaneous Doubling Time
			# instantaneous_growth_rate_A = read_stacked_columns(
			# 	all_cells_A, 'Mass', 'instantaneous_growth_rate',
			# 	ignore_exception=True).squeeze()
			# instantaneous_dt_A = np.log(2) / instantaneous_growth_rate_A / 60.0

			# # Disable line clipping
			# for ax in fig.get_axes():
			# 	for artist in ax.get_children():
			# 		try:
			# 			artist.set_clip_on(False)
			# 		except AttributeError:
			# 			pass
			# # Save figure
			# fig.subplots_adjust(hspace=0.5)
			# print(f"\nSeed: {seed_index}")
			# print("Total number of plots made: ", plot_num - 1)
			# # plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
			# if VARIANT_INDEX_B != -1:
			# 	variants_str = f"{VARIANT_INDEX_A}_vs_{VARIANT_INDEX_B}"
			# else:
			# 	variants_str = f"{VARIANT_INDEX_A}"
			# exportFigure(
			# 	plt, plotOutDir,
			# 	plotOutFileName + f"_seed_{seed_index}_var_{variants_str}",
			# 	metadata)
			# plt.close("all")

		# Save number of time steps
		num_time_steps_matrix_file = os.path.join(
			plotOutDir,
			f"variant_{VARIANT_INDEX_A}_num_time_steps.csv")
		np.savetxt(
			num_time_steps_matrix_file, num_time_steps_matrix, delimiter=",",
			header="seed_index,num_time_steps",
			comments="", fmt="%s")

		# Save time matrix
		time_matrix_file = os.path.join(
			plotOutDir,
			f"variant_{VARIANT_INDEX_A}_time_matrix.csv")
		np.savetxt(
			time_matrix_file, time_matrix, delimiter=",",
			header="seed_index," + ",".join([f"time_step_{i}" for i in range(max_num_time_steps)]),
			comments="", fmt="%s")

		# Save RNAP matrix
		rnap_matrix_file = os.path.join(
			plotOutDir,
			f"variant_{VARIANT_INDEX_A}_rnap_matrix.csv")
		np.savetxt(
			rnap_matrix_file, rnap_matrix, delimiter=",",
			header="seed_index," + ",".join([f"time_step_{i}" for i in range(max_num_time_steps)]),
			comments="", fmt="%s")

		# Save ribosome matrix
		ribosome_matrix_file = os.path.join(
			plotOutDir,
			f"variant_{VARIANT_INDEX_A}_ribosome_matrix.csv")
		np.savetxt(
			ribosome_matrix_file, ribosome_matrix, delimiter=",",
			header="seed_index," + ",".join([f"time_step_{i}" for i in range(max_num_time_steps)]),
			comments="", fmt="%s")

		# Save mass matrix
		mass_matrix_file = os.path.join(
			plotOutDir,
			f"variant_{VARIANT_INDEX_A}_mass_matrix.csv")
		np.savetxt(
			mass_matrix_file, mass_matrix, delimiter=",",
			header="seed_index," + ",".join([f"time_step_{i}" for i in range(max_num_time_steps)]),
			comments="", fmt="%s")

		print(f"Number of times resized: {num_times_resized}")


if __name__ == "__main__":
	Plot().cli()
