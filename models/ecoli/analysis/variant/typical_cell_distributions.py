"""
Analyze how production machinery counts change over time and across at most two
variants.
"""

### TODO: subdivide by cell cycle fractions instead of time steps

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

START_GEN_INDEX = 8
END_GEN_INDEX = 24 # Not inclusive
GEN_RANGE = np.arange(START_GEN_INDEX, END_GEN_INDEX)

VARIANT_INDEX_A = 0

SELECTED_SEED_INDEXES = range(50)

NUM_STEPS_PER_CELL_CYCLE = 500

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		plotOutDir = os.path.join(plotOutDir, plotOutFileName)
		if not os.path.exists(plotOutDir):
			os.makedirs(plotOutDir)

		# If the data files do not exist, generate them
		if not os.path.exists(os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_num_time_steps.csv")):
			print("Making matrices for variant", VARIANT_INDEX_A)
			total_cells_included = 0

			for i, seed_index in enumerate(SELECTED_SEED_INDEXES):
				for j, GEN_INDEX in enumerate(GEN_RANGE):
					all_cells_A = self.ap.get_cells(
						variant=[VARIANT_INDEX_A],
						seed=[seed_index],
						generation=np.arange(GEN_INDEX,GEN_INDEX + 1),
						only_successful=True)
					if len(all_cells_A) == 0:
						continue
					sim_dir = all_cells_A[0]
					simOutDir = os.path.join(sim_dir, 'simOut')

					time_A = read_stacked_columns(
						all_cells_A, 'Main', 'time', ignore_exception=True)
					initial_time_absolute_A = time_A[0]
					time_A = time_A - initial_time_absolute_A

					num_time_points = min(len(time_A), NUM_STEPS_PER_CELL_CYCLE)
					time_indices_to_include = np.linspace(0, len(time_A) - 1, num_time_points, dtype=int)

					if total_cells_included == 0:
						# Set up data tables
						num_time_steps_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 2),
							-1.0, dtype=float)

						time_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)

						rnap_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)

						ribosome_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)

						mass_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)

						unnormalized_rnap_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)

						unnormalized_ribosome_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)
						
						unnormalized_mass_matrix = np.full((
							len(SELECTED_SEED_INDEXES) * len(GEN_RANGE), 1 + NUM_STEPS_PER_CELL_CYCLE),
							-1.0, dtype=float)

					num_time_steps_matrix[total_cells_included, 0] = seed_index
					num_time_steps_matrix[total_cells_included, 1] = len(time_A)

					time_matrix[total_cells_included, 0] = seed_index
					time_A = time_A.flatten()
					time_matrix[total_cells_included, 1:(1 + NUM_STEPS_PER_CELL_CYCLE)] = time_A[time_indices_to_include]

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
					unnormalized_rnap_counts_A = total_rnap_counts_A.copy()
					# Normalize to the first time point
					total_rnap_counts_A = total_rnap_counts_A / total_rnap_counts_A[0]

					rnap_matrix[total_cells_included, 0] = seed_index
					total_rnap_counts_A = total_rnap_counts_A.flatten()
					rnap_matrix[total_cells_included, 1:(1 +
						NUM_STEPS_PER_CELL_CYCLE)] = total_rnap_counts_A[time_indices_to_include]

					unnormalized_rnap_matrix[total_cells_included, 0] = seed_index
					unnormalized_rnap_counts_A = unnormalized_rnap_counts_A.flatten()
					unnormalized_rnap_matrix[total_cells_included, 1:(1 +
						NUM_STEPS_PER_CELL_CYCLE)] = unnormalized_rnap_counts_A[time_indices_to_include]

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
					unnormalized_ribosome_counts_A = total_ribosome_counts_A.copy()
					# Normalize to the first time point
					total_ribosome_counts_A = total_ribosome_counts_A / total_ribosome_counts_A[0]

					ribosome_matrix[total_cells_included, 0] = seed_index
					total_ribosome_counts_A = total_ribosome_counts_A.flatten()
					ribosome_matrix[total_cells_included, 1:(1 + NUM_STEPS_PER_CELL_CYCLE)] = total_ribosome_counts_A[time_indices_to_include]

					unnormalized_ribosome_matrix[total_cells_included, 0] = seed_index
					unnormalized_ribosome_counts_A = unnormalized_ribosome_counts_A.flatten()
					unnormalized_ribosome_matrix[total_cells_included, 1:(1 +
						NUM_STEPS_PER_CELL_CYCLE)] = unnormalized_ribosome_counts_A[time_indices_to_include]

					# Cell Mass
					cell_mass_A = read_stacked_columns(
						all_cells_A, 'Mass', 'cellMass', ignore_exception=True).squeeze()
					unnormalized_cell_mass_A = cell_mass_A.copy()
					# Normalize to the first time point
					cell_mass_A = cell_mass_A / cell_mass_A[0]

					mass_matrix[total_cells_included, 0] = seed_index
					cell_mass_A = cell_mass_A.flatten()
					mass_matrix[total_cells_included, 1:(1 +
						NUM_STEPS_PER_CELL_CYCLE)] = cell_mass_A[time_indices_to_include]
					
					unnormalized_mass_matrix[total_cells_included, 0] = seed_index
					unnormalized_cell_mass_A = unnormalized_cell_mass_A.flatten()
					unnormalized_mass_matrix[total_cells_included, 1:(1 +
						NUM_STEPS_PER_CELL_CYCLE)] = unnormalized_cell_mass_A[time_indices_to_include]

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
					total_cells_included += 1

			# Save number of time steps
			num_time_steps_matrix = num_time_steps_matrix[:total_cells_included, :]
			num_time_steps_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_num_time_steps.csv")
			np.savetxt(
				num_time_steps_matrix_file, num_time_steps_matrix, delimiter=",",
				header="seed_index,num_time_steps",
				comments="", fmt="%s")

			# Save time matrix
			time_matrix = time_matrix[:total_cells_included, :]
			time_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_time_matrix.csv")
			np.savetxt(
				time_matrix_file, time_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

			# Save RNAP matrix
			rnap_matrix = rnap_matrix[:total_cells_included, :]
			rnap_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_rnap_matrix.csv")
			np.savetxt(
				rnap_matrix_file, rnap_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

			# Save unnormalized RNAP matrix
			unnormalized_rnap_matrix = unnormalized_rnap_matrix[:total_cells_included, :]
			unnormalized_rnap_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_unnormalized_rnap_matrix.csv")
			np.savetxt(
				unnormalized_rnap_matrix_file, unnormalized_rnap_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

			# Save ribosome matrix
			ribosome_matrix = ribosome_matrix[:total_cells_included, :]
			ribosome_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_ribosome_matrix.csv")
			np.savetxt(
				ribosome_matrix_file, ribosome_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

			# Save unnormalized ribosome matrix
			unnormalized_ribosome_matrix = unnormalized_ribosome_matrix[:total_cells_included, :]
			unnormalized_ribosome_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_unnormalized_ribosome_matrix.csv")
			np.savetxt(
				unnormalized_ribosome_matrix_file, unnormalized_ribosome_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

			# Save mass matrix
			mass_matrix = mass_matrix[:total_cells_included, :]
			mass_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_mass_matrix.csv")
			np.savetxt(
				mass_matrix_file, mass_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

			print(f"Total number of cells included: {total_cells_included}")
			print(f"Expected: {len(SELECTED_SEED_INDEXES) * len(GEN_RANGE)}")

			# Save unnormalized mass matrix
			unnormalized_mass_matrix = unnormalized_mass_matrix[:total_cells_included, :]
			unnormalized_mass_matrix_file = os.path.join(
				plotOutDir,
				f"variant_{VARIANT_INDEX_A}_unnormalized_mass_matrix.csv")
			np.savetxt(
				unnormalized_mass_matrix_file, unnormalized_mass_matrix, delimiter=",",
				header="seed_index," + ",".join([
					f"cell_cycle_step_{i}_of_{NUM_STEPS_PER_CELL_CYCLE}"
					for i in range(NUM_STEPS_PER_CELL_CYCLE)]),
				comments="", fmt="%s")

		else:
			# Load in the matrices
			num_time_steps_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_num_time_steps.csv"),
				delimiter=",", skiprows=1)
			time_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_time_matrix.csv"),
				delimiter=",", skiprows=1)
			rnap_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_rnap_matrix.csv"),
				delimiter=",", skiprows=1)
			unnormalized_rnap_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_unnormalized_rnap_matrix.csv"),
				delimiter=",", skiprows=1)
			ribosome_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_ribosome_matrix.csv"),
				delimiter=",", skiprows=1)
			unnormalized_ribosome_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_unnormalized_ribosome_matrix.csv"),
				delimiter=",", skiprows=1)
			mass_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_mass_matrix.csv"),
				delimiter=",", skiprows=1)
			unnormalized_mass_matrix = np.loadtxt(
				os.path.join(plotOutDir, f"variant_{VARIANT_INDEX_A}_unnormalized_mass_matrix.csv"),
				delimiter=",", skiprows=1)

		cell_cycle_steps_to_plot_distrib = np.arange(0, NUM_STEPS_PER_CELL_CYCLE, 100)
		cell_cycle_steps_to_plot_distrib = np.append(
			cell_cycle_steps_to_plot_distrib, NUM_STEPS_PER_CELL_CYCLE - 1)
		num_cells = num_time_steps_matrix.shape[0]

		# Plot the distributions
		for cell_cycle_step in cell_cycle_steps_to_plot_distrib:

			# Get the RNAP counts for this cell cycle step
			rnap_counts = rnap_matrix[:, 1 + cell_cycle_step]
			unnormalized_rnap_counts = unnormalized_rnap_matrix[:, 1 + cell_cycle_step]

			# Get the ribosome counts for this cell cycle step
			ribosome_counts = ribosome_matrix[:, 1 + cell_cycle_step]
			unnormalized_ribosome_counts = unnormalized_ribosome_matrix[:, 1 + cell_cycle_step]

			# Get the mass for this cell cycle step
			cell_mass = mass_matrix[:, 1 + cell_cycle_step]
			unnormalized_cell_mass = unnormalized_mass_matrix[:, 1 + cell_cycle_step]

			# Get the number of time steps
			num_time_steps = num_time_steps_matrix[:, 1]

			# Plotting
			n_bins = 20
			plt.figure(figsize=(12, 8))

			plt.subplot(2, 2, 1)
			plt.hist(rnap_counts, bins=n_bins, color=poster_colors["poster_green"], alpha=0.7)
			plt.title(f"RNAP Counts at Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE}")
			plt.xlabel("RNAP Counts")
			plt.ylabel("Frequency")
			# Add a normal distribution fit
			mean_rnap = np.mean(rnap_counts)
			std_rnap = np.std(rnap_counts)
			x = np.linspace(mean_rnap - 3*std_rnap, mean_rnap + 3*std_rnap, 100)
			y = (1 / (std_rnap * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean_rnap) / std_rnap) ** 2)
			plt.plot(x, y * len(rnap_counts) * (max(rnap_counts) - min(rnap_counts)) / 50, color='black', linestyle='--', label='Normal Fit')
			plt.legend()

			plt.subplot(2, 2, 2)
			plt.hist(ribosome_counts, bins=n_bins, color=poster_colors["poster_blue"], alpha=0.7)
			plt.title(f"Ribosome Counts at Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE}")
			plt.xlabel("Ribosome Counts")
			plt.ylabel("Frequency")
			# Add a normal distribution fit
			mean_ribosome = np.mean(ribosome_counts)
			std_ribosome = np.std(ribosome_counts)
			x = np.linspace(mean_ribosome - 3*std_ribosome, mean_ribosome + 3*std_ribosome, 100)
			y = (1 / (std_ribosome * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean_ribosome) / std_ribosome) ** 2)
			plt.plot(x, y * len(ribosome_counts) * (max(ribosome_counts) - min(ribosome_counts)) / 50, color='black', linestyle='--', label='Normal Fit')
			plt.legend()

			plt.subplot(2, 2, 3)
			plt.hist(cell_mass, bins=n_bins, color=poster_colors["poster_purple"], alpha=0.7)
			plt.title(f"Cell Mass at Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE}")
			plt.xlabel("Cell Mass (fg)")
			plt.ylabel("Frequency")
			# Add a normal distribution fit
			mean_mass = np.mean(cell_mass)
			std_mass = np.std(cell_mass)
			x = np.linspace(mean_mass - 3*std_mass, mean_mass + 3*std_mass, 100)
			y = (1 / (std_mass * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean_mass) / std_mass) ** 2)
			plt.plot(x, y * len(cell_mass) * (max(cell_mass) - min(cell_mass)) / 50, color='black', linestyle='--', label='Normal Fit')
			plt.legend()

			plt.subplot(2, 2, 4)
			plt.hist(num_time_steps_matrix[:, 1], bins=n_bins, color=poster_colors["poster_gold"], alpha=0.7)
			plt.title(f"Time Steps Per Cell Distribution, {num_cells} Cells")
			plt.xlabel("Number of Time Steps")
			plt.ylabel("Frequency")
			# Add a normal distribution fit
			mean_time_steps = np.mean(num_time_steps_matrix[:, 1])
			std_time_steps = np.std(num_time_steps_matrix[:, 1])
			x = np.linspace(mean_time_steps - 3*std_time_steps, mean_time_steps + 3*std_time_steps, 100)
			y = (1 / (std_time_steps * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean_time_steps) / std_time_steps) ** 2)
			plt.plot(x, y * len(num_time_steps_matrix[:, 1]) * (max(num_time_steps_matrix[:, 1]) - min(num_time_steps_matrix[:, 1])) / 50, color='black', linestyle='--', label='Normal Fit')
			plt.legend()

			plt.tight_layout()
			exportFigure(
				plt, plotOutDir,
				f"variant_{VARIANT_INDEX_A}_distributions_cell_cycle_step_{cell_cycle_step}",
				metadata)
			plt.close()

			# Make qq plots to show how good the normal fits are
			import scipy.stats as stats
			from scipy.stats import shapiro
			plt.figure(figsize=(12, 8))

			plt.subplot(2, 2, 1)
			stats.probplot(rnap_counts, dist="norm", plot=plt)
			rnap_stat, rnap_p = shapiro(rnap_counts)
			plt.title(f"RNAP, Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE}, S-W stat {rnap_stat} p-value {rnap_p:.4f}")
			plt.xlabel("Theoretical Quantiles")
			plt.ylabel("Sample Quantiles")
			plt.title(plt.gca().get_title(), fontsize=8)

			plt.subplot(2, 2, 2)
			stats.probplot(ribosome_counts, dist="norm", plot=plt)
			ribosome_stat, ribosome_p = shapiro(ribosome_counts)
			plt.title(f"Ribo, Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE}, S-W stat {ribosome_stat} p-value {ribosome_p:.4f}")
			plt.xlabel("Theoretical Quantiles")
			plt.ylabel("Sample Quantiles")
			plt.title(plt.gca().get_title(), fontsize=8)

			plt.subplot(2, 2, 3)
			stats.probplot(cell_mass, dist="norm", plot=plt)
			mass_stat, mass_p = shapiro(cell_mass)
			plt.title(f"Mass, Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE}, S-W stat {mass_stat} p-value {mass_p:.4f}")
			plt.xlabel("Theoretical Quantiles")
			plt.ylabel("Sample Quantiles")
			plt.title(plt.gca().get_title(), fontsize=8)

			plt.tight_layout()
			exportFigure(
				plt, plotOutDir,
				f"variant_{VARIANT_INDEX_A}_qq_plots_cell_cycle_step_{cell_cycle_step}",
				metadata)
			plt.close()

			# Plot a scatterplot of RNAP counts vs time steps / 60 and ribosome counts vs time steps / 60
			plt.figure(figsize=(12, 6))
			plt.subplot(1, 3, 1)
			plt.scatter(num_time_steps / 60, unnormalized_rnap_counts, color=poster_colors["poster_green"], alpha=0.7)
			plt.title(f"RNAP Counts vs N Time Steps / 60 (Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE})")
			plt.xlabel("Time (min)")
			plt.ylabel("RNAP Counts")
			from scipy.stats import linregress
			slope, intercept, r_value, p_value, std_err = linregress(num_time_steps / 60, unnormalized_rnap_counts)
			plt.plot(num_time_steps / 60, intercept + slope * (num_time_steps / 60), color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}, R²={r_value**2:.2f}')
			plt.legend()
			plt.title(plt.gca().get_title(), fontsize=8)

			plt.subplot(1, 3, 2)
			plt.scatter(num_time_steps / 60, unnormalized_ribosome_counts, color=poster_colors["poster_blue"], alpha=0.7)
			plt.title(f"Ribosome Counts vs N Time Steps / 60 (Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE})")
			plt.xlabel("Time (min)")
			plt.ylabel("Ribosome Counts")
			slope, intercept, r_value, p_value, std_err = linregress(num_time_steps / 60, unnormalized_ribosome_counts)
			plt.plot(num_time_steps / 60, intercept + slope * (num_time_steps / 60), color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}, R²={r_value**2:.2f}')
			plt.legend()
			plt.title(plt.gca().get_title(), fontsize=8)

			plt.subplot(1, 3, 3)
			plt.scatter(num_time_steps / 60, unnormalized_cell_mass, color=poster_colors["poster_purple"], alpha=0.7)
			plt.title(f"Cell Mass vs N Time Steps / 60 (Cell Cycle Step {cell_cycle_step} of {NUM_STEPS_PER_CELL_CYCLE})")
			plt.xlabel("Time (min)")
			plt.ylabel("Cell Mass (fg)")
			slope, intercept, r_value, p_value, std_err = linregress(num_time_steps / 60, unnormalized_cell_mass)
			plt.plot(num_time_steps / 60, intercept + slope * (num_time_steps / 60), color='red', label=f'Fit: y={slope:.2f}x+{intercept:.2f}, R²={r_value**2:.2f}')
			plt.legend()
			plt.title(plt.gca().get_title(), fontsize=8)

			plt.tight_layout()
			exportFigure(
				plt, plotOutDir,
				f"variant_{VARIANT_INDEX_A}_scatter_rnap_ribo_cell_cycle_step_{cell_cycle_step}",
				metadata)
			plt.close()


if __name__ == "__main__":
	Plot().cli()
