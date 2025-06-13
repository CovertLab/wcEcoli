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
	exportFigure, read_stacked_columns, read_stacked_bulk_molecules,
	read_bulk_molecule_counts)

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

poster_colors = {
    "light_gray": (0.75, 0.75, 0.75),
    "poster_green": (66/255, 170/255, 154/255),
    "poster_blue": (27/255, 132/255, 198/255),
    "poster_purple": (188/255, 140/255, 191/255),
    "poster_gold": (221/255, 203/255, 119/255),
    "poster_light_blue": (136/255, 205/255, 240/255),
    "poster_red": (202/255, 0/255, 32/255),
}

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		n_total_gens = self.ap.n_generation
		start_gen_index = 7
		end_gen_index = 11

		# For SHEROCK
		selected_variant_indexes = [
			# 0, # Expression factor: 0, Translation efficiency: 0, Doubling time: 53
			23, # Expression factor: 8.5, Translation efficiency: 0, Doubling time: 50
			8, # Expression factor: 8.5, Translation Efficiency: 2.20, Doubling time: 66
			1, # Expression factor: 8.5, Translation Efficiency: 3.25, Doubling time: 78
		]

		var_to_color = {
			0: poster_colors["poster_gold"],
			1: poster_colors["poster_green"],
			8: poster_colors["poster_purple"],
			23: poster_colors["poster_blue"],
		}

		# FOR LOCAL TESTING
		selected_variant_indexes = [0, 20, 21]
		var_to_color = {
			0: poster_colors["poster_gold"],
			20: poster_colors["poster_green"],
			21: poster_colors["poster_purple"],
			}

		selected_seed_indexes = [0]

		all_cells_dict = {}
		all_time_dict = {}
		all_time_no_first_dict = {}
		all_num_time_steps_dict = {}
		for seed_index in selected_seed_indexes:
			all_cells_dict[seed_index] = {}
			all_time_dict[seed_index] = {}
			all_time_no_first_dict[seed_index] = {}
			all_num_time_steps_dict[seed_index] = {}

		# Loop through variant indexes
		for i, variant_index in enumerate(selected_variant_indexes):
			for seed_index in selected_seed_indexes:
				all_cells = self.ap.get_cells(
					variant=[variant_index],
					seed=[seed_index],
					generation=np.arange(start_gen_index, end_gen_index),
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

				all_cells_dict[seed_index][variant_index] = all_cells
				all_time_dict[seed_index][variant_index] = time
				all_time_no_first_dict[seed_index][variant_index] = time_no_first
				all_num_time_steps_dict[seed_index][variant_index] = num_time_steps



		plot_types = ["counts", "ratios"]

		for plot_type in plot_types:

			# TODO: loop over seed indexes later
			seed_index = selected_seed_indexes[0]
			variants_set = all_cells_dict[seed_index].keys()
			all_cells_set = all_cells_dict[seed_index]
			time_set = all_time_dict[seed_index]
			time_no_first_set = all_time_no_first_dict[seed_index]
			num_time_steps_set = all_num_time_steps_dict[seed_index]

			# Create figure
			total_plots = 30 # TODO: Modularize and get rid of this magic number
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.figure(figsize = (12, total_plots*3))


			plot_num = 1
			ax1 = plt.subplot(total_plots, 1, plot_num)

			# Doubling Time
			for variant_index in variants_set:
				cell_paths = all_cells_set[variant_index]
				time = time_set[variant_index]
				num_time_steps = num_time_steps_set[variant_index]
				color = var_to_color[variant_index]
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
				color = var_to_color[variant_index]
				cell_mass = read_stacked_columns(
					cell_paths, 'Mass', 'cellMass',
					ignore_exception=True).squeeze()
				if plot_type == "counts":
					plt.plot(time / 60.0, cell_mass, color=color, linewidth=1,
							 label=f"Variant {variant_index}")
					plt.ylabel("Cell Mass (fg)", fontsize=8)
				elif plot_type == "ratios":
					num_time_steps = num_time_steps_set[variant_index]
					initial_cell_mass = read_stacked_columns(
						cell_paths, 'Mass', 'cellMass',
						ignore_exception=True, fun=lambda x: x[0]).squeeze()
					initial_cell_mass_vec = np.repeat(initial_cell_mass, num_time_steps)
					cell_mass_ratio = cell_mass / initial_cell_mass_vec
					plt.plot(time / 60.0, cell_mass_ratio, color=color, linewidth=1,
							 label=f"Variant {variant_index}")
					plt.axhline(y=2, color='k', linestyle='--', linewidth=0.5)
					plt.ylabel("Initial Cell Mass Ratio", fontsize=8)
			plt.xlabel("Time (min)", fontsize=8)
			plt.legend(fontsize=6, loc='upper left', ncol=2)
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








if __name__ == "__main__":
	Plot().cli()
