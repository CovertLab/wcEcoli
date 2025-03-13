"""
Plots the doubling times for all seeds and all generations (after some
threshold number of burn-in generations) as a histogram.

Notes
-----
More strictly speaking, this is the division time or cell-cycle time, for which
are only doubling in an average sense.

There are many hard-coded values in here for the bounds and bins.  This is to
standardize the output across sets of simulations.

"""
# TODO: update the docstring above

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)

from models.ecoli.analysis import cohortAnalysisPlot

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		if "init_sims" in metadata:
			num_seeds = metadata['init_sims']
		else:
			num_seeds = metadata['total_init_sims']

		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False

		plt.figure(figsize=(8.5, 10))
		total_plots = num_seeds
		plot_num = 1

		all_dts_dict = {}
		all_dts = np.array([])
		all_avg_dts = np.array([])

		all_dts_last_8 = np.array([])
		all_avg_dts_last_8 = np.array([])

		for seed in range(num_seeds):
			print("seed: ", seed)

			cell_paths = self.ap.get_cells(seed=[seed])

			if seed == 0:
				ax1 = plt.subplot(total_plots, 1, plot_num)
			else:
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			time = read_stacked_columns(
				cell_paths, 'Main', 'time', ignore_exception=True)
			dt = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			num_time_steps = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: len(x)).squeeze()

			# Create a new numpy array where each dt[i] is repeated num_time_steps[i] times
			dt_to_plot = np.repeat(dt, num_time_steps)
			plt.plot(time / 60., dt_to_plot)

			# Add a line for average doubling time
			avg_dt = np.mean(dt)
			plt.axhline(y=avg_dt, color='r', linestyle='--')

			plt.xlabel("Time (min)")
			plt.ylabel("Doubling Time (min)", fontsize="small")
			plt.title("Doubling Time")

			all_dts_dict[seed] = dt
			all_dts = np.append(all_dts, dt)
			all_avg_dts = np.append(all_avg_dts, avg_dt)

			all_dts_last_8 = np.append(all_dts_last_8, dt[-8:])
			all_avg_dts_last_8 = np.append(all_avg_dts_last_8, np.mean(dt[-8:]))

			plot_num += 1

		plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Make a boxplot of the dts for each seed
		plt.figure(figsize=(8.5, 10))
		plt.boxplot(all_dts_dict.values())
		plt.xlabel("Seed")
		plt.ylabel("Doubling Time (min)", fontsize="small")
		plt.title("Doubling Time")
		exportFigure(plt, plotOutDir, plotOutFileName + "_boxplot", metadata)
		plt.close("all")

		# Make a histogram of the dts for all seeds
		plt.figure(figsize=(8.5, 10))
		plt.hist(all_dts, bins=50, range=(0, 100))
		avg_dt = np.mean(all_dts)
		plt.axvline(x=avg_dt, color='r', linestyle='--')
		plt.xlabel("Doubling Time (min)", fontsize="small")
		plt.ylabel("Frequency", fontsize="small")
		plt.title("Doubling Time")
		exportFigure(plt, plotOutDir, plotOutFileName + "_hist", metadata)
		plt.close("all")

		# Make a scatterplot of seed vs average for that seed for all seeds
		plt.figure(figsize=(8.5, 10))
		plt.scatter(range(num_seeds), all_avg_dts)
		plt.axhline(y=avg_dt, color='r', linestyle='--')
		plt.xlabel("Seed")
		plt.ylabel("Average Doubling Time (min)", fontsize="small")
		plt.title("Doubling Time")
		exportFigure(plt, plotOutDir, plotOutFileName + "_scatter", metadata)

		# Make a histogram where one color shows the distribution of the first 4 values in each seed and the other shows the distribution of the last 4 values in each seed
		plt.figure(figsize=(8.5, 10))
		first_four = np.array([all_dts_dict[seed][:4] for seed in all_dts_dict])
		last_four = np.array([all_dts_dict[seed][-4:] for seed in all_dts_dict])
		plt.hist(first_four.flatten(), bins=50, range=(0, 100), alpha=0.5, color='b', label='First 4')
		plt.hist(last_four.flatten(), bins=50, range=(0, 100), alpha=0.5, color='r', label='Last 4')

		avg_first = np.mean(first_four.flatten())
		avg_last = np.mean(last_four.flatten())
		plt.axvline(x=avg_first, color='b', linestyle='--')
		plt.axvline(x=avg_last, color='r', linestyle='--')

		plt.xlabel("Doubling Time (min)", fontsize="small")
		plt.ylabel("Frequency", fontsize="small")
		plt.title("Doubling Time")
		plt.legend()
		exportFigure(plt, plotOutDir, plotOutFileName + "_hist_first_last", metadata)
		plt.close("all")

		# Make a scatterplot of seed vs dt for the last 8 generations for all seeds
		plt.figure(figsize=(8.5, 10))
		plt.scatter(range(num_seeds), all_avg_dts_last_8)
		plt.axhline(y=np.mean(all_avg_dts_last_8), color='r', linestyle='--')
		plt.xlabel("Seed")
		plt.ylabel("Average Doubling Time (min)", fontsize="small")
		plt.title("Doubling Time")
		exportFigure(plt, plotOutDir, plotOutFileName + "_scatter_last_8", metadata)

		# Make a histogram of the dts for all seeds for the last 8 generations
		plt.figure(figsize=(8.5, 10))
		plt.hist(all_dts_last_8, bins=25, range=(35, 75))
		avg_dt_last_8 = np.mean(all_dts_last_8)
		plt.axvline(x=avg_dt_last_8, color='r', linestyle='--')
		plt.xlabel("Doubling Time (min)", fontsize="small")
		plt.ylabel("Frequency", fontsize="small")
		plt.title("Doubling Time")
		exportFigure(plt, plotOutDir, plotOutFileName + "_hist_last_8", metadata)

		# Make a histogram of the dts for all seeds for the last 8 generations, but add a line for the average of the last 8 generations for each seed
		plt.figure(figsize=(8.5, 8.5))
		plt.hist(all_dts_last_8, bins=25, range=(35, 75))
		for seed in all_dts_dict:
			plt.axvline(x=np.mean(all_dts_dict[seed][-8:]), color='black', linestyle='--', label='Seed Average (8 cells)')
		plt.axvline(x=avg_dt_last_8, color='r', linestyle='--', label="Overall Average (128 Cells)")
		plt.xlabel("Doubling Time (min)", fontsize="small")
		plt.ylabel("Frequency", fontsize="small")
		plt.title("Doubling Time")
		# Get rid of duplicate labels in legend
		handles, labels = plt.gca().get_legend_handles_labels()
		by_label = dict(zip(labels, handles))
		plt.legend(by_label.values(), by_label.keys())
		exportFigure(plt, plotOutDir, plotOutFileName + "_hist_last_8_with_avg", metadata)
		plt.close("all")

		# Make a histogram of the dts for all seeds for the last 8 generations, but add a line for the average of the last 8 generations for each seed
		plt.figure(figsize=(8.5, 10))
		for seed in all_dts_dict:
			plt.hist(all_dts_dict[seed][-8:], bins=25, range=(35, 75), alpha=0.5)
			plt.axvline(x=np.mean(all_dts_dict[seed][-8:]), color='r', linestyle='--')
		plt.xlabel("Doubling Time (min)", fontsize="small")
		plt.ylabel("Frequency", fontsize="small")
		plt.title("Doubling Time")
		exportFigure(plt, plotOutDir, plotOutFileName + "_sep_hist_last_8_with_avg", metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
