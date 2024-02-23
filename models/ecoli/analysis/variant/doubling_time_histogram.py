"""
Plot histograms of doubling times and compares the distributions across
different variants. Also optionally plots percentage of simulations that reach
the simultation-specified number of generations without failing.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.analysis.plotting_tools import (
	DEFAULT_MATPLOTLIB_COLORS as COLORS)

FONT_SIZE = 9
DOUBLING_TIME_BOUNDS_MINUTES = [0, 180]
N_BINS = 36

# Set True to exclude cells that hit time limit
EXCLUDE_TIMEOUT_CELLS = True

# Set True to plot completion rate
PLOT_COMPLETION_RATES = True

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 4

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		max_cell_length = 180  # mins
		if not EXCLUDE_TIMEOUT_CELLS:
			max_cell_length += 1000000  # Arbitrary large number

		# Data extraction
		doubling_times = {}
		completion_rates = {}
		n_total_gens = self.ap.n_generation
		variant_indexes = self.ap.get_variants()

		# Loop through all variant indexes
		for variant_index in variant_indexes:
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Get mask to exclude cells that timed out
			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', max_cell_length,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()

			# Get doubling times from cells with this variant index
			dt = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant_index] = dt[exclude_timeout_cell_mask]

			# Count the number of simulations that reach the last generation
			num_last_gen = len(
				self.ap.get_cells(
					variant=[variant_index], generation=[n_total_gens - 1],
					only_successful=True)
				)

			# Count the number of simulations at generation zero
			num_zero_gen = len(
				self.ap.get_cells(
					variant=[variant_index], generation=[0],
					only_successful=True)
				)

			# Calculate completion rate
			completion_rates[variant_index] = num_last_gen / num_zero_gen

		n_variants = len(doubling_times)

		fig = plt.figure(figsize=(8, 2*(n_variants + 1)))
		bins = np.linspace(
			DOUBLING_TIME_BOUNDS_MINUTES[0],
			DOUBLING_TIME_BOUNDS_MINUTES[1],
			N_BINS + 1
			)

		# First subplot overlays the histograms of all variants
		ax0 = fig.add_subplot(n_variants + 1, 1, 1)
		subplots = []

		for i, (variant_index, dt) in enumerate(doubling_times.items()):
			color = COLORS[i % len(COLORS)]
			ax0.hist(
				dt, bins=bins, color=color, alpha=0.5,
				label=f'Var {variant_index} (n={len(dt)}, {np.mean(dt):.1f} $\pm$ {np.std(dt):.1f})')

			# Add vertical line at mean doubling time
			ax0.axvline(np.mean(dt), color=color, ls='--', lw=3, alpha=0.5)

			# Later subplots show the histograms for each variant separately
			ax = fig.add_subplot(n_variants + 1, 1, i + 2, sharex=ax0)
			ax.hist(
				dt, bins=bins, color=color, alpha=0.5,
				)
			ax.tick_params(labelsize=FONT_SIZE)
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			subplots.append(ax)

		# Set ylims of all subplots to be equal to the first subplot
		ax0_ylim = ax0.get_ylim()
		for ax in subplots:
			ax.set_ylim(ax0_ylim)

		ax0.legend()

		ax0.set_xlim(*DOUBLING_TIME_BOUNDS_MINUTES)
		ax0.set_xticks(
			np.linspace(
				DOUBLING_TIME_BOUNDS_MINUTES[0],
				DOUBLING_TIME_BOUNDS_MINUTES[1],
				10
				)
			)
		ax0.tick_params(labelsize=FONT_SIZE)

		ax0.set_xlabel('Doubling time (min)', fontsize=FONT_SIZE)
		ax0.spines["top"].set_visible(False)
		ax0.spines["right"].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Optionally plot bar plots for completion rates
		if PLOT_COMPLETION_RATES:
			fig = plt.figure(figsize=(8, 3))
			ax = fig.add_subplot(1, 1, 1)

			for i, (variant_index, completion_rate) in enumerate(
					completion_rates.items()):
				color = COLORS[i % len(COLORS)]
				ax.bar(variant_index, completion_rate,
						color=color, alpha=0.5, label=f'Var {variant_index}')

			self.remove_border(ax)
			ax.set_xlabel('Variant indexes', fontsize=FONT_SIZE)
			ax.set_ylabel('Completion rates', fontsize=FONT_SIZE)
			ax.tick_params(labelsize=FONT_SIZE)
			ax.legend()

			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + '_completion_rates', metadata)
			plt.close('all')

if __name__ == "__main__":
	Plot().cli()
