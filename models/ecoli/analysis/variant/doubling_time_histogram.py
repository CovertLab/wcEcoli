"""
Plot doubling time for all generations, early generations, and/or late (i.e. not early) generations
Plot percentage of simulations that reach a given generation number for each variant
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, stacked_cell_max_mask
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS, labeled_indexable_hist

exclude_timeout_cells = 1 # 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_early_gens = 1 # 1 to plot early (before MIN_LATE_CELL_INDEX), and late generationss in addition to all generations

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index
MAX_CELL_INDEX = 8 # do not include any generation >= this index
COUNT_INDEX = 7 # Count number of sims that reach this generation (remember index 7 corresponds to generation 8)
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized
MAX_CELL_LENGTH = 180
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def bar(self, ax, data, xlabel, ylabel, xlim=None):
		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			ax.bar(variant, variant_data, color=color, label=f'Var {variant}')

		if xlim:
			ax.set_xlim(xlim)
		self.remove_border(ax)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		print("Running analysis script with exclude_timeout_cells=", exclude_timeout_cells,
			  " and exclude_early_gens=",exclude_early_gens)

		# Data extraction
		print("---Data Extraction---")
		doubling_times = {}
		reached_count_gen = {}
		generations = {}

		variants = self.ap.get_variants()
		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			exclude_timeout_cell_mask = stacked_cell_max_mask(all_cells, 'Main', 'time', MAX_CELL_LENGTH,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			all_cells_gens = np.array([int(c.split("/")[-2][-6:]) for c in all_cells])[exclude_timeout_cell_mask]
			generations[variant] = all_cells_gens

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
									  fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[exclude_timeout_cell_mask]

			# Count the number of simulations that reach generation COUNT_INDEX + 1
			reached_count_gen[variant] = np.count_nonzero(all_cells_gens == COUNT_INDEX)/np.count_nonzero(all_cells_gens == 0)
			
		# Plotting
		print("---Plotting---")
		std_xlim = [30,185]
		std_dt_xlab = 'Doubling Time (min)'
		std_bar_xlab = 'Variant'
		std_bar_ylab = 'Percentage of Sims that Reached Generation ' + str(COUNT_INDEX+1)

		data_start = [0]
		data_end = [MAX_CELL_INDEX]
		plot_label = ['_all_gens']
		if exclude_early_gens:
			data_start += [0,MIN_LATE_CELL_INDEX]
			data_end += [MIN_LATE_CELL_INDEX,MAX_CELL_INDEX]
			plot_label += ['_early_gens', '_late_gens']

		for j in range(len(data_start)):
			_, axes = plt.subplots(2, 1, figsize=(10, 10))
			labeled_indexable_hist(self, axes[0], doubling_times, generations, data_start[j], data_end[j], COLORS, std_dt_xlab)
			self.bar(axes[1], reached_count_gen, std_bar_xlab, std_bar_ylab)
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName+plot_label[j], metadata)

			axes[0].set_xlim(std_xlim)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_label[j] + '_trimmed', metadata)

if __name__ == "__main__":
	Plot().cli()
