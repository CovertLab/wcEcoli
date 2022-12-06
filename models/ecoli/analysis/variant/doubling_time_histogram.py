"""
Plot doubling time for all generations, early generations, and late (i.e. not early) generations
Plot percentage of simulations that reach a given generation number for each variant
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS


FONT_SIZE=9
MAX_CELL_LENGTH = 180
#MAX_CELL_LENGTH += 1 # comment out this line to filter sims that reach the max time of 180 min
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized
COUNT_GENERATION = 7 # count number of sims that reached this generation per variant

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def hist(self, ax, data, xlabel, bin_width=1., xlim=None, sf=1):
		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			bins = max(1, int(np.ceil((variant_data.max() - variant_data.min()) / bin_width)))
			mean = variant_data.mean()
			std = variant_data.std()
			ax.hist(variant_data, bins, color=color, alpha=0.5,
				label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
			ax.axvline(mean, color=color, linestyle='--', linewidth=1)

		if xlim:
			ax.set_xlim(xlim)
		self.remove_border(ax)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def bar(self, ax, data, xlabel, ylabel, bin_width=1., xlim=None, sf=1):
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
		doubling_times = {}
		doubling_times_early_gens = {}
		doubling_times_late_gens = {}

		reached_count_gen = {}

		print("---Data Extraction---")
		for variant in self.ap.get_variants():
			print("Variant: ", variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt < MAX_CELL_LENGTH]

			all_cells_gens = [int(c.split("/")[-2][-6:]) for c in all_cells]
			reached_count_gen[variant] = all_cells_gens.count(COUNT_GENERATION)/all_cells_gens.count(0)

			if len(all_cells) >= MIN_LATE_CELL_INDEX:
				early_cell_index = [i for i,v in enumerate(all_cells_gens) if v < MIN_LATE_CELL_INDEX]
				late_cell_index = [i for i,v in enumerate(all_cells_gens) if v >= MIN_LATE_CELL_INDEX]

				dt_early_cells = dt[early_cell_index]
				dt_late_cells = dt[late_cell_index]

				doubling_times_early_gens[variant] = dt_early_cells[dt_early_cells < MAX_CELL_LENGTH ]
				doubling_times_late_gens[variant] = dt_late_cells[dt_late_cells < MAX_CELL_LENGTH ]


		print("---Plotting---")
		# LATE GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], doubling_times_late_gens, 'Doubling Time (min)')
		self.bar(axes[1], reached_count_gen, 'Variant','Percentage of Sims that Reached Generation ' + str(COUNT_GENERATION+1))
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_lategens', metadata)

		axes[0].set_xlim([30, 185])
		exportFigure(plt, plotOutDir, plotOutFileName + '_lategens_trimmed', metadata)

		# EARLY GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], doubling_times_early_gens, 'Doubling Time (min)')
		self.bar(axes[1], reached_count_gen, 'Variant','Percentage of Sims that Reached Generation ' + str(COUNT_GENERATION+1))
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_earlygens', metadata)

		axes[0].set_xlim([30, 185])
		exportFigure(plt, plotOutDir, plotOutFileName + '_earlygens_trimmed', metadata)

		# ALL GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], doubling_times, 'Doubling Time (min)')
		self.bar(axes[1], reached_count_gen, 'Variant','Percentage of Sims that Reached Generation ' + str(COUNT_GENERATION+1))
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_allgens', metadata)

		axes[0].set_xlim([30, 185])
		exportFigure(plt, plotOutDir, plotOutFileName + '_allgens_trimmed', metadata)


if __name__ == "__main__":
	Plot().cli()
