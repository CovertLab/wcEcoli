"""
Plot doubling time and instantaneous growth rate for all generations, early generations, and late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS


FONT_SIZE=9
MAX_CELL_LENGTH = 180
MAX_CELL_LENGTH += 1 # comment out this line to filter sims that reach the max time of 180 min
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized


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

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		doubling_times = {}
		doubling_times_early_gens = {}
		doubling_times_late_gens = {}

		growth_rates = {}
		growth_rates_early_gens = {}
		growth_rates_late_gens = {}

		def downsample(x):
			"""Average every n_downsample points to one value to smooth and downsample"""
			n_downsample = 100
			if (extra_points := x.shape[0] % n_downsample) != 0:
				x = x[:-extra_points]
			return x.reshape(-1, n_downsample).mean(1).reshape(-1, 1)

		for variant in self.ap.get_variants():
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt <= MAX_CELL_LENGTH]
			growth_rates[variant] = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'Mass', 'instantaneous_growth_rate',
				remove_first=True, fun=downsample).squeeze() * 3600.

			if len(all_cells) >= MIN_LATE_CELL_INDEX:
				all_cells_gens = [int(c.split("/")[-2][-6:]) for c in all_cells]
				early_cell_index = [i for i, v in enumerate(all_cells_gens) if v < MIN_LATE_CELL_INDEX]
				late_cell_index = [i for i, v in enumerate(all_cells_gens) if v >= MIN_LATE_CELL_INDEX]

				early_cells = all_cells[early_cell_index]
				late_cells = all_cells[late_cell_index]

				dt_early_cells = dt[early_cell_index]
				dt_late_cells = dt[late_cell_index]

				doubling_times_early_gens[variant] = dt_early_cells[dt_early_cells < MAX_CELL_LENGTH]
				doubling_times_late_gens[variant] = dt_late_cells[dt_late_cells < MAX_CELL_LENGTH]

				growth_rates_early_gens[variant] = read_stacked_columns(early_cells[dt_early_cells < MAX_CELL_LENGTH], 'Mass',
															 'instantaneous_growth_rate',
															 remove_first=True, fun=downsample).squeeze() * 3600.

				growth_rates_late_gens[variant] = read_stacked_columns(late_cells[dt_late_cells < MAX_CELL_LENGTH], 'Mass',
																		'instantaneous_growth_rate',
																		remove_first=True,
																		fun=downsample).squeeze() * 3600.
		# LATE GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], doubling_times_late_gens, 'Doubling Time (min)')
		self.hist(axes[1], growth_rates_late_gens, 'Growth rates (1/hr)', bin_width=0.05, sf=2)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_lategens', metadata)

		axes[0].set_xlim([30, 185])
		axes[1].set_xlim([0, 2.5])
		exportFigure(plt, plotOutDir, plotOutFileName + '_lategens_trimmed', metadata)

		# EARLY GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], doubling_times_early_gens, 'Doubling Time (min)')
		self.hist(axes[1], growth_rates_early_gens, 'Growth rates (1/hr)', bin_width=0.05, sf=2)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_earlygens', metadata)

		axes[0].set_xlim([30, 185])
		axes[1].set_xlim([0, 2.5])
		exportFigure(plt, plotOutDir, plotOutFileName + '_earlygens_trimmed', metadata)

		# ALL GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], doubling_times, 'Doubling Time (min)')
		self.hist(axes[1], growth_rates, 'Growth rates (1/hr)', bin_width=0.05, sf=2)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_allgens', metadata)

		axes[0].set_xlim([30, 185])
		axes[1].set_xlim([0, 2.5])
		exportFigure(plt, plotOutDir, plotOutFileName + '_allgens_trimmed', metadata)


if __name__ == "__main__":
	Plot().cli()
