"""
Compares timetraces of instantaneous doubling times across different variants.
"""

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
from wholecell.analysis.plotting_tools import (
	DEFAULT_MATPLOTLIB_COLORS as COLORS)
from wholecell.utils import units


WINDOW_SIZE = 60
WINDOW_INTERVALS = 2

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		n_total_seeds = self.ap.n_seed
		n_total_gens = self.ap.n_generation
		variants = self.ap.get_variants()
		variant_to_all_timetraces = {}
		variant_to_swa_timetrace = {}
		variant_to_t_max = {}

		for variant in variants:
			all_doubling_times = []
			t_max = 0

			# Get timetraces of doubling times from each seed in this variant
			for seed in np.arange(n_total_seeds):
				cell_paths = self.ap.get_cells(
					variant=[variant], seed=[seed], only_successful=True)

				# Only get data for sims that completed all generations
				if len(cell_paths) < n_total_gens:
					continue

				time = read_stacked_columns(cell_paths, 'Main', 'time')
				growth_rates = (1 / units.s) * read_stacked_columns(
					cell_paths, 'Mass', 'instantaneous_growth_rate')
				doubling_times = (
					(1 / growth_rates) * np.log(2)).asNumber(units.min)

				all_doubling_times.append((time[1:], doubling_times[1:]))
				if time[-1] > t_max:
					t_max = time[-1]

			# Calculate sliding window averages of concentrations across
			# multiple seeds
			n_timetraces = len(all_doubling_times)
			swa_time = []
			swa_conc = []

			for t_window_min in np.arange(0, t_max - WINDOW_SIZE, WINDOW_INTERVALS):
				t_window_max = t_window_min + WINDOW_SIZE
				conc_this_window = []

				for (t, conc) in all_doubling_times:
					mask = np.logical_and(t_window_min <= t, t < t_window_max)

					# Skip seed if no data exists for this window
					if not np.any(mask):
						continue
					conc_this_window.append(conc[mask].mean())

				# Skip window if no data exists for more than half of the
				# successful seeds
				if len(conc_this_window) < n_timetraces / 2:
					continue

				swa_time.append((t_window_min + t_window_max) / 2)
				swa_conc.append(np.mean(conc_this_window))

			swa_doubling_time = (np.array(swa_time), np.array(swa_conc))

			# Store values for this variant in dictionary
			variant_to_all_timetraces[variant] = all_doubling_times
			variant_to_swa_timetrace[variant] = swa_doubling_time
			variant_to_t_max[variant] = t_max

		plt.figure(figsize=(10, 4))
		ax1 = plt.subplot(1, 1, 1)

		# Plot timetraces of doubling times from each variant on one plot
		for i, all_timetraces in enumerate(variant_to_all_timetraces.values()):
			color = COLORS[i % len(COLORS)]
			for (t, conc) in all_timetraces:
				ax1.plot(t / 60, conc, clip_on=False, c=color, lw=0.5, alpha=0.1)

		# Plot sliding window time averages with thicker lines
		for i, (variant, swa_timetrace) in enumerate(
				variant_to_swa_timetrace.items()):
			color = COLORS[i % len(COLORS)]
			ax1.plot(
				swa_timetrace[0] / 60, swa_timetrace[1], clip_on=False, c=color,
				lw=3, label=f'variant {variant}')

		ax1.set_xlabel('Time (min)')
		ax1.set_ylabel('Doubling times (min)')
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 10))
		ax1.spines["left"].set_position(("outward", 10))
		ax1.set_xlim([0, max([t for t in variant_to_t_max.values()]) / 60])
		ax1.set_ylim([0, 200])
		ax1.legend(loc=1)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
