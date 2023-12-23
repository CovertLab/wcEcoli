"""
Compares timetraces of instantaneous doubling times between two sets of
simulations.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader, TableReaderError
from wholecell.utils import units


WINDOW_SIZE = 60
WINDOW_INTERVALS = 2

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		if metadata['variant'] != 'timelines' or ap1.get_variants()[0] != 28 or ap2.get_variants()[0] != 28:
			print('Skipping analysis -- this analysis only runs on timeline'
				  'variants with amino acid upshifts')

		def read_sims(ap):
			all_doubling_times = []
			t_max = 0

			# Get paths into each seed
			for seed in np.arange(ap.n_seed):
				cell_paths = ap.get_cells(seed=[seed], only_successful=True)

				# Only get data for sims that completed all generations
				if len(cell_paths) < ap.n_generation:
					continue

				time = read_stacked_columns(cell_paths, 'Main', 'time')

				# Get instantaneous doubling times
				growth_rates = (1 / units.s) * read_stacked_columns(
					cell_paths, 'Mass', 'instantaneous_growth_rate')
				doubling_times = (
					(1 / growth_rates) * np.log(2)).asNumber(units.min)

				all_doubling_times.append((time[1:], doubling_times[1:]))
				if time[-1] > t_max:
					t_max = time[-1]

			# Calculate sliding window averages of doubling times across
			# multiple seeds
			n_timetraces = len(all_doubling_times)
			swa_time = []
			swa_doubling_times = []

			for t_window_min in np.arange(0, t_max - WINDOW_SIZE, WINDOW_INTERVALS):
				t_window_max = t_window_min + WINDOW_SIZE
				dt_this_window = []

				for (t, dt) in all_doubling_times:
					mask = np.logical_and(t_window_min <= t, t < t_window_max)

					# Skip seed if no data exists for this window
					if not np.any(mask):
						continue
					dt_this_window.append(dt[mask].mean())

				# Skip window if no data exists for more than half of the
				# successful seeds
				if len(dt_this_window) < n_timetraces / 2:
					continue

				swa_time.append((t_window_min + t_window_max) / 2)
				swa_doubling_times.append(np.mean(dt_this_window))

			swa_doubling_times = (np.array(swa_time), np.array(swa_doubling_times))

			return all_doubling_times, swa_doubling_times, t_max

		all_dt1, swa_dt1, t_max1 = read_sims(ap1)
		all_dt2, swa_dt2, t_max2 = read_sims(ap2)

		plt.figure(figsize=(10, 4))

		# Plot timetraces of doubling times from two sims on same plot
		ax1 = plt.subplot(1, 1, 1)
		for (t, conc) in all_dt1:
			ax1.plot(t / 60, conc, c='C0', lw=0.5, alpha=0.1)
		for (t, conc) in all_dt2:
			ax1.plot(t / 60, conc, c='C1', lw=0.5, alpha=0.1)

		# Plot sliding window time averages with thicker lines
		ax1.plot(
			swa_dt1[0] / 60, swa_dt1[1], c='C0', lw=3, label='reference')
		ax1.plot(
			swa_dt2[0] / 60, swa_dt2[1], c='C1', lw=3, label='input')
		ax1.set_xlabel('Time (min)')
		ax1.set_ylabel('Doubling times (min)')
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 10))
		ax1.spines["left"].set_position(("outward", 10))
		ax1.set_xlim([0, max(t_max1, t_max2) / 60])
		ax1.set_ylim([0, 100])
		ax1.legend(loc=1)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

	def setup(self, inputDir: str) -> Tuple[
		AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data

if __name__ == "__main__":
	Plot().cli()
