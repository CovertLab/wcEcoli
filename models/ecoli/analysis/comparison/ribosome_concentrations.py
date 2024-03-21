"""
Compares timetraces of ribosome concentrations between two sets of simulations.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
import os

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns, read_stacked_bulk_molecules)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


WINDOW_SIZE = 60
WINDOW_INTERVALS = 2

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, validation_data2 = self.setup(input_sim_dir)

		s30_full_complex_id = sim_data1.molecule_ids.s30_full_complex
		s50_full_complex_id = sim_data1.molecule_ids.s50_full_complex
		n_avogadro = sim_data1.constants.n_avogadro

		def read_sims(ap):
			all_ribosome_timetraces = []
			t_max = 0

			# Get timetraces of ribosome concentrations from each seed
			for seed in np.arange(ap.n_seed):
				cell_paths = ap.get_cells(seed=[seed], only_successful=True)

				# Only get data for sims that completed all generations
				if len(cell_paths) < ap.n_generation:
					continue

				time = read_stacked_columns(
					cell_paths, 'Main', 'time', remove_first=True).squeeze()

				# Get index of active ribosomes in the unique molecule counts
				# reader
				unique_molecule_counts_reader = TableReader(
					os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
				unique_molecule_ids = unique_molecule_counts_reader.readAttribute(
					'uniqueMoleculeIds')
				active_ribosome_idx = unique_molecule_ids.index('active_ribosome')

				# Get counts of active ribosomes
				active_ribosome_counts = read_stacked_columns(
					cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
					remove_first=True,
					ignore_exception=True)[:, active_ribosome_idx]

				# Get counts of inactive ribosomal subunits
				ribosome_subunit_counts = read_stacked_bulk_molecules(
					cell_paths, ([s30_full_complex_id, s50_full_complex_id],),
					remove_first=True, ignore_exception=True)[0]

				# Get cell volume
				cell_volume = (units.L) * 1e-15 * read_stacked_columns(
					cell_paths, 'Mass', 'cellVolume', remove_first=True,
					ignore_exception=True
					).flatten()

				# Calculate concentrations of all ribosomes in uM
				ribosome_conc = ((1 / (n_avogadro * cell_volume)) * (
					active_ribosome_counts + ribosome_subunit_counts.min(axis=1)
					)).asNumber(units.mol/units.L) * 1e6

				all_ribosome_timetraces.append((time[1:], ribosome_conc[1:]))
				if time[-1] > t_max:
					t_max = time[-1]

			# Calculate sliding window averages of concentrations across
			# multiple seeds
			n_timetraces = len(all_ribosome_timetraces)
			swa_time = []
			swa_conc = []

			for t_window_min in np.arange(0, t_max - WINDOW_SIZE, WINDOW_INTERVALS):
				t_window_max = t_window_min + WINDOW_SIZE
				conc_this_window = []

				for (t, conc) in all_ribosome_timetraces:
					mask = np.logical_and(t_window_min <= t, t < t_window_max)

					# Skip seed if no data exists for this window
					if not np.any(mask):
						continue
					conc_this_window.append(conc[mask].mean())

				# Skip window if no data exists for more than half of the
				# successful seeds
				if len(conc_this_window) < n_timetraces/2:
					continue

				swa_time.append((t_window_min + t_window_max) / 2)
				swa_conc.append(np.mean(conc_this_window))

			swa_ribosome_timetrace = (np.array(swa_time), np.array(swa_conc))

			return all_ribosome_timetraces, swa_ribosome_timetrace, t_max

		all_tt1, swa_tt1, t_max1 = read_sims(ap1)
		all_tt2, swa_tt2, t_max2 = read_sims(ap2)

		plt.figure(figsize=(10, 3))

		# Plot timetraces of ribosome concentrations from two sims on same plot
		ax1 = plt.subplot(1, 1, 1)
		for (t, conc) in all_tt1:
			ax1.plot(t / 60, conc, clip_on=False, c='C0', lw=0.5, alpha=0.1)
		for (t, conc) in all_tt2:
			ax1.plot(t / 60, conc, clip_on=False, c='C1', lw=0.5, alpha=0.1)

		# Plot sliding window time averages with thicker lines
		ax1.plot(
			swa_tt1[0] / 60, swa_tt1[1], clip_on=False, c='C0', lw=3,
			label='reference')
		ax1.plot(
			swa_tt2[0] / 60, swa_tt2[1], clip_on=False, c='C1', lw=3,
			label='input')
		ax1.set_xlabel('Time (min)')
		ax1.set_ylabel('Ribosome concentration ($\mu$M)')
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 10))
		ax1.spines["left"].set_position(("outward", 10))
		ax1.set_xlim([0, max(t_max1, t_max2) / 60])
		ax1.set_ylim([0, 60])
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
