"""
Comparison of protein counts from genes belonging to the same operon between two
sets of simulations.
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
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader

LOW_EXP_TU_ID = 'TU00199[c]'
HIGH_EXP_TU_ID = 'TU0-1123[c]'
GEN_RANGE = np.arange(0, 32)
GEN_TICKS = [0, 8, 16, 24]


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.

		# noinspection PyUnusedLocal
		ap1, _, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		cistron_id_to_monomer_id = {
			monomer['cistron_id']: monomer['id'] for monomer
			in sim_data2.process.translation.monomer_data
			}
		low_exp_operon_cistron_ids = [
			sim_data2.process.transcription.cistron_data['id'][i] for i
			in sim_data2.process.transcription.rna_id_to_cistron_indexes(LOW_EXP_TU_ID)]
		low_exp_operon_monomer_ids = [
			cistron_id_to_monomer_id[cistron_id] for cistron_id
			in low_exp_operon_cistron_ids]
		high_exp_operon_cistron_ids = [
			sim_data2.process.transcription.cistron_data['id'][i] for i
			in sim_data2.process.transcription.rna_id_to_cistron_indexes(HIGH_EXP_TU_ID)]
		high_exp_operon_monomer_ids = [
			cistron_id_to_monomer_id[cistron_id] for cistron_id
			in high_exp_operon_cistron_ids]

		def read_monomer_trace(ap):
			cell_paths = ap.get_cells(seed=[0], generation=GEN_RANGE)

			simOutDir = os.path.join(cell_paths[0], "simOut")
			monomer_counts_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))
			all_monomer_ids = monomer_counts_reader.readAttribute('monomerIds')

			# Get indexes of monomers translated from the TU
			low_exp_monomer_indexes = [
				all_monomer_ids.index(monomer_id) for monomer_id
				in low_exp_operon_monomer_ids]
			high_exp_monomer_indexes = [
				all_monomer_ids.index(monomer_id) for monomer_id
				in high_exp_operon_monomer_ids]

			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time')
			gen_start_time = (read_stacked_columns(
				cell_paths, 'Main', 'time', fun=lambda x: x[0]) / 60).flatten()
			all_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts')
			low_exp_monomer_counts = all_monomer_counts[:, low_exp_monomer_indexes]
			high_exp_monomer_counts = all_monomer_counts[:, high_exp_monomer_indexes]

			return time, gen_start_time, low_exp_monomer_counts, high_exp_monomer_counts

		t1, gen_start_time1, lc1, hc1 = read_monomer_trace(ap1)
		t2, gen_start_time2, lc2, hc2 = read_monomer_trace(ap2)

		plt.figure(figsize=(6, 6))
		gen_ticks1 = [gen_start_time1[x] for x in GEN_TICKS] + [t1.flatten()[-1] / 60]
		gen_ticks2 = [gen_start_time2[x] for x in GEN_TICKS] + [t2.flatten()[-1] / 60]
		gen_labels = GEN_TICKS + [32]

		ax0 = plt.subplot(2, 2, 1)
		ax0.plot(t1 / 60, lc1, clip_on=False)
		ax0.set_ylabel('Protein Counts')
		ax0.spines["top"].set_visible(False)
		ax0.spines["right"].set_visible(False)
		ax0.spines["bottom"].set_position(("outward", 10))
		ax0.spines["left"].set_position(("outward", 10))
		ax0.set_xlim([gen_ticks1[0], gen_ticks1[-1]])
		ax0.set_ylim([0, 250])
		ax0.set_yticks([0, 250])
		ax0.set_xticks(gen_ticks1)
		ax0.set_xticklabels(gen_labels)

		ax1 = plt.subplot(2, 2, 2)
		ax1.plot(t2 / 60, lc2, clip_on=False)
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 10))
		ax1.spines["left"].set_position(("outward", 10))
		ax1.set_xlim([gen_ticks2[0], gen_ticks2[-1]])
		ax1.set_ylim([0, 250])
		ax1.set_yticks([0, 250])
		ax1.set_xticks(gen_ticks2)
		ax1.set_xticklabels(gen_labels)

		ax2 = plt.subplot(2, 2, 3)
		ax2.plot(t1 / 60, hc1, clip_on=False)
		ax2.set_ylabel('Protein Counts')
		ax2.set_xlabel('Generations')
		ax2.spines["top"].set_visible(False)
		ax2.spines["right"].set_visible(False)
		ax2.spines["bottom"].set_position(("outward", 10))
		ax2.spines["left"].set_position(("outward", 10))
		ax2.set_xlim([gen_ticks1[0], gen_ticks1[-1]])
		ax2.set_ylim([0, 800])
		ax2.set_yticks([0, 800])
		ax2.set_xticks(gen_ticks1)
		ax2.set_xticklabels(gen_labels)

		ax3 = plt.subplot(2, 2, 4)
		ax3.plot(t2 / 60, hc2, clip_on=False)
		ax3.set_xlabel('Generations')
		ax3.spines["top"].set_visible(False)
		ax3.spines["right"].set_visible(False)
		ax3.spines["bottom"].set_position(("outward", 10))
		ax3.spines["left"].set_position(("outward", 10))
		ax3.set_xlim([gen_ticks2[0], gen_ticks2[-1]])
		ax3.set_ylim([0, 800])
		ax3.set_yticks([0, 800])
		ax3.set_xticks(gen_ticks2)
		ax3.set_xticklabels(gen_labels)

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
