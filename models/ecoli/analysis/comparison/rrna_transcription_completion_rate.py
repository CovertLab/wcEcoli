"""
Comparison plot to compare the completion rates of rRNA transcription events.
Plots the number of completed rRNA transcription events vs the number of rRNA
transcription events that were prematurely terminated by collisions with
replication forks.
"""

import os
from typing import Tuple

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


TU_ID_TO_RRNA_OPERON_ID = {
	'TU0-1181[c]': 'rrnA',
	'TU0-1182[c]': 'rrnB',
	'TU0-1183[c]': 'rrnC',
	'TU0-1191[c]': 'rrnD',
	'TU0-1186[c]': 'rrnE',
	'TU0-1187[c]': 'rrnG',
	'TU0-1189[c]': 'rrnH',
	}

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		def read_sims(ap, sim_data):
			# Get cell paths
			cell_paths = ap.get_cells(
				generation=np.arange(0, ap.n_generation),
				only_successful=True)
			n_cells = len(cell_paths)

			# Get rna ID attributes from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			transcript_elongation_reader = TableReader(
				os.path.join(sim_out_dir, 'TranscriptElongationListener'))
			te_rna_ids = transcript_elongation_reader.readAttribute('rnaIds')
			rnap_data_reader = TableReader(
				os.path.join(sim_out_dir, 'RnapData'))
			rd_rna_ids = rnap_data_reader.readAttribute('rnaIds')

			# Get indexes of rRNA operons in both listeners
			te_rrna_indexes = np.array([
				te_rna_ids.index(key) for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])
			rd_rna_indexes = np.array([
				rd_rna_ids.index(key) for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])

			# Get number of completed and failed rRNA transcript events
			completed_events = read_stacked_columns(
				cell_paths, 'TranscriptElongationListener',
				'countRnaSynthesized', ignore_exception=True,
				remove_first=True, fun=lambda x: x[:, te_rrna_indexes]
				)
			incomplete_events = read_stacked_columns(
				cell_paths, 'RnapData',
				'incomplete_transcription_event', ignore_exception=True,
				remove_first=True, fun=lambda x: x[:, te_rrna_indexes]
				)

			n_complete = completed_events.sum(axis=0) / n_cells
			n_incomplete = incomplete_events.sum(axis=0) / n_cells

			return n_complete, n_incomplete

		n1_complete, n1_incomplete = read_sims(ap1, sim_data1)
		n2_complete, n2_incomplete = read_sims(ap2, sim_data2)

		fig = plt.figure(figsize=(9, 3))
		gs = fig.add_gridspec(1, 2, width_ratios=(6.5, 1))

		# Plot bar plots for each rRNA operon
		ax0 = fig.add_subplot(gs[0, 0])
		for i, (c1, c2, ic1, ic2) in enumerate(
				zip(n1_complete, n2_complete, n1_incomplete, n2_incomplete)):
			ax0.bar(
				1.5 * i - 0.25, c1, width=0.5, alpha=0.5, color='#555555')
			ax0.bar(
				1.5 * i - 0.25, ic1, width=0.5, alpha=0.5, color='C3',
				bottom=c1)
			ax0.bar(
				1.5 * i + 0.25, c2, width=0.5, alpha=0.8, color='#555555')
			ax0.bar(
				1.5 * i + 0.25, ic2, width=0.5, alpha=0.8, color='C3',
				bottom=c2)

		ax0.set_xticks(1.5 * np.arange(7))
		ax0.set_xticklabels([v for v in TU_ID_TO_RRNA_OPERON_ID.values()])
		ax0.set_xlim([-0.8, 9.8])
		ax0.set_ylabel('Numbers of transcription events')
		ax0.spines['top'].set_visible(False)
		ax0.spines['right'].set_visible(False)

		# Plot bar plots for totals across all rRNA operons
		ax1 = fig.add_subplot(gs[0, 1])

		ax1.bar(
			-0.25, n1_complete.sum(), width=0.5, alpha=0.5, color='#555555',
			label='complete / reference')
		ax1.bar(
			-0.25, n1_incomplete.sum(), width=0.5, alpha=0.5, color='C3',
			bottom=n1_complete.sum(), label='incomplete / reference')
		ax1.bar(
			0.25, n2_complete.sum(), width=0.5, alpha=0.8, color='#555555',
			label='complete / input')
		ax1.bar(
			0.25, n2_incomplete.sum(), width=0.5, alpha=0.8, color='C3',
			bottom=n2_complete.sum(), label='incomplete / input')

		ax1.set_xticks([0])
		ax1.set_xticklabels(['all operons'])
		ax1.set_xlim([-0.8, 0.8])
		ax1.spines['top'].set_visible(False)
		ax1.spines['right'].set_visible(False)

		ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})

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
