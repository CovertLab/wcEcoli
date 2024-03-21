"""
Comparison plot to compare the averaged densities of RNAPs on genes encoding
rRNAs.
"""

import os
from typing import Tuple

from matplotlib import pyplot as plt
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
from wholecell.utils import units


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

		#if ap1.n_generation <= 4 or ap2.n_generation <= 4:
		#	print('Not enough generations to run analysis.')
		#	return

		rnap_footprint_size = sim_data1.process.transcription.active_rnap_footprint_size.asNumber(
			units.nt)

		def read_sims(ap, sim_data):
			# Get RNAP elongation rate used in simulations
			nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
			rnap_elong_rate = sim_data.process.transcription.rnaPolymeraseElongationRateDict[
				nutrients].asNumber(units.nt/units.s)

			# Get cell paths
			cell_paths = ap.get_cells(
				generation=np.arange(0, ap.n_generation),
				only_successful=True)

			# Get rna ID attributes from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rnap_data_reader = TableReader(
				os.path.join(sim_out_dir, 'RnapData'))
			rnap_data_rna_ids = rnap_data_reader.readAttribute('rnaIds')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			rna_synth_prob_rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')

			# Get indexes of rRNA operons in both listeners
			rnap_data_rrna_indexes = np.array([
				rnap_data_rna_ids.index(key)
				for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])
			rna_synth_prob_rrna_indexes = np.array([
				rna_synth_prob_rna_ids.index(key)
				for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])

			# Get number of rRNA initiation events and copy numbers
			initiation_events = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True, remove_first=True,
				fun=lambda x: x[:, rnap_data_rrna_indexes]
				)
			copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'promoter_copy_number',
				ignore_exception=True, remove_first=True,
				fun=lambda x: x[:, rna_synth_prob_rrna_indexes])

			# Calculate total number of initiation events per copy
			initiation_events_per_copy_each_rrna = (
				initiation_events / copy_numbers
				).sum(axis=0)

			# Get mask for rRNAs with zero initiation events throughout all
			# sims - these rRNAs are assumed to be nonfunctional (most likely
			# knocked out)
			nonfunctional_rrna_mask = (
				initiation_events_per_copy_each_rrna == 0)

			# Calculate initiation events per copy amongst all functional rRNAs
			initiation_events_per_copy_all_functional_rrnas = (
				initiation_events[:, ~nonfunctional_rrna_mask].sum(axis=1)
				/ copy_numbers[:, ~nonfunctional_rrna_mask].sum(axis=1)
				).sum()

			# Get total simulation time
			generation_lengths = read_stacked_columns(
				cell_paths, 'Main', 'time',
				ignore_exception=True,
				fun=lambda x: x[-1] - x[0])
			total_simulation_time = generation_lengths.sum()

			# Calculate average distance between adjacent RNAPs for each operon
			# and for all operons
			with np.errstate(divide='ignore'):
				avg_rrna_distance_each_rrna = (
					rnap_elong_rate * total_simulation_time
					/ initiation_events_per_copy_each_rrna)
			avg_rrna_distance_all_rrnas = (
				rnap_elong_rate * total_simulation_time
				/ initiation_events_per_copy_all_functional_rrnas)

			return avg_rrna_distance_each_rrna, avg_rrna_distance_all_rrnas

		d1, dt1 = read_sims(ap1, sim_data1)
		d2, dt2 = read_sims(ap2, sim_data2)

		fig = plt.figure(figsize=(9, 3))
		gs = fig.add_gridspec(1, 2, width_ratios=(6.5, 1))

		# Plot bar plots for each rRNA operon
		ax0 = fig.add_subplot(gs[0, 0])
		for i, (c1, c2) in enumerate(zip(d1, d2)):
			if c1 != np.inf:
				ax0.bar(
					1.5*i - 0.25, c1, width=0.5, alpha=0.5, color='C0',
					label='reference')
			# If value is equal to infinity (nonfunctional rRNA with zero
			# expression), mark with a red "x" at y=0
			else:
				ax0.scatter(
					1.5*i - 0.25, 0, marker='x', color='red', clip_on=False)
			if c2 != np.inf:
				ax0.bar(
					1.5*i + 0.25, c2, width=0.5, alpha=0.5, color='C1',
					label='input')
			else:
				ax0.scatter(
					1.5*i + 0.25, 0, marker='x', color='red', clip_on=False)

		ax0.axhline(y=rnap_footprint_size, ls='--', color='red')
		ax0.set_xticks(1.5*np.arange(7))
		ax0.set_xticklabels([v for v in TU_ID_TO_RRNA_OPERON_ID.values()])
		ax0.set_xlim([-0.8, 9.8])
		ax0.set_ylabel('Average distance\nbetween adjacent RNAPs [nt]')
		ax0.spines['top'].set_visible(False)
		ax0.spines['right'].set_visible(False)

		# Plot bar plots for averaged densities across all rRNA operons
		ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)

		ax1.bar(
			-0.25, dt1, width=0.5, alpha=0.5, color='C0', label='reference')
		ax1.bar(
			0.25, dt2, width=0.5, alpha=0.5, color='C1', label='input')

		ax1.axhline(y=rnap_footprint_size, ls='--', color='red')
		ax1.set_xticks([0])
		ax1.set_xticklabels([
			'all operons'])
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
