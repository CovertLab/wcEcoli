"""
Comparison plot to compare the percentages of rRNAs that are successfully
complexed into ribosomes before they are degraded between two sets of
simulations.
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



class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		ap1, sim_data1, _ = self.setup(reference_sim_dir)
		ap2, sim_data2, _ = self.setup(input_sim_dir)

		def read_sims(ap, sim_data):
			# Get IDs of rRNA cistrons
			cistron_data = sim_data.process.transcription.cistron_data
			rRNA_ids = cistron_data['id'][cistron_data['is_rRNA']]
			rRNA_23s_ids = cistron_data['id'][cistron_data['is_23S_rRNA']]
			rRNA_16s_ids = cistron_data['id'][cistron_data['is_16S_rRNA']]
			rRNA_5s_ids = cistron_data['id'][cistron_data['is_5S_rRNA']]

			# Get indexes of rRNAs in the RNA degradation listener
			rna_degradation_reader = TableReader(
				os.path.join(ap.get_cells()[0], 'simOut', 'RnaDegradationListener'))
			cistron_ids = rna_degradation_reader.readAttribute('cistron_ids')
			cistron_id_to_index = {
				cistron_id: i for (i, cistron_id) in enumerate(cistron_ids)}

			rRNA_indexes = np.array([
				cistron_id_to_index[rRNA_id] for rRNA_id in rRNA_ids
				])
			rRNA_id_to_rRNA_index = {
				cistron_id: i for (i, cistron_id) in enumerate(rRNA_ids)}

			rRNA_23s_indexes = np.array([
				rRNA_id_to_rRNA_index[rRNA_id] for rRNA_id in rRNA_23s_ids
				])
			rRNA_16s_indexes = np.array([
				rRNA_id_to_rRNA_index[rRNA_id] for rRNA_id in rRNA_16s_ids
				])
			rRNA_5s_indexes = np.array([
				rRNA_id_to_rRNA_index[rRNA_id] for rRNA_id in rRNA_5s_ids
				])

			# Get ids of complexation reactions for ribosomal subunits
			s50_subunit_id = sim_data.molecule_ids.s50_full_complex
			s30_subunit_id = sim_data.molecule_ids.s30_full_complex
			complex_id_to_rxn_id = {
				complex_id: rxn_id for (complex_id, rxn_id)
				in zip(sim_data.process.complexation.ids_complexes,
					   sim_data.process.complexation.ids_reactions)
				}
			s50_rxn_id = complex_id_to_rxn_id[s50_subunit_id]
			s30_rxn_id = complex_id_to_rxn_id[s30_subunit_id]

			# Get indexes of the reactions in the complexation listener
			complexation_reader = TableReader(
				os.path.join(ap.get_cells()[0], 'simOut', 'ComplexationListener'))
			complexation_rxn_ids = complexation_reader.readAttribute('reactionIDs')
			rxn_id_to_index = {
				rxn_id: i for (i, rxn_id) in enumerate(complexation_rxn_ids)}

			subunit_rxn_indexes = np.array([
				rxn_id_to_index[s50_rxn_id], rxn_id_to_index[s30_rxn_id]])

			# Get total counts of rRNA molecules degraded throughout all
			# successful sims
			cell_paths = ap.get_cells(only_successful=True)
			n_cells = len(cell_paths)
			counts_rna_degraded = read_stacked_columns(
				cell_paths, 'RnaDegradationListener',
				'count_RNA_degraded_per_cistron',
				ignore_exception=True, fun=lambda x: x[:, rRNA_indexes])
			total_23s_rRNA_degraded = counts_rna_degraded[:, rRNA_23s_indexes].sum()
			total_16s_rRNA_degraded = counts_rna_degraded[:, rRNA_16s_indexes].sum()
			total_5s_rRNA_degraded = counts_rna_degraded[:, rRNA_5s_indexes].sum()

			# Get total counts of subunit complexation events through all sims
			complexation_events = read_stacked_columns(
				cell_paths, 'ComplexationListener',
				'complexationEvents',
				ignore_exception=True, fun=lambda x: x[:, subunit_rxn_indexes])
			total_50s_complexed = complexation_events[:, 0].sum()
			total_30s_complexed = complexation_events[:, 1].sum()

			# Calculate the mean squared deviation between the two counts
			msd_excess_subunits_complexed = np.sqrt((
				(complexation_events[:, 0] - complexation_events[:, 1])**2).mean())

			# Divide by number of cells
			return {
				'avg_23s_rRNA_degraded': total_23s_rRNA_degraded / n_cells,
				'avg_16s_rRNA_degraded': total_16s_rRNA_degraded / n_cells,
				'avg_5s_rRNA_degraded': total_5s_rRNA_degraded / n_cells,
				'avg_50s_complexed': total_50s_complexed / n_cells,
				'avg_30s_complexed': total_30s_complexed / n_cells,
				'msd_excess_subunits_complexed': msd_excess_subunits_complexed,
			}

		data1 = read_sims(ap1, sim_data1)
		data2 = read_sims(ap2, sim_data2)

		fig = plt.figure(figsize=(12, 3))
		gs = fig.add_gridspec(1, 3, width_ratios=(3, 1, 1))

		# Plot bar plots for each rRNA
		ax0 = fig.add_subplot(gs[0, 0])
		ax0.bar(0, data1['avg_30s_complexed'], width=0.7, alpha=0.5,
			   color='#555555', label='complexed')
		ax0.bar(0, data1['avg_16s_rRNA_degraded'], width=0.7, alpha=0.5,
			   color='C3', bottom=data1['avg_30s_complexed'], label='degraded')
		ax0.bar(1, data2['avg_30s_complexed'], width=0.7, alpha=0.5,
			   color='#555555')
		ax0.bar(1, data2['avg_16s_rRNA_degraded'], width=0.7, alpha=0.5,
			   color='C3', bottom=data2['avg_30s_complexed'])

		ax0.bar(2.5, data1['avg_50s_complexed'], width=0.7, alpha=0.5,
			   color='#555555')
		ax0.bar(2.5, data1['avg_23s_rRNA_degraded'], width=0.7, alpha=0.5,
			   color='C3', bottom=data1['avg_50s_complexed'])
		ax0.bar(3.5, data2['avg_50s_complexed'], width=0.7, alpha=0.5,
			   color='#555555')
		ax0.bar(3.5, data2['avg_23s_rRNA_degraded'], width=0.7, alpha=0.5,
			   color='C3', bottom=data2['avg_50s_complexed'])

		ax0.bar(5, data1['avg_50s_complexed'], width=0.7, alpha=0.5,
			   color='#555555')
		ax0.bar(5, data1['avg_5s_rRNA_degraded'], width=0.7, alpha=0.5,
			   color='C3', bottom=data1['avg_50s_complexed'])
		ax0.bar(6, data2['avg_50s_complexed'], width=0.7, alpha=0.5,
			   color='#555555')
		ax0.bar(6, data2['avg_5s_rRNA_degraded'], width=0.7, alpha=0.5,
			   color='C3', bottom=data2['avg_50s_complexed'])

		ax0.set_xticks([0, 1, 2.5, 3.5, 5, 6])
		ax0.set_xticklabels([
			'reference,\n16S', 'input,\n16S',
			'reference,\n23S', 'input,\n23S',
			'reference,\n5S', 'input,\n5S'])
		ax0.set_xlim([-0.8, 6.8])
		ax0.set_ylabel('counts of rRNAs\ncomplexed/degraded')
		ax0.spines['top'].set_visible(False)
		ax0.spines['right'].set_visible(False)

		# Plot bar plots for total number of rRNAs
		ax1 = fig.add_subplot(gs[0, 1])
		avg_complexed_counts1 = data1['avg_30s_complexed'] + 2*data1['avg_50s_complexed']
		avg_complexed_counts2 = data2['avg_30s_complexed'] + 2*data2['avg_50s_complexed']
		avg_degraded_counts1 = (
			data1['avg_16s_rRNA_degraded']
			+ data1['avg_23s_rRNA_degraded']
			+ data1['avg_5s_rRNA_degraded'])
		avg_degraded_counts2 = (
			data2['avg_16s_rRNA_degraded']
			+ data2['avg_23s_rRNA_degraded']
			+ data2['avg_5s_rRNA_degraded'])

		ax1.bar(0, avg_complexed_counts1, width=0.7, alpha=0.5,
				color='#555555', label='complexed')
		ax1.bar(0, avg_degraded_counts1, width=0.7, alpha=0.5,
				color='C3', bottom=avg_complexed_counts1, label='degraded')
		ax1.bar(1, avg_complexed_counts2, width=0.7, alpha=0.5,
				color='#555555')
		ax1.bar(1, avg_degraded_counts2, width=0.7, alpha=0.5,
				color='C3', bottom=avg_complexed_counts2)

		ax1.set_xticks([0, 1])
		ax1.set_xticklabels([
			'reference,\nall rRNA', 'input,\nall rRNA'])
		ax1.set_xlim([-0.8, 1.8])
		ax1.spines['top'].set_visible(False)
		ax1.spines['right'].set_visible(False)

		ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})

		ax2 = fig.add_subplot(gs[0, 2])

		ax2.bar(0, data1['msd_excess_subunits_complexed'], width=0.7, alpha=0.5,
				color='#555555')
		ax2.bar(1, data2['msd_excess_subunits_complexed'], width=0.7, alpha=0.5,
				color='#555555')

		ax2.set_xticks([0, 1])
		ax2.set_xticklabels(['reference', 'input'])
		ax2.set_xlim([-0.8, 1.8])
		ax2.set_ylabel('MSD, excess subunits produced')
		ax2.spines['top'].set_visible(False)
		ax2.spines['right'].set_visible(False)

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
