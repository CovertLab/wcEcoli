"""
Template for comparison analysis plots
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

class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# TODO(Albert): we could make a histogram with a bar plot for the diff aas?
		#  in the operon version, also indicating how much comes from rRNA operons?
		ap1, _, _ = self.setup(reference_sim_dir)
		ap2, sim_data, _ = self.setup(input_sim_dir)

		transcription = sim_data.process.transcription
		uncharged_tRNA_cistron_ids = transcription.uncharged_trna_names
		charged_tRNA_cistron_ids = transcription.charged_trna_names

		def read_sims(ap):
			cell_paths = ap.get_cells(seed=[0])
			simOutDir = os.path.join(cell_paths[0], "simOut")
			bulk_reader = TableReader(os.path.join(simOutDir, 'BulkMolecules'))
			bulk_names = bulk_reader.readAttribute('objectNames')
			uncharged_tRNA_idxs = np.isin(bulk_names, uncharged_tRNA_cistron_ids)
			charged_tRNA_idxs = np.isin(bulk_names, charged_tRNA_cistron_ids)

			# Load data
			uncharged_tRNA_counts = np.mean(read_stacked_columns(
				cell_paths, 'BulkMolecules', 'counts')[:, uncharged_tRNA_idxs], axis=0)
			charged_tRNA_counts = np.mean(read_stacked_columns(
				cell_paths, 'BulkMolecules', 'counts')[:, charged_tRNA_idxs], axis=0)
			total_tRNA_counts = uncharged_tRNA_counts + charged_tRNA_counts
			return total_tRNA_counts

		total_tRNA_counts_off = read_sims(ap1)
		total_tRNA_counts_on = read_sims(ap2)

		## Get tRNAs of interest
		high_ratio = np.array(total_tRNA_counts_on/total_tRNA_counts_off > 1.5, dtype=bool)
		low_ratio = np.array(total_tRNA_counts_off / total_tRNA_counts_on > 1.5, dtype=bool)
		# These assume we are running a sim with rRNA operons
		tRNA_has_rRNA_mask = transcription.rna_data['is_rRNA'][transcription.rna_data['includes_tRNA']]
		tRNA_from_rRNA_mask = transcription.tRNA_cistron_tu_mapping_matrix.dot(tRNA_has_rRNA_mask).astype(bool)

		special_tRNAs = np.logical_or(high_ratio, low_ratio)
		special_tRNAs = np.logical_or(special_tRNAs, tRNA_from_rRNA_mask)

		# Get tRNAs of the same amino acid as the previous special tRNAs
		aas_of_tRNAs = np.array([np.where(tRNA)[0][0] for tRNA in transcription.aa_from_trna.T])
		off_ratio_aas = list(set(aas_of_tRNAs[special_tRNAs]))
		off_ratio_associated_tRNAs = np.isin(aas_of_tRNAs, np.array(off_ratio_aas))

		special_tRNAs = np.logical_or(special_tRNAs, off_ratio_associated_tRNAs)

		## Plot counts of each cistron
		plt.figure(figsize=(16, 16))
		ax1 = plt.subplot()
		ax1.scatter(total_tRNA_counts_off[~tRNA_from_rRNA_mask], total_tRNA_counts_on[~tRNA_from_rRNA_mask], clip_on=False,
					c='b')
		ax1.scatter(total_tRNA_counts_off[tRNA_from_rRNA_mask], total_tRNA_counts_on[tRNA_from_rRNA_mask], clip_on=False,
					c='r')

		for idx in np.where(special_tRNAs)[0]:
			ax1.text(total_tRNA_counts_off[idx], total_tRNA_counts_on[idx], uncharged_tRNA_cistron_ids[idx])

		ax1.set_title('tRNA Counts w/ and w/o rRNA operons')
		ax1.plot([0, 4000], [0, 4000])

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


if __name__ == '__main__':
	Plot().cli()
