"""
Compares time-averaged and normalized tRNA cistron expression (total of charged
and uncharged) sets of sims through a scatterplot and labels interesting tRNAs.
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_bulk_molecules)
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader

BOUNDS=[0, 0.04]
RATIO_CUTOFF=1.5

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
			cell_paths = ap.get_cells()

			uncharged_tRNA_counts = np.mean(read_stacked_bulk_molecules(
				cell_paths, uncharged_tRNA_cistron_ids, remove_first=True
				)[0], axis=0)
			charged_tRNA_counts = np.mean(read_stacked_bulk_molecules(
				cell_paths, charged_tRNA_cistron_ids, remove_first=True
				)[0], axis=0)
			total_tRNA_counts = uncharged_tRNA_counts + charged_tRNA_counts

			return total_tRNA_counts / np.sum(total_tRNA_counts)

		total_tRNA_counts1 = read_sims(ap1)
		total_tRNA_counts2 = read_sims(ap2)

		## Get tRNAs of interest
		large_ratio = np.abs(np.log(total_tRNA_counts1 / total_tRNA_counts2)) > np.log(1.5)

		# These assume we are running a sim with rRNA operons
		tRNA_has_rRNA_mask = transcription.rna_data['is_rRNA'][transcription.rna_data['includes_tRNA']]
		tRNA_from_rRNA_mask = transcription.tRNA_cistron_tu_mapping_matrix.dot(tRNA_has_rRNA_mask).astype(bool)
		special_tRNAs = np.logical_or(large_ratio, tRNA_from_rRNA_mask)

		# Get tRNAs of the same amino acid as the previous special tRNAs
		aas_of_tRNAs = np.array([np.where(tRNA)[0][0] for tRNA in transcription.aa_from_trna.T])
		off_ratio_aas = list(set(aas_of_tRNAs[special_tRNAs]))
		off_ratio_associated_tRNAs = np.isin(aas_of_tRNAs, np.array(off_ratio_aas))
		special_tRNAs = np.logical_or(special_tRNAs, off_ratio_associated_tRNAs)

		## Plot counts of each cistron
		plt.figure(figsize=(16, 16))
		ax1 = plt.subplot()
		ax1.scatter(total_tRNA_counts1[~tRNA_from_rRNA_mask], total_tRNA_counts2[~tRNA_from_rRNA_mask], clip_on=False,
					c='b', label="Pure tRNA")
		if tRNA_from_rRNA_mask.any():
			ax1.scatter(total_tRNA_counts1[tRNA_from_rRNA_mask], total_tRNA_counts2[tRNA_from_rRNA_mask], clip_on=False,
					c='r', label="tRNA from rRNA operon")

		for idx in np.where(special_tRNAs)[0]:
			ax1.text(total_tRNA_counts1[idx], total_tRNA_counts2[idx], uncharged_tRNA_cistron_ids[idx],
					 fontsize=15)

		ax1.set_title('Normalized tRNA counts of two sets of sims', fontsize=25)
		ax1.plot(BOUNDS, BOUNDS)

		ax1.set_xlabel('Normalized tRNA counts, sim 1', fontsize=15, labelpad=15)
		ax1.set_ylabel('Normalized tRNA counts, sim 2', fontsize=15, labelpad=15)
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 15))
		ax1.spines["left"].set_position(("outward", 15))
		ax1.set_xlim(BOUNDS)
		ax1.set_ylim(BOUNDS)
		ax1.legend(loc=2, prop={'size': 15})



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
