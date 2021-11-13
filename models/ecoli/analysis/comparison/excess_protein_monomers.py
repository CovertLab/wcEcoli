"""
Template for comparison analysis plots
"""

from functools import reduce
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


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, inputDir1, plotOutDir, plotOutFileName, inputDir2, unused, metadata):
		# noinspection PyUnusedLocal
		ap1, sim_data1, _ = self.setup(inputDir1)
		# noinspection PyUnusedLocal
		ap2, sim_data2, _ = self.setup(inputDir2)

		# Load from sim_data
		all_subunit_ids = sim_data1.process.complexation.molecule_names
		all_complex_ids = sim_data1.process.complexation.ids_complexes
		stoich_matrix_monomers = sim_data1.process.complexation.stoich_matrix_monomers()
		monomer_id_to_cistron_id = {
			monomer['id']: monomer['cistron_id']
			for monomer in sim_data1.process.translation.monomer_data
			}
		cistron_id_to_rna_indexes = sim_data2.process.transcription.cistron_id_to_rna_indexes

		monomer_counts_reader = TableReader(
			os.path.join(ap1.get_cells()[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)}

		complex_id_to_subunit_ids = {}
		complex_id_to_subunit_stoichs = {}

		for i, complex_id in enumerate(all_complex_ids):
			subunit_mask = stoich_matrix_monomers[:, i] < 0
			subunit_indexes = np.where(subunit_mask)[0]
			complex_id_to_subunit_stoichs[complex_id] = -stoich_matrix_monomers[:, i][subunit_mask]
			complex_id_to_subunit_ids[complex_id] = [
				all_subunit_ids[i] for i in subunit_indexes]

		def is_cotranscribed(subunit_ids):
			"""
			Returns True if the cistrons that each encode for the list of
			protein monomer subunits are ever found on the same transcription
			unit.
			"""
			# Get IDs of cistrons that encode for each monomer
			cistron_ids = [
				monomer_id_to_cistron_id[monomer_id]
				for monomer_id in subunit_ids]

			# Get indexes of transcription units that contain each cistron
			rna_indexes = [
				cistron_id_to_rna_indexes(cistron_id)
				for cistron_id in cistron_ids]

			# Find any overlapping transcription unit indexes
			return len(reduce(np.intersect1d, rna_indexes)) > 0


		def read_sims(ap):
			all_monomer_counts = read_stacked_columns(
				ap.get_cells(), 'MonomerCounts', 'monomerCounts', remove_first=True)

			# Take averages across timepoints
			all_monomer_counts_mean = all_monomer_counts.mean(axis=0)

			complex_id_to_excess_monomer_index = {}
			complex_id_to_is_cotranscribed = {}

			# Loop through each protein complex
			for i, complex_id in enumerate(all_complex_ids):
				# Get stoichiometries and IDs of each monomer subunit
				subunit_stoichs = complex_id_to_subunit_stoichs[complex_id]
				subunit_ids = complex_id_to_subunit_ids[complex_id]

				# Skip homogeneous protein complexes
				if len(subunit_ids) == 1:
					continue

				# Get indexes of monomer subunits in the monomer counts table
				try:
					monomer_indexes = np.array([
						monomer_id_to_index[subunit_id] for subunit_id in subunit_ids
						])
				except KeyError:
					# Skip complexes with non-protein monomers or monomers that
					# were excluded from the model
					continue

				# Get counts of monomers and rescale with stoichiometries
				monomer_counts = all_monomer_counts_mean[monomer_indexes]
				rescaled_monomer_counts = monomer_counts / subunit_stoichs

				# Skip complexes with zero counts
				if np.all(rescaled_monomer_counts == 0):
					continue

				# Calculate excess monomer index and determine if monomers are
				# cotranscribed
				complex_id_to_excess_monomer_index[complex_id] = (
					(rescaled_monomer_counts.max() - rescaled_monomer_counts.min())
					/ rescaled_monomer_counts.sum())
				complex_id_to_is_cotranscribed[complex_id] = is_cotranscribed(subunit_ids)

			return complex_id_to_excess_monomer_index, complex_id_to_is_cotranscribed

		excess_monomer_index1, is_cotranscribed1 = read_sims(ap1)
		excess_monomer_index2, is_cotranscribed2 = read_sims(ap2)

		complexes_to_plot = [
			complex_id for complex_id in excess_monomer_index1.keys()
			if (complex_id in excess_monomer_index2) and (is_cotranscribed1[complex_id])
			]

		# Overlay swarm and violin plots
		plt.figure(figsize=(6, 6))

		for complex_id in complexes_to_plot:
			if excess_monomer_index2[complex_id] > excess_monomer_index1[complex_id]:
				c = 'k'
				print(f'{complex_id}\t{excess_monomer_index1[complex_id]}\t{excess_monomer_index2[complex_id]}')
			else:
				c = 'grey'
			plt.plot(
				[0, 1],
				[excess_monomer_index1[complex_id], excess_monomer_index2[complex_id]],
				color = c)

		plt.xticks([])
		plt.ylabel('Excess monomer index')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		import ipdb; ipdb.set_trace()


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
