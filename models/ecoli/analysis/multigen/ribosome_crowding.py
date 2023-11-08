"""
Comparison of target translation probabilities vs actual translation
probabilities for mRNAs whose translation probabilities exceeded the limit
set by the physical size and the elongation rates of ribosomes.
"""

import os
import pickle

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Listeners used
		ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))

		# Load data
		time_min = read_stacked_columns(cell_paths, 'Main', 'time') / 60
		monomer_ids = ribosome_reader.readAttribute('monomerIds')
		target_prob_translation_per_transcript = read_stacked_columns(
			cell_paths, 'RibosomeData', 'target_prob_translation_per_transcript')
		actual_prob_translation_per_transcript = read_stacked_columns(
			cell_paths, 'RibosomeData', 'actual_prob_translation_per_transcript')

		# Get indexes of proteins corresponding to mRNAs that were overcrowded
		# at some point in the sim
		overcrowded_monomer_indexes = np.where(
			(target_prob_translation_per_transcript
			 - actual_prob_translation_per_transcript).max(axis=0) > 0)[0]
		n_overcrowded_monomers = len(overcrowded_monomer_indexes)

		# Determine the gene ids corresponding to these proteins
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		monomer_to_mRNA_id_dict = dict(zip(monomer_sim_data['id'],
										   monomer_sim_data['cistron_id']))
		mRNA_to_gene_id_dict = dict(zip(mRNA_sim_data['id'],
										   mRNA_sim_data['gene_id']))
		overcrowded_gene_ids = [mRNA_to_gene_id_dict.get(monomer_to_mRNA_id_dict.get(
			monomer_ids[monomer_index])) for monomer_index in
			overcrowded_monomer_indexes]

		# Plot the target vs actual rna synthesis probabilites of these mRNAs
		plt.figure(figsize=(6, 1.5*n_overcrowded_monomers))

		for i, monomer_index in enumerate(overcrowded_monomer_indexes):
			target_prob_this_monomer = target_prob_translation_per_transcript[
				:, monomer_index]
			actual_prob_this_monomer = actual_prob_translation_per_transcript[
				:, monomer_index]

			ax = plt.subplot(n_overcrowded_monomers, 1, i + 1)
			ax.plot(time_min, target_prob_this_monomer, label='target')
			ax.plot(time_min, actual_prob_this_monomer, label='actual')
			ax.set_ylabel(f'{overcrowded_gene_ids[i]}\ntranslation probs')

			if i == 0:
				ax.set_title(f'Total number of proteins '
							 f'corresponding to overcrowded mRNAs: '
							 f'{n_overcrowded_monomers}')
				ax.legend(loc=1)

			if i == n_overcrowded_monomers - 1:
				ax.set_xlabel('Time [min]')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == '__main__':
	Plot().cli()
