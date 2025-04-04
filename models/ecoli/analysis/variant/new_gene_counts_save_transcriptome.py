"""
Save the average RNA counts for each variant index as a column in a csv file.
"""

import os
import pickle

import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16
# IGNORE_FIRST_N_GENS = 1

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		n_total_gens = self.ap.n_generation

		selected_variant_indexes = self.ap.get_variants()

		n_variants = len(selected_variant_indexes)

		# Initialize everything
		all_cells = self.ap.get_cells(
			variant=[selected_variant_indexes[0]],
			generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
			only_successful=True)

		# Load tables and attributes for RNA cistrons
		RNA_reader = TableReader(
			os.path.join(all_cells[0], 'simOut', 'RNACounts'))
		RNA_cistron_ids = np.array(RNA_reader.readAttribute('mRNA_cistron_ids'))

		# Determine which RNAs correspond to essential genes
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data}
		gene_ids = [
			cistron_id_to_gene_id[cistron_id] for cistron_id in RNA_cistron_ids]

		#TODO: rRNA have to be done separately?

		essential_genes = validation_data.essential_genes.essential_genes
		is_essential = np.array(
			[True if gene_id in essential_genes else False for gene_id in gene_ids])

		is_essential = np.reshape(is_essential, (len(is_essential), 1))
		RNA_cistron_ids = np.reshape(RNA_cistron_ids, (len(RNA_cistron_ids), 1))
		# Initialize a pandas dataframe to store the data
		all_avg_RNA_cistron_counts = np.zeros((len(RNA_cistron_ids), n_variants))

		# Loop through variant indexes
		for i, variant_index in enumerate(selected_variant_indexes):
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Read columns
			# remove_first=True because countsToMolar is 0 at first time step
			RNA_cistron_counts = read_stacked_columns(
				all_cells, 'RNACounts', 'mRNA_cistron_counts',
				remove_first=True, ignore_exception=True)

			RNA_cistron_counts_avg = RNA_cistron_counts.mean(axis=0)

			all_avg_RNA_cistron_counts[:, i] = np.copy(RNA_cistron_counts_avg)

		all_labeled_avg_RNA_cistron_counts = np.hstack((RNA_cistron_ids, is_essential, all_avg_RNA_cistron_counts))
		header = np.array(['RNA_cistron_ids'] + ['is_essential'] + selected_variant_indexes)
		all_labeled_avg_RNA_cistron_counts = np.vstack((header, all_labeled_avg_RNA_cistron_counts))

		np.savetxt(
			os.path.join(plotOutDir, 'all_labeled_avg_RNA_cistron_counts.csv'),
			all_labeled_avg_RNA_cistron_counts, delimiter=',', fmt='%s')

if __name__ == "__main__":
	Plot().cli()
