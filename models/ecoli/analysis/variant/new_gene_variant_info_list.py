import pickle
import os

import csv
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = 16


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


		# Get list of cistron IDs from sim_data
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_ids = cistron_data['id']
		# Filter list for cistron IDs with associated protein ids
		cistron_id_to_protein_id = {
			protein['cistron_id']: protein['id']
			for protein in sim_data.process.translation.monomer_data}
		mRNA_cistron_ids = [
			cistron_id for cistron_id in cistron_ids if cistron_id in cistron_id_to_protein_id]
		# Get IDs of associated monomers and genes
		monomer_ids = [
			cistron_id_to_protein_id.get(cistron_id, None)
			for cistron_id in mRNA_cistron_ids]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data}
		gene_ids = [
			cistron_id_to_gene_id[cistron_id] for cistron_id in mRNA_cistron_ids]

		# Load tables and attributes for proteins
		monomer_reader = TableReader(
			os.path.join(all_cells[0], 'simOut', 'MonomerCounts'))
		monomer_ids_table = np.array(monomer_reader.readAttribute('monomerIds'))
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids_table)}
		monomer_indexes = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids_table])

		doubling_times = np.zeros(len(selected_variant_indexes))


		essential_genes = validation_data.essential_genes.essential_genes
		is_essential = np.array(
			[True if gene_id in essential_genes else False for gene_id in gene_ids])
		total_genes = len(gene_ids)
		total_essential_genes = len(essential_genes)
		is_essential = np.reshape(is_essential, (len(is_essential), 1))
		gene_ids = np.reshape(gene_ids, (len(gene_ids), 1))
		total_timesteps = np.zeros((len(monomer_ids), 1))
		total_individual_cells = np.zeros((len(monomer_ids), 1))


		monomer_ids = np.reshape(monomer_ids, (len(monomer_ids), 1))
		# Initialize a pandas dataframe to store the data
		all_avg_monomer_counts = np.zeros((len(monomer_ids), 1))

		variant_table = np.array(['']+['']+['']+['']+['']+['']+['']+[''])
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
			monomer_counts = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts',
				remove_first=True, ignore_exception=True)

			monomer_counts_avg = monomer_counts.mean(axis=0)
			monomer_counts_avg = monomer_counts_avg.reshape(-1, 1)


			all_avg_monomer_counts = np.copy(monomer_counts_avg)

			# count whether each gene's monomer exists in each
			# timestep or not
			monomer_zeros = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts',
			     ignore_exception=True, fun=lambda x: (x == 0).sum(axis=0))[:, monomer_indexes]
			total_timesteps[0:len(monomer_ids)-1,0] = monomer_counts.shape[0]

			# count number of cells where monomer goes to zero
			cell_has_zero_monomer = read_stacked_columns(
				all_cells, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True, fun=lambda x: (x == 0).any(axis=0))[:, monomer_indexes]
			cells_with_zero = cell_has_zero_monomer.sum(axis=0)

			total_individual_cells[0:len(monomer_ids)-1,0] = all_cells.shape[0]

			dt = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[i] = np.mean(dt)

		# and col names as selected_variant_indexes

			gene_ids = np.array(gene_ids).reshape(-1, 1)
			monomer_ids = np.array(monomer_ids).reshape(-1, 1)
			is_essential = np.array(is_essential).reshape(-1, 1)
			total_timesteps = np.full((len(monomer_ids), 1), monomer_counts.shape[0])
			cells_with_zero = cells_with_zero.reshape(-1, 1)
			monomer_zeros= monomer_zeros.reshape(-1, 1)
			total_individual_cells = np.full((len(monomer_ids), 1), len(all_cells))


			variant_table_i = np.hstack((gene_ids, monomer_ids, is_essential, monomer_zeros, total_timesteps, cells_with_zero, total_individual_cells, all_avg_monomer_counts ))
			header1 = np.array([f'variant number {i}']  + ['doubling time:'] + [doubling_times[i]] + ['']+['']+['']+['']+[''])
			header2 = np.array(['gene ids'] + ['monomer_ids'] + ['is_essential'] + ['number of timesteps monomer not present'] + ['total timesteps'] + ['number of cells where monomer disappears'] + ['total individual cells'] + ['avg monomer count'])

			variant_table = np.vstack((variant_table, header1, header2, variant_table_i))



		np.savetxt(
    		os.path.join(plotOutDir, 'variant_table.csv'), variant_table, delimiter=',', fmt='%s')
if __name__ == "__main__":
	Plot().cli()
