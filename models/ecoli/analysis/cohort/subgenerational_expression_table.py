"""
Generates a table of genes that are subgenerationally expressed, with their
expression frequencies and average/maximum mRNA/protein counts.
"""

import pickle
import os
import csv
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.io.tablereader import TableReader


IGNORE_FIRST_N_GENS = 8


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return

		# Get list of cistron IDs from sim_data
		cistron_data = sim_data.process.transcription.cistron_data
		cistron_ids = cistron_data['id']

		# Filter list for cistron IDs with associated protein ids
		cistron_id_to_protein_id = {
			protein['cistron_id']: protein['id']
			for protein in sim_data.process.translation.monomer_data
			}
		mRNA_cistron_ids = [
			cistron_id for cistron_id in cistron_ids
			if cistron_id in cistron_id_to_protein_id]

		# Get IDs of associated monomers and genes
		monomer_ids = [
			cistron_id_to_protein_id.get(cistron_id, None)
			for cistron_id in mRNA_cistron_ids]
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id'] for cistron in cistron_data
			}
		gene_ids = [
			cistron_id_to_gene_id[cistron_id]
			for cistron_id in mRNA_cistron_ids]

		# Get subcolumn for mRNA cistron IDs in RNA counts table
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		
		print(f"Processing {len(cell_paths)} cells...")
		
		simOutDir = os.path.join(cell_paths[0], 'simOut')
		rna_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_ids_rna_counts_table = rna_counts_reader.readAttribute(
			'mRNA_cistron_ids')

		# Get indexes of mRNA cistrons in this subcolumn
		mRNA_cistron_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(mRNA_cistron_ids_rna_counts_table)
			}
		mRNA_cistron_indexes = np.array([
			mRNA_cistron_id_to_index[cistron_id] for cistron_id
			in mRNA_cistron_ids
			])

		# Initialize arrays for accumulating statistics
		n_genes = len(mRNA_cistron_ids)
		n_cells = len(cell_paths)
		mRNA_exists_count = np.zeros(n_genes, dtype=int)
		max_mRNA_counts = np.zeros(n_genes, dtype=int)
		
		# Process RNA counts cell by cell to save memory
		print("Processing mRNA counts...")
		for i, cell_path in enumerate(cell_paths):
			if i % 100 == 0:
				print(f"  Cell {i}/{n_cells}")
			try:
				simOutDir = os.path.join(cell_path, 'simOut')
				rna_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
				mRNA_counts = rna_reader.readColumn('mRNA_cistron_counts')[:, mRNA_cistron_indexes]
				
				# Update existence count
				mRNA_exists_count += (mRNA_counts.sum(axis=0) > 0).astype(int)
				
				# Update max counts
				cell_max = mRNA_counts.max(axis=0)
				max_mRNA_counts = np.maximum(max_mRNA_counts, cell_max)
				
			except Exception as e:
				print(f"  Warning: Could not read cell {cell_path}: {e}")
				continue

		# Calculate probability
		p_mRNA_exists_in_gen = mRNA_exists_count / n_cells

		# Get subcolumn for monomer IDs in monomer counts table
		simOutDir = os.path.join(cell_paths[0], 'simOut')
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, 'MonomerCounts'))
		monomer_ids_monomer_counts_table = monomer_counts_reader.readAttribute(
			'monomerIds')

		# Get indexes of monomers in this subcolumn
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids_monomer_counts_table)
			}
		monomer_indexes = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids
			])

		# Process monomer counts cell by cell to save memory
		max_monomer_counts = np.zeros(n_genes, dtype=int)
		print("Processing monomer counts...")
		for i, cell_path in enumerate(cell_paths):
			if i % 100 == 0:
				print(f"  Cell {i}/{n_cells}")
			try:
				simOutDir = os.path.join(cell_path, 'simOut')
				monomer_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))
				monomer_counts = monomer_reader.readColumn('monomerCounts')[:, monomer_indexes]
				
				# Update max counts
				cell_max = monomer_counts.max(axis=0)
				max_monomer_counts = np.maximum(max_monomer_counts, cell_max)
				
			except Exception as e:
				print(f"  Warning: Could not read cell {cell_path}: {e}")
				continue

		# Get indexes of subgenerationally expressed genes
		subgen_exp_indexes = np.where(np.logical_and(
			p_mRNA_exists_in_gen > 0, p_mRNA_exists_in_gen < 1
			))[0]

		print(f"Found {len(subgen_exp_indexes)} subgenerationally expressed genes")

		# Write data to table
		output_file = os.path.join(plotOutDir, plotOutFileName + '.tsv')
		print(f"Writing results to {output_file}")
		with open(output_file, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'gene_name', 'cistron_name', 'protein_name',
				'p_expressed', 'max_mRNA_count', 'max_protein_count'
				])

			for i in subgen_exp_indexes:
				writer.writerow([
					gene_ids[i], mRNA_cistron_ids[i], monomer_ids[i][:-3],
					p_mRNA_exists_in_gen[i], max_mRNA_counts[i],
					max_monomer_counts[i]
					])
		
		print("Done!")


if __name__ == '__main__':
	Plot().cli()