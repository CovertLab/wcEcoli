"""
Generates a table of genes that are subgenerationally expressed, with their
expression frequencies and average/maximum mRNA/protein counts for each variant
index.
"""

import pickle
import os

import csv
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = 16

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		variants = self.ap.get_variants()

		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return

		# Get ids and indexes of cistrons and monomers
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
		# Get subcolumn for mRNA cistron IDs in RNA counts table
		cell_paths_initialize = self.ap.get_cells(
			variant=[min(variants)],
			generation=np.arange(0, 1),
			only_successful=True)
		simOutDir = os.path.join(cell_paths_initialize[0], 'simOut')
		rna_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_ids_rna_counts_table = rna_counts_reader.readAttribute(
			'mRNA_cistron_ids')
		# Get indexes of mRNA cistrons in this subcolumn
		mRNA_cistron_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(mRNA_cistron_ids_rna_counts_table)}
		mRNA_cistron_indexes = np.array([
			mRNA_cistron_id_to_index[cistron_id] for cistron_id
			in mRNA_cistron_ids])
		# Get subcolumn for monomer IDs in monomer counts table
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, 'MonomerCounts'))
		monomer_ids_monomer_counts_table = monomer_counts_reader.readAttribute(
			'monomerIds')
		# Get indexes of monomers in this subcolumn
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids_monomer_counts_table)}
		monomer_indexes = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids])

		essential_genes = validation_data.essential_genes.essential_genes
		is_essential = np.array([True if gene_id in essential_genes else False for gene_id in gene_ids])
		total_genes = len(gene_ids)
		total_essential_genes = len(essential_genes)

		for variant in variants:
			print("Variant: ", variant)

			cell_paths = self.ap.get_cells(
				variant=[variant],
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
				only_successful=True)

			# Get boolean matrix for whether each gene's mRNA exists in each
			# generation or not
			mRNA_exists_in_gen = read_stacked_columns(
				cell_paths, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True, fun=lambda x: x.sum(axis=0) > 0)[:, mRNA_cistron_indexes]

			# Divide by total number of cells to get probability
			p_mRNA_exists_in_gen = (
				mRNA_exists_in_gen.sum(axis=0) / mRNA_exists_in_gen.shape[0])

			# Get maximum counts of mRNAs for each gene across all timepoints
			max_mRNA_counts = read_stacked_columns(
				cell_paths, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True).max(axis=0)[mRNA_cistron_indexes]

			# Get maximum counts of monomers for each gene across all timepoints
			max_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True).max(axis=0)[monomer_indexes]

			# Get indexes of subgenerationally expressed genes
			subgen_exp_indexes = np.where(np.logical_and(
				p_mRNA_exists_in_gen > 0, p_mRNA_exists_in_gen < 1
				))[0]

			# Get number of subgenerationally expressed genes
			n_subgen_genes = len(subgen_exp_indexes)

			# Count number of subgenerationally expressed genes that are essential
			n_subgen_essential_genes = np.sum(is_essential[subgen_exp_indexes])

			# Write data to table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '.tsv'), 'a') as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant', 'n_subgen_genes', 'n_subgen_essential_genes',
						'fraction_genes_subgen',
						'fraction_essential_genes_subgen'])

				writer.writerow([
					variant, n_subgen_genes, n_subgen_essential_genes,
					n_subgen_genes / total_genes,
					n_subgen_essential_genes / total_essential_genes])

			# Write gene data to detailed table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_detailed.tsv'), 'a') as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant','gene_name', 'cistron_name', 'protein_name',
						'p_expressed', 'max_mRNA_count', 'max_protein_count',
						'is_essential'])

				for i in subgen_exp_indexes:
					writer.writerow([
						variant, gene_ids[i], mRNA_cistron_ids[i],
						monomer_ids[i][:-3], p_mRNA_exists_in_gen[i],
						max_mRNA_counts[i], max_monomer_counts[i],
						is_essential[i]])

if __name__ == "__main__":
	Plot().cli()
