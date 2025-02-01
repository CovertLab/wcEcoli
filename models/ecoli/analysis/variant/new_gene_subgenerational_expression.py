"""
Generates a tables of genes that are subgenerationally and not expressed at the
mRNA and monomer level, with their expression frequencies and average/maximum
mRNA/protein counts for each variant index.
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

			if variant == min(variants):
				file_mode = "w"
			else:
				file_mode = "a"

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

			# Get indexes of subgenerationally expressed genes
			subgen_exp_indexes = np.where(np.logical_and(
				p_mRNA_exists_in_gen > 0, p_mRNA_exists_in_gen < 1
				))[0]

			# Get number of subgenerationally expressed genes
			n_subgen_genes = len(subgen_exp_indexes)

			# Count number of subgenerationally expressed genes that are essential
			n_subgen_essential_genes = np.sum(is_essential[subgen_exp_indexes])

			# Get indexes of subgenerationally expressed genes
			zero_exp_indexes = np.where(p_mRNA_exists_in_gen == 0)[0]

			# Get number of 0 expressed genes
			n_zero_genes = len(zero_exp_indexes)

			# Count number of subgenerationally expressed genes that are essential
			n_zero_essential_genes = np.sum(is_essential[zero_exp_indexes])

			# Get maximum counts of monomers for each gene across all timepoints
			max_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True).max(axis=0)[monomer_indexes]

			# Get boolean matrix for whether each gene's monomer exists in each
			# generation or not
			monomer_exists_in_gen = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True, fun=lambda x: x.sum(axis=0) > 0)[:, monomer_indexes]

			# Divide by total number of cells to get probability
			p_monomer_exists_in_gen = (
					monomer_exists_in_gen.sum(axis=0) / monomer_exists_in_gen.shape[0])

			# Get indexes of genes that are subgenerational at the monomer level
			monomer_subgen_exp_indexes = np.where(np.logical_and(
				p_monomer_exists_in_gen > 0, p_monomer_exists_in_gen < 1
				))[0]

			# Get number of genes that are subgenerational at the monomer level
			n_monomer_subgen_genes = len(monomer_subgen_exp_indexes)

			# Count number of essential genes that are subgenerational at the monomer level
			n_monomer_subgen_essential_genes = np.sum(is_essential[monomer_subgen_exp_indexes])

			# Get indexes of zero expressed genes
			monomer_zero_exp_indexes = np.where(p_monomer_exists_in_gen == 0)[0]

			# Get number of 0 expressed genes
			n_monomer_zero_genes = len(monomer_zero_exp_indexes)

			# Count number of subgenerationally expressed genes that are essential
			n_monomer_zero_essential_genes = np.sum(is_essential[monomer_zero_exp_indexes])

			# Write data to table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant', 'n_subgen_genes', 'n_subgen_essential_genes',
						'fraction_genes_subgen',
						'fraction_essential_genes_subgen',
						'n_zero_genes', 'n_zero_essential_genes',
						'fraction_genes_zero',
						'fraction_essential_genes_zero'
						])

				writer.writerow([
					variant, n_subgen_genes, n_subgen_essential_genes,
					n_subgen_genes / total_genes,
					n_subgen_essential_genes / total_essential_genes,
					n_zero_genes, n_zero_essential_genes,
					n_zero_genes / total_genes,
					n_zero_essential_genes / total_essential_genes])

			# Write gene data to detailed table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_detailed.tsv'), file_mode) as f:
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

			# Write gene zero expression data to detailed table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_zero_detailed.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant','gene_name', 'cistron_name', 'protein_name',
						'p_expressed', 'max_mRNA_count', 'max_protein_count',
						'is_essential'])

				for i in zero_exp_indexes:
					writer.writerow([
						variant, gene_ids[i], mRNA_cistron_ids[i],
						monomer_ids[i][:-3], p_mRNA_exists_in_gen[i],
						max_mRNA_counts[i], max_monomer_counts[i],
						is_essential[i]])

			# Write monomer data to table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_monomer.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant', 'n_monomer_subgen_genes', 'n_monomer_subgen_essential_genes',
						'fraction_genes_monomer_subgen',
						'fraction_essential_genes_monomer_subgen',
						'n_monomer_zero_genes', 'n_monomer_zero_essential_genes',
						'fraction_genes_zero_subgen',
						'fraction_essential_genes_zero_subgen'])

				writer.writerow([
					variant, n_monomer_subgen_genes, n_monomer_subgen_essential_genes,
					n_monomer_subgen_genes / total_genes,
					n_monomer_subgen_essential_genes / total_essential_genes,
					n_monomer_zero_genes, n_monomer_zero_essential_genes,
					n_monomer_zero_genes / total_genes,
					n_monomer_zero_essential_genes / total_essential_genes])

			# Write monomer data to detailed table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_monomer_detailed.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant','gene_name', 'cistron_name', 'protein_name',
						'p_mRNA_expressed', 'p_monomer_expressed',
						'max_mRNA_count', 'max_protein_count',
						'is_essential'])

				for i in monomer_subgen_exp_indexes:
					writer.writerow([
						variant, gene_ids[i], mRNA_cistron_ids[i],
						monomer_ids[i][:-3], p_mRNA_exists_in_gen[i],
						p_monomer_exists_in_gen[i], max_mRNA_counts[i],
						max_monomer_counts[i], is_essential[i]])

			# Write zero monomer data to detailed table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_zero_monomer_detailed.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant', 'gene_name', 'cistron_name', 'protein_name',
						'p_mRNA_expressed', 'p_monomer_expressed',
						'max_mRNA_count', 'max_protein_count',
						'is_essential'])

				for i in monomer_zero_exp_indexes:
					writer.writerow([
						variant, gene_ids[i], mRNA_cistron_ids[i],
						monomer_ids[i][:-3], p_mRNA_exists_in_gen[i],
						p_monomer_exists_in_gen[i], max_mRNA_counts[i],
						max_monomer_counts[i], is_essential[i]])

			# Get list of monomers whose counts were 0 at any time point
			monomer_ever_nonexistent_by_gen = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True, fun=lambda x: np.count_nonzero(x == 0, axis = 0))[:, monomer_indexes]
			monomer_ever_nonexistent = np.sum(monomer_ever_nonexistent_by_gen, axis = 0) > 0

			num_time_steps_monomer_not_present_for = np.sum(monomer_ever_nonexistent_by_gen, axis = 0)
			num_gens_monomer_ever_nonexistent = np.sum(monomer_ever_nonexistent_by_gen > 0, axis = 0)

			# Get indexes of genes that are ever nonexistent at the monomer level
			monomer_ever_nonexistent_indexes = np.where(monomer_ever_nonexistent > 0)[0]

			# Get number of genes that are ever nonexistent at the monomer level
			n_monomer_ever_nonexistent_genes = len(monomer_ever_nonexistent_indexes)

			# Count number of essential genes that are ever nonexistent at the monomer level
			n_monomer_ever_nonexistent_essential_genes = np.sum(is_essential[monomer_ever_nonexistent_indexes])

			# Write monomer data to table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_monomer_ever_nonexistent.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant', 'n_monomer_ever_nonexistent_genes',
						'n_monomer_ever_nonexistent_essential_genes',
						'fraction_genes_monomer_ever_nonexistent',
						'fraction_essential_genes_monomer_ever_nonexistent',])

				writer.writerow([
					variant, n_monomer_ever_nonexistent_genes, n_monomer_ever_nonexistent_essential_genes,
					n_monomer_ever_nonexistent_genes / total_genes,
					n_monomer_ever_nonexistent_essential_genes / total_essential_genes])

			# Write ever non-existent monomer data to detailed table for this variant
			with open(os.path.join(
					plotOutDir, plotOutFileName + '_monomer_ever_nonexistent_detailed.tsv'), file_mode) as f:
				writer = csv.writer(f, delimiter='\t')
				if variant == min(variants):
					writer.writerow([
						'variant','gene_name', 'cistron_name', 'protein_name',
						'p_mRNA_expressed', 'p_monomer_expressed',
						'max_mRNA_count', 'max_protein_count',
						'num_time_steps_monomer_not_present_for',
						'num_gens_monomer_ever_nonexistent',
						'is_essential'])

				for i in monomer_ever_nonexistent_indexes:
					writer.writerow([
						variant, gene_ids[i], mRNA_cistron_ids[i],
						monomer_ids[i][:-3], p_mRNA_exists_in_gen[i],
						p_monomer_exists_in_gen[i], max_mRNA_counts[i],
						max_monomer_counts[i],
						num_time_steps_monomer_not_present_for[i],
						num_gens_monomer_ever_nonexistent[i], is_essential[i]])


if __name__ == "__main__":
	Plot().cli()
