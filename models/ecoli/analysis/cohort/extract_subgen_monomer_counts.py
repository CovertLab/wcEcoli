"""
Template for cohort analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from matplotlib import cm

import csv

from wholecell.utils import units
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure, stacked_cell_identification,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

IGNORE_FIRST_N_GENS = 20
SEED_RANGE = np.arange(0, 40)


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
			# Ignore data from predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed = SEED_RANGE,
			only_successful=True)
		
		print('Analyzing %d cells...' % len(cell_paths))

		# There are 4346 mRNA ids with counts
		RNA_reader = TableReader(
				os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
		mRNA_ids = RNA_reader.readAttribute('mRNA_cistron_ids')
		RNA_reader.close()

		mRNA_id_to_index = {
			cistron_id: i for (i, cistron_id)
			in enumerate(mRNA_ids)
		}

		# There are 4539 mRNA ids total w/ gene names
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id']
			for cistron in sim_data.process.transcription.cistron_data
		}

		# There are 4310 mRNA ids with associated protein/monomer ids
		protein_id_to_cistron_id = {
			protein['id']:protein['cistron_id']
			for protein in sim_data.process.translation.monomer_data
		}

		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()


		# Get indexes of monomers in this subcolumn
		monomer_id_to_index = {
			monomer_id: i for (i, monomer_id)
			in enumerate(monomer_ids)
		}

		monomer_indices = np.array([
			monomer_id_to_index[monomer_id] for monomer_id in monomer_ids
		])

		# order cistrons in the order of monomers ids
		cistron_ids_in_order = np.array([
			protein_id_to_cistron_id[monomer_id] for monomer_id in monomer_ids
		])

		# order gene names in the order of monomers ids
		gene_ids_in_order = np.array([
			cistron_id_to_gene_id[cistron_id] for cistron_id in cistron_ids_in_order
		])

		# # Get indices of cistron_ids_in_order
		# mRNA_ids_indices = np.array([
		# 	mRNA_id_to_index[cistron_id] for cistron_id
		# 	in cistron_ids_in_order
		# ])

		# Monomer exists per gen
		monomer_exists_in_gen = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts',
			ignore_exception=True, fun=lambda x: x.sum(axis=0) > 0)[
							:, monomer_indices]
		print("monomer exists in gen shape: ")
		print(monomer_exists_in_gen.shape)
		# Divide by total number of cells to get probability
		p_monomer_exists_in_gen = (
			monomer_exists_in_gen.sum(axis=0) / monomer_exists_in_gen.shape[0]
			)
		subgenerational_monomer_mask = (
				(p_monomer_exists_in_gen > 0)
				& (p_monomer_exists_in_gen < 1)
		)

		sub_gen_p_monomer_exists_in_gen = np.array(p_monomer_exists_in_gen)[subgenerational_monomer_mask]

		print("sub_gen_p_monomer_exists_in_gen: ")
		print(sub_gen_p_monomer_exists_in_gen.shape)

		expression_status_array = np.full(
			p_monomer_exists_in_gen.shape,
			'always_expressed' , dtype='<U20')

		expression_status_array[subgenerational_monomer_mask] = 'subgen'
		expression_status_array[p_monomer_exists_in_gen == 0] = 'never_expressed'

		sub_gen_monomer_indices = np.array(monomer_indices)[subgenerational_monomer_mask]
		sub_gen_monomer_ids = np.array(monomer_ids)[subgenerational_monomer_mask]
		sub_gen_cistron_ids = np.array(cistron_ids_in_order)[subgenerational_monomer_mask]
		sub_gen_gene_ids = np.array(gene_ids_in_order)[subgenerational_monomer_mask]


		max_monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				ignore_exception=True).max(axis=0)[sub_gen_monomer_indices]
		print("max monomer counts shape: ")
		print(max_monomer_counts.shape)
		
		# import ipdb; ipdb.set_trace()

		# max_monomer_counts_all = np.zeros(monomer_indices.shape)
		# mean_monomer_counts_all = np.zeros(monomer_indices.shape)
		
		# max_mRNA_counts_all = np.zeros(mRNA_ids_indices.shape)
		# mean_mRNA_counts_all = np.zeros(mRNA_ids_indices.shape)

		# Use for loop to avoid memory issues (128 seeds and 32 gens = 4096 cells)

		# for i, monomer_index in enumerate(monomer_indices):
		# 	# Get maximum counts of monomers for each gene across all timepoints
		# 	max_monomer_counts = read_stacked_columns(
		# 		cell_paths, 'MonomerCounts', 'monomerCounts',
		# 		ignore_exception=True).max(axis=0)[monomer_index]
		# 	print("max monomer counts shape: ")
		# 	print(max_monomer_counts.shape)
		# 	mean_monomer_counts = read_stacked_columns(
		# 		cell_paths, 'MonomerCounts', 'monomerCounts',
		# 		ignore_exception=True).mean(axis=0)[monomer_index]
		# 	print("mean monomer counts shape: ")
		# 	print(mean_monomer_counts.shape)
						
		# 	# Get maximum counts of mRNAs for each gene across all timepoints
		# 	mRNA_index = mRNA_ids_indices[i]
		# 	max_mRNA_counts = read_stacked_columns(
		# 		cell_paths, 'RNACounts', 'mRNA_cistron_counts',
		# 		ignore_exception=True).max(axis=0)[mRNA_index]
		# 	print("max mRNA counts shape: ")
		# 	print(max_mRNA_counts.shape)
			
		# 	mean_mRNA_counts = read_stacked_columns(
		# 		cell_paths, 'RNACounts', 'mRNA_cistron_counts', 
		# 		ignore_exception=True).mean(axis=0)[mRNA_index]
		# 	print("mean mRNA counts shape: ")
		# 	print(mean_mRNA_counts.shape)

		# 	max_monomer_counts_all[i] = max_monomer_counts
		# 	mean_monomer_counts_all[i] = mean_monomer_counts
		# 	max_mRNA_counts_all[i] = max_mRNA_counts
		# 	mean_mRNA_counts_all[i] = mean_mRNA_counts


		# Write data to table
		with open(os.path.join(plotOutDir, plotOutFileName + '_subgen_only.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'gene_name', 'cistron_name', 'monomer_name',
				'prob_monomer_expressed', 
				# 'max_mRNA_count', 
				# 'mean_mRNA_count',
				'max_protein_count', 
				# 'mean_protein_count', 
				# 'expression_status'
			])

			for i,_ in enumerate(sub_gen_monomer_indices):
				writer.writerow([
					sub_gen_gene_ids[i], sub_gen_cistron_ids[i], sub_gen_monomer_ids[i],
					sub_gen_p_monomer_exists_in_gen[i], 
					# max_mRNA_counts_all[i], 
					# mean_mRNA_counts_all[i],
					max_monomer_counts[i],
					# mean_monomer_counts_all[i], 
					
				])


if __name__ == '__main__':
	Plot().cli()
