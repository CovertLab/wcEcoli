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
import pandas as pd
import seaborn as sns

from wholecell.utils import units
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure, stacked_cell_identification,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

IGNORE_FIRST_N_GENS = 8


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
			# Ignore data from predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)

		# There are 4346 mRNA ids with counts
		RNA_reader = TableReader(
				os.path.join(cell_paths[0], 'simOut', 'RNACounts'))
		mRNA_ids = RNA_reader.readAttribute('mRNA_cistron_ids')

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

		# Get indices of cistron_ids_in_order
		mRNA_ids_indices = np.array([
			mRNA_id_to_index[cistron_id] for cistron_id
			in cistron_ids_in_order
		])

		# Get maximum counts of monomers for each gene across all timepoints
		max_monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts',
			ignore_exception=True).max(axis=0)[monomer_indices]

		mean_monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts',
			ignore_exception=True).mean(axis=0)[monomer_indices]

		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts')[:, monomer_indices]

		# Get maximum counts of mRNAs for each gene across all timepoints
		max_mRNA_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts',
			ignore_exception=True).max(axis=0)[mRNA_ids_indices]

		mean_mRNA_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts',
			ignore_exception=True).mean(axis=0)[mRNA_ids_indices]

		def extract_doubling_times(cell_paths):
			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
			# Determine doubling time
			doubling_times = read_stacked_columns(cell_paths, 'Main', 'time', fun=lambda x: (x[-1] - x[0])).squeeze().astype(int)
			end_generation_times = np.cumsum(doubling_times) + time[0] #
			start_generation_indices = np.searchsorted(time, end_generation_times[:-1], side = 'left').astype(int)
			start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(len(doubling_times))
			end_generation_indices = start_generation_indices + doubling_times
			return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices


		def subgen_monomer_status_per_seed(monomer_counts, end_generation_indices, doubling_times):

			monomer_expressed_bool = (monomer_counts > 0) * 1

			monomer_expressed_bool_by_generation = np.split(monomer_expressed_bool, end_generation_indices, axis=0)

			monomer_expressed_by_generation = np.array([
				np.sum(monomer_expressed_bool_by_generation[i], axis=0)
				for i in range(len(doubling_times))
			])

			monomer_expressed_per_generation_bool = (monomer_expressed_by_generation > 0) * 1

			return monomer_expressed_per_generation_bool

		subgen_matrix_list = []

		for seed in self.ap.get_seeds():
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
				only_successful=True)

			if not np.all([self.ap.get_successful(cell) for cell in cell_paths_per_seed]):
				continue

			_, doubling_times, _, _, end_generation_indices = extract_doubling_times(
				cell_paths_per_seed)
			monomer_expressed_per_generation_bool = subgen_monomer_status_per_seed(monomer_counts, end_generation_indices, doubling_times)
			subgen_matrix_list.append(monomer_expressed_per_generation_bool)

		all_seeds_subgen_status_bool = np.vstack(subgen_matrix_list)

		frequency_of_monomer_over_all_generations = np.sum(all_seeds_subgen_status_bool, axis=0) / len(
			all_seeds_subgen_status_bool)

		subgenerational_monomer_mask = (
				(frequency_of_monomer_over_all_generations > 0)
				& (frequency_of_monomer_over_all_generations < 1)
		)

		expression_status_array = np.full(
			frequency_of_monomer_over_all_generations.shape,
			'always_expressed' , dtype='<U20')

		expression_status_array[subgenerational_monomer_mask] = 'subgen'
		expression_status_array[frequency_of_monomer_over_all_generations == 0] = 'never_expressed'

		# Write data to table
		with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'gene_name', 'cistron_name', 'monomer_name',
				'prob_monomer_expressed', 'max_mRNA_count', 'mean_mRNA_count',
				'max_protein_count', 'mean_protein_count', 'expression_status'
			])

			for i in monomer_indices:
				writer.writerow([
					gene_ids_in_order[i], cistron_ids_in_order[i], monomer_ids[i],
					frequency_of_monomer_over_all_generations[i], max_mRNA_counts[i], mean_mRNA_counts[i],
					max_monomer_counts[i], mean_monomer_counts[i], expression_status_array[i]
				])


if __name__ == '__main__':
	Plot().cli()
