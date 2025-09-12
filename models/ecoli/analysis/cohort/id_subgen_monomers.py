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

		# Extract all monomer and cistron ids in simulation
		monomer_ids = sim_data.process.translation.monomer_data['id']
		cistron_ids = sim_data.process.transcription.cistron_data['id']

		# Filter list of cistron IDs to only keep ones with associated monomer ids
		monomer_id_to_cistron_id = {
			monomer['id']: monomer['cistron_id']
			for monomer in monomer_ids
		}

		mRNA_cistron_ids = [
			cistron_id for cistron_id in cistron_ids
			if cistron_id in monomer_id_to_cistron_id.values()]


		# Get subcolumn for mRNA cistron IDs in RNA counts table to extract cistron indices
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)

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

		# Get subcolumn for monomer IDs in monomer counts table to extract monomer indices
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

		cistron_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts')[:, cistron_indexes]


		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts')[:, monomer_indexes]

		# Get maximum counts of monomers for each gene across all timepoints
		max_monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts',
			ignore_exception=True).max(axis=0)[monomer_indexes]

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


		def subgen_monomer_status(monomer_counts, end_generation_indices, doubling_times):

			monomer_expressed_bool = (monomer_counts > 0) * 1

			monomer_expressed_bool_by_generation = np.split(monomer_expressed_bool, end_generation_indices, axis=0)

			monomer_expressed_by_generation = np.array([
				np.sum(monomer_expressed_bool_by_generation[i], axis=0)
				for i in range(len(doubling_times))
			])

			monomer_expressed_per_generation_bool = (monomer_expressed_by_generation > 0) * 1

			frequency_of_monomer_over_all_generations = np.sum(monomer_expressed_per_generation_bool, axis=0) / len(
				cell_cycle_length)

			subgenerational_monomer_mask = (
					(frequency_of_monomer_over_all_generations > 0)
					& (frequency_of_monomer_over_all_generations < 1)
			)

			subgenerational_monomer_ids = np.array(monomer_ids)[subgenerational_monomer_mask]
			not_expressed_monomer_ids = np.array(monomer_ids)[frequency_of_monomer_over_all_generations == 0]
			always_expressed_monomer_ids = np.array(monomer_ids)[frequency_of_monomer_over_all_generations == 1]

			return subgenerational_monomer_ids, not_expressed_monomer_ids, always_expressed_monomer_ids


		_, doubling_times, _, _, end_generation_indices = extract_doubling_times(cell_paths)

		subgenerational_monomer_ids, not_expressed_monomer_ids, always_expressed_monomer_ids = subgen_monomer_status(monomer_counts, end_generation_indices, doubling_times)

		import ipdb;
		ipdb.set_trace()

if __name__ == '__main__':
	Plot().cli()
