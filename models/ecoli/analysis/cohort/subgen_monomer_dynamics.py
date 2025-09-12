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
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

IGNORE_FIRST_N_GENS = 8

monomers_of_interest = ['GLYCDEH-MONOMER[c]',  # gldA
						'BETAGALACTOSID-MONOMER[c]',  # lacZ
						'RIBULOKIN-MONOMER[c]',  # araB
						'BAES-MONOMER[i]',  # baeS
						'G6504-MONOMER[m]',  # gfcE
						'EG11250-MONOMER[c]',  # chpS
						'EG11222-MONOMER[c]',  # alkA
						'G7263-MONOMER[c]',  # murQ
						'EG11249-MONOMER[c]',  # mazF const
						'EG10466-MONOMER[c]'  # hupA const
						]


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
			doubling_times = read_stacked_columns(cell_paths, 'Main', 'time',
												  fun=lambda x: (x[-1] - x[0])).squeeze().astype(int)
			end_generation_times = np.cumsum(doubling_times) + time[0]  #
			start_generation_indices = np.searchsorted(time, end_generation_times[:-1], side='left').astype(int)
			start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(len(doubling_times))
			end_generation_indices = start_generation_indices + doubling_times
			return time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices

		def plot_monomer_complex_dynamics(monomer_counts, monomers_of_interest,
										  time, end_generation_times,
										  seed):

			# Protein Counts Plot
			num_groups = len(monomers_of_interest)
			cols = min(num_groups, 4)  # Number of columns for the grid
			rows = (num_groups + cols - 1) // cols  # Calculate the number of rows
			fig_width = cols * 10
			fig_height = fig_width / 2
			fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(fig_width, fig_height), sharex=False)

			if cols == 1:
				axes = [axes]

			# Flatten the axes array for easier iteration
			for i, ax in enumerate(axes):
				if i < len(monomers_of_interest):  # Check if there's a monomer for this subplot
					monomer_id = monomers_of_interest[i]
					monomer_index = monomer_id_to_index[monomer_id]
					interest_monomer_counts = monomer_counts[:, monomer_index]
					ax.plot(time / 60, interest_monomer_counts, color='royalblue', linewidth=6,
							label=monomer_id)

					ax.set_xlabel("Time (min)", fontsize=20);
					ax.set_ylabel("Monomer counts", fontsize=20)
					ax.set_title(monomer_id)
					ax.tick_params(axis='x', labelsize=20)
					ax.tick_params(axis='y', labelsize=20)
					ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
					for x in end_generation_times / 60:
						ax.axvline(x=x,
								   color='grey',
								   linestyle='dashed')
			# Remove any empty subplots
			for i in range(num_groups, rows * cols):
				fig.delaxes(axes.flat[i])

			sns.despine();
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + f'_monomer_dynamics_{seed}', metadata)

		for seed in self.ap.get_seeds():
			cell_paths = self.ap.get_cells(
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
				only_successful=True)

			if not np.all([self.ap.get_successful(cell) for cell in cell_paths]):
				continue

			time, doubling_times, end_generation_times, start_generation_indices, end_generation_indices = extract_doubling_times(
				cell_paths)

			if monomers_of_interest is not None:
				plot_monomer_complex_dynamics(monomer_counts, monomers_of_interest,
											  time, end_generation_times,
											  seed)

		extract_doubling_times(cell_paths)

		import ipdb;
		ipdb.set_trace()


if __name__ == '__main__':
	Plot().cli()
