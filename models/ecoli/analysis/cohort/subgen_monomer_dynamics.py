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
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

IGNORE_FIRST_N_GENS = 8
SEEDS = np.arange(7, 10)

monomers_of_interest = ['GLYCDEH-MONOMER[c]',  # gldA
						'BETAGALACTOSID-MONOMER[c]',  # lacZ
						'RIBULOKIN-MONOMER[c]',  # araB
						'BAES-MONOMER[i]',  # baeS
						'G6504-MONOMER[o]',  # gfcE
						'EG11250-MONOMER[c]',  # chpS
						'EG11222-MONOMER[c]',  # alkA
						'G7263-MONOMER[c]',  # murQ
						'EG11249-MONOMER[c]',  # mazF const
						'EG10466-MONOMER[c]'  # hupA const
						]

monomers_of_interest_name_dict = {'GLYCDEH-MONOMER[c]': 'gldA',
						'BETAGALACTOSID-MONOMER[c]': 'lacZ',
						'RIBULOKIN-MONOMER[c]': 'araB',
						'BAES-MONOMER[i]': 'baeS',
						'G6504-MONOMER[o]': 'gfcE',
						'EG11250-MONOMER[c]': 'chpS',
						'EG11222-MONOMER[c]': 'alkA',
						'G7263-MONOMER[c]': 'murQ',
						'EG11249-MONOMER[c]': 'mazF',
						'EG10466-MONOMER[c]': 'hupA'
								  }


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
			monomer_id_to_index[monomer_id] for monomer_id in monomers_of_interest
		])


		# order cistrons in the order of monomers ids
		cistron_ids_in_order = np.array([
			protein_id_to_cistron_id[monomer_id] for monomer_id in monomers_of_interest
		])

		gene_names_in_order = np.array([
			monomers_of_interest_name_dict[monomer_id] for monomer_id in monomers_of_interest
		])

		# Get indices of cistron_ids_in_order
		mRNA_ids_indices = np.array([
			mRNA_id_to_index[cistron_id] for cistron_id
			in cistron_ids_in_order
		])

		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts')[:, monomer_indices]


		cistron_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts')[:, mRNA_ids_indices]


		def extract_doubling_times(cell_paths):
			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
			# Determine doubling time
			doubling_times = read_stacked_columns(cell_paths, 'Main', 'time', fun=lambda x: (x[-1] - x[0])).squeeze().astype(int)
			end_generation_times = np.cumsum(doubling_times) + time[0] #
			start_generation_indices = np.searchsorted(time, end_generation_times[:-1], side = 'left').astype(int)
			start_generation_indices = np.insert(start_generation_indices, 0, 0) + np.arange(len(doubling_times))
			end_generation_indices = start_generation_indices + doubling_times
			return time.astype(int), doubling_times, end_generation_times, start_generation_indices, end_generation_indices

		def plot_counts_dynamics(counts, color, molecules_of_interest, molecule_type, time, end_generation_times, seed):

			# Counts Plot
			num_groups = len(molecules_of_interest)
			cols = min(num_groups, 4)  # Number of columns for the grid
			rows = (num_groups + cols - 1) // cols  # Calculate the number of rows
			fig_width = cols * 10
			fig_height = fig_width / 2
			fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(fig_width, fig_height), sharex=False)

			if cols == 1:
				axes = [axes]
			axes = axes.flatten()

			# Flatten the axes array for easier iteration
			for i, ax in enumerate(axes):
				if i < len(molecules_of_interest):  # Check if there's a molecule for this subplot
					molecule_id = molecules_of_interest[i]

					interest_counts = counts[time, i]

					ax.plot(time / 60, interest_counts, color=color, linewidth=6,
							label=molecule_id)

					ax.set_xlabel('Time (min)', fontsize=30);
					ax.set_ylabel(f'{molecule_type} counts', fontsize=30)
					ax.set_title(molecule_id, fontsize=30)
					ax.tick_params(axis='x', labelsize=30)
					ax.tick_params(axis='y', labelsize=30)

					ax.spines['right'].set_visible(False)
					ax.spines['top'].set_visible(False)

					for x in end_generation_times / 60:
						ax.axvline(x=x,
								   color='grey',
								   linestyle='dashed')
			# Remove any empty subplots
			for i in range(num_groups, rows * cols):
				fig.delaxes(axes.flat[i])

			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + f'_{molecule_type}_dynamics_{seed}', metadata)


		for seed in SEEDS:
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
				only_successful=True)

			if not np.all([self.ap.get_successful(cell) for cell in cell_paths_per_seed]):
				continue

			time, _, end_generation_times, _, _ = extract_doubling_times(
				cell_paths_per_seed)

			if monomers_of_interest is not None:

				plot_counts_dynamics(monomer_counts, 'mediumseagreen', gene_names_in_order, 'monomer', time,
									 end_generation_times, seed)


				plot_counts_dynamics(cistron_counts, 'mediumseagreen', gene_names_in_order, 'mRNA', time,
									 end_generation_times, seed)




if __name__ == '__main__':
	Plot().cli()
