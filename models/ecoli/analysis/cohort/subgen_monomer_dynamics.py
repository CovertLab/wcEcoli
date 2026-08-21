"""
Per-generation mRNA and monomer count traces for selected subgen and anaerobic genes

x is normalized so each generation spans unit width:
  x(t) = gen_number + (t - t_start) / (t_end - t_start)

Writes one figure pair (monomer + mRNA) per plotted seed:
  <plotOutFileName>_{monomer,mRNA}_dynamics_seed<SSSSSS>_<color>.pdf
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
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.analysis.analysis_tools import (exportFigure, stacked_cell_identification,
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS

# not so random subset of seeds 
SEEDS = np.arange(20, 25)
COLOR_LINE = 'mediumseagreen'  # 'skyblue'

#   'anaerobic'  -- anaerobic-respiration complexes, one subunit each
#   'curated10'  -- the 10-gene panel shared by subgen_peak_counts.py /
#                   protein_distribution.py
PANEL = 'anaerobic'
monomers_of_interest, monomers_of_interest_name_dict = sc.curated_panel(PANEL)


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
			# Ignore data from predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed = SEEDS,
			only_successful=True)
		
		cell_paths = self.ap.get_cells(
        generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), 
        seed=SEEDS,
        only_successful=True)

		if len(cell_paths) == 0 :
			return
		

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


		def extract_doubling_times(cell_paths):
			# Load data
			time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
			# Determine doubling time
			doubling_times = read_stacked_columns(cell_paths, 'Main', 'time', fun=lambda x: (x[-1] - x[0])).squeeze().astype(int)
			end_generation_times = np.cumsum(doubling_times) + time[0] #
			# Per-generation row boundaries from actual row counts. read_stacked_columns
			# stacks one cell (= one generation) per block, so labeling each row by its
			# source cell segments generations exactly, with no assumption that a row is
			# 1 s and no boundary-duplicate drift.
			cell_ids = stacked_cell_identification(cell_paths, 'Main', 'time').squeeze().astype(int)
			n_cells = len(cell_paths)
			start_generation_indices = np.searchsorted(cell_ids, np.arange(n_cells), side='left')
			# Inclusive last row of each generation; callers slice time[start:end + 1].
			end_generation_indices = np.append(start_generation_indices[1:], len(cell_ids)) - 1
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

					interest_counts = counts[:, i]

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
			exportFigure(plt, plotOutDir, plotOutFileName + f'_{molecule_type}_dynamics_{seed}_{color}', metadata)
			
		def plot_counts_dynamics_gen(counts, color, molecules_of_interest, molecule_type, time, start_generation_indices,
                         end_generation_indices, end_generation_times, seed):
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
			

			for i, ax in enumerate(axes):
				if i < len(molecules_of_interest):  # Check if there's a molecule for this subplot
					molecule_id = molecules_of_interest[i]
					interest_counts = counts[:, i]

					gen_number = 0

					# Iterate over each generation and plot the segment separately
					for start_idx, end_idx in zip(start_generation_indices, end_generation_indices):
						
						# Select the time and counts data for the current generation
						# Note: We use end_idx + 1 to include the last point of the generation
						gen_time = time[start_idx : end_idx + 1]
						gen_counts = interest_counts[start_idx : end_idx + 1]
						time_start = gen_time[0]
						time_end = gen_time[-1]
						generation_duration = time_end - time_start
						normalized_gen_time = (gen_time - time_start)/generation_duration
						x_axis_data = gen_number + normalized_gen_time

						# Plot the segment. The use of multiple ax.plot calls creates disconnected lines.
						ax.plot(x_axis_data, gen_counts, color=color, linewidth=6, label=molecule_id if start_idx == 0 else None)
						gen_number += 1

					ax.set_xlabel('Generation number', fontsize=30);
					ax.set_ylabel(f'{molecule_type} counts', fontsize=30)
					ax.set_title(molecule_id, fontsize=30)
					ax.tick_params(axis='x', labelsize=30)
					ax.tick_params(axis='y', labelsize=30)
					max_gen = len(start_generation_indices)
					ax.set_xticks(range(max_gen + 1))
					ax.set_xlim(0, max_gen) # Set limit from 0 to the total number of generations

					ax.spines['right'].set_visible(False)
					ax.spines['top'].set_visible(False)

					# Generation boundaries are already the integer x ticks above:
					# x is normalized to generation number, so every tick IS a
					# boundary. The old code drew axvlines at
					# `end_generation_times / 60` -- minutes on a generation-number
					# axis, which put them far off the right of the plot (or all at
					# ~0 early on). subgen_monomer_dynamics_def5.py omits them for
					# the same reason.

			# Remove any empty subplots
			for i in range(num_groups, rows * cols):
				fig.delaxes(axes.flat[i])

			plt.tight_layout()
			# `_seed%06d`, matching subgenerationalTranscription_def5.py. The old
			# infix was `_gen_{seed}`, which read as a generation number when it has
			# always been the seed.
			exportFigure(plt, plotOutDir,
				plotOutFileName
					+ f'_{molecule_type}_dynamics_seed{seed:06d}_{color}',
				metadata)
			

		# Strict-successful lineages (completed every generation and no cell at
		# the 180-min doubling cap).
		success = sc.compute_lineage_success(self.ap, self.ap.n_generation)
		plotted = [s for s in SEEDS if s in success['successful_seeds']]
		print('Plotting %d of %d requested seeds (%s); the rest are not '
			'strict-successful lineages. This is the only reason the figure count '
			'differs between cohorts.'
			% (len(plotted), len(SEEDS),
				', '.join(str(s) for s in plotted) or 'none'))

		for seed in SEEDS:
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation), seed=[seed],
				only_successful=True)

			if seed not in success['successful_seeds']:
				continue
			
			monomer_counts = read_stacked_columns(cell_paths_per_seed, 'MonomerCounts', 'monomerCounts')[:, monomer_indices]
			cistron_counts = read_stacked_columns(cell_paths_per_seed, 'RNACounts', 'mRNA_cistron_counts')[:, mRNA_ids_indices]

			time, doubling_times, end_generation_times,start_generation_indices, end_generation_indices = extract_doubling_times(
				cell_paths_per_seed)

			if monomers_of_interest is not None:

				plot_counts_dynamics_gen(monomer_counts, COLOR_LINE, gene_names_in_order, 'monomer', time, start_generation_indices,
											end_generation_indices, end_generation_times, seed)

				plot_counts_dynamics_gen(cistron_counts, COLOR_LINE, gene_names_in_order, 'mRNA', time, start_generation_indices,
										end_generation_indices, end_generation_times, seed)




if __name__ == '__main__':
	Plot().cli()