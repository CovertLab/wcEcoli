"""
Per-generation mRNA and monomer count traces, for two gene sets.

x is normalized so each generation spans unit width:
  x(t) = gen_number + (t - t_start) / (t_end - t_start)

Both gene sets are plotted in one run, and the set name is in every output
filename:

  anaerobic     the fixed anaerobic-respiration panel, one subunit per complex
                (sc.curated_panel('anaerobic')), over SEEDS_ANAEROBIC.
  subgen_top    the N_SUBGEN_TO_PLOT subgenerational genes with the highest
                Definition-5 mean, taken from the extraction cache, over the
                first N_SEEDS_TO_PLOT strict-successful seeds.

Writes one figure pair (monomer + mRNA) per gene set per plotted seed:
  <plotOutFileName>_{monomer,mRNA}_<set>_seed<SSSSSS>.pdf

Definition 5 always means the CI form (sc.classify_def5_ci). The `subgen_top`
set needs subgen_extract.py to have run; if its cache is missing that set is
skipped with a message and the `anaerobic` figures are still produced.

Time-series traces always come from simOut, since dynamics need the full
within-generation trace; only the gene SELECTION for `subgen_top` is cached.
"""

import os
import pickle

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_helper_functions as sc
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, stacked_cell_identification)
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS

# not so random subset of seeds, for the fixed panel
SEEDS_ANAEROBIC = np.arange(20, 25)
# How many subgen genes to plot (highest Def-5 mean first), and over how many seeds
N_SUBGEN_TO_PLOT = 12
N_SEEDS_TO_PLOT = 3
COLOR_LINE = 'mediumseagreen'  # 'skyblue'

#   'anaerobic'  -- anaerobic-respiration complexes, one subunit each
#   'curated10'  -- the 10-gene panel shared by subgen_peak_counts.py /
#                   subgen_protein_distribution.py
PANEL = 'anaerobic'


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		# Ignore data from predefined number of generations per seed
		if self.ap.n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return

		# Strict-successful seeds
		success = sc.compute_seed_success(self.ap, self.ap.n_generation)
		successful = success['successful_seeds']
		if not len(successful):
			print('No strict-successful seeds found. Skipping.')
			return

		gene_sets = [self._anaerobic_gene_set(sim_data, successful)]
		subgen_set = self._subgen_gene_set(plotOutDir, successful)
		if subgen_set is not None:
			gene_sets.append(subgen_set)

		for set_name, names, cistron_ids, monomer_ids, seeds in gene_sets:
			if not len(seeds):
				print('[%s] No strict-successful seeds to plot. Skipping.'
					% set_name)
				continue
			print('[%s] %d genes over seeds %s: %s'
				% (set_name, len(names),
					', '.join(str(s) for s in seeds), ', '.join(names)))
			self._plot_gene_set(plotOutDir, plotOutFileName, metadata,
				set_name, names, cistron_ids, monomer_ids, seeds)

	# gene sets

	def _anaerobic_gene_set(self, sim_data, successful):
		"""The fixed curated panel, ordered by monomer id."""
		monomers_of_interest, name_dict = sc.curated_panel(PANEL)
		protein_id_to_cistron_id = {
			protein['id']: protein['cistron_id']
			for protein in sim_data.process.translation.monomer_data
			}
		# order cistrons in the order of monomers ids
		cistron_ids = [
			protein_id_to_cistron_id[m] for m in monomers_of_interest]
		names = [name_dict[m] for m in monomers_of_interest]
		seeds = [int(s) for s in SEEDS_ANAEROBIC if s in successful]
		if len(seeds) < len(SEEDS_ANAEROBIC):
			print('[%s] Plotting %d of %d requested seeds; the rest are not '
				'strict-successful seeds. This is the only reason the figure '
				'count differs between cohorts.'
				% (PANEL, len(seeds), len(SEEDS_ANAEROBIC)))
		return (PANEL, names, cistron_ids, list(monomers_of_interest), seeds)

	def _subgen_gene_set(self, plotOutDir, successful):
		"""The top-N def5_CI subgen genes, or None if the extraction is absent.

		This is the only part of the script that needs subgen_extract.py; the
		panel above stands on its own, so a missing cache costs half the
		figures rather than all of them.
		"""
		try:
			clf = sc.canonical_def5_classification(plotOutDir)
		except FileNotFoundError:
			print('[subgen_top] No extraction in %s, so the subgen gene set is '
				'skipped. Run subgen_extract.py to include it.' % plotOutDir)
			return None
		if clf['n_seeds'] == 0:
			print('[subgen_top] No successful seeds in the extraction. Skipping.')
			return None
		gene_ids, cistron_ids, monomer_ids = sc.load_raw_genes(plotOutDir)
		cat = clf['stats']['cat']
		mean = clf['stats']['mean']

		# Select the N subgen genes with the highest Def-5 mean.
		subgen_idx = np.where(cat == 'subgen')[0]
		if subgen_idx.size == 0:
			print('[subgen_top] No subgenerational genes under Definition 5. '
				'Skipping.')
			return None
		subgen_idx = subgen_idx[np.argsort(mean[subgen_idx])[::-1]]
		chosen = subgen_idx[:N_SUBGEN_TO_PLOT]
		seeds = sorted(int(s) for s in successful)[:N_SEEDS_TO_PLOT]
		return ('subgen_top',
			[gene_ids[i] for i in chosen],
			[cistron_ids[i] for i in chosen],
			[monomer_ids[i] for i in chosen],
			seeds)

	# plotting

	def _plot_gene_set(self, plotOutDir, plotOutFileName, metadata, set_name,
			names, cistron_ids, monomer_ids, seeds):
		probe = self.ap.get_cells(
			seed=[seeds[0]],
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(probe) == 0:
			print('[%s] No cells for probe seed %d. Skipping.'
				% (set_name, seeds[0]))
			return

		# Map the chosen genes into the listener subcolumns (from a real cell).
		mRNA_indices, monomer_indices = self._index_maps(
			probe[0], cistron_ids, monomer_ids)

		for seed in seeds:
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(
					IGNORE_FIRST_N_GENS, self.ap.n_generation),
				seed=[seed], only_successful=True)
			if len(cell_paths_per_seed) == 0:
				continue

			monomer_counts = read_stacked_columns(
				cell_paths_per_seed, 'MonomerCounts', 'monomerCounts'
				)[:, monomer_indices]
			cistron_counts = read_stacked_columns(
				cell_paths_per_seed, 'RNACounts', 'mRNA_cistron_counts'
				)[:, mRNA_indices]
			(time, start_generation_indices,
				end_generation_indices) = self._extract_doubling_times(
				cell_paths_per_seed)

			for counts, molecule_type in (
					(monomer_counts, 'monomer'), (cistron_counts, 'mRNA')):
				self._plot_counts_dynamics_gen(
					plotOutDir, plotOutFileName, metadata, set_name, counts,
					names, molecule_type, time, start_generation_indices,
					end_generation_indices, seed)

	def _index_maps(self, cell_path, cistron_ids, monomer_ids):
		"""Column indices for these genes in RNACounts and MonomerCounts."""
		rna_reader = TableReader(os.path.join(cell_path, 'simOut', 'RNACounts'))
		mRNA_ids = rna_reader.readAttribute('mRNA_cistron_ids')
		rna_reader.close()
		mRNA_id_to_index = {c: i for i, c in enumerate(mRNA_ids)}

		monomer_reader = TableReader(
			os.path.join(cell_path, 'simOut', 'MonomerCounts'))
		monomer_ids_table = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()
		monomer_id_to_index = {m: i for i, m in enumerate(monomer_ids_table)}

		return (np.array([mRNA_id_to_index[c] for c in cistron_ids]),
			np.array([monomer_id_to_index[m] for m in monomer_ids]))

	def _extract_doubling_times(self, cell_paths):
		time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
		# Per-generation row boundaries from actual row counts.
		# read_stacked_columns stacks one cell (= one generation) per block, so
		# labeling each row by its source cell segments generations exactly,
		# with no assumption that a row is 1 s and no boundary-duplicate drift.
		cell_ids = stacked_cell_identification(
			cell_paths, 'Main', 'time').squeeze().astype(int)
		n_cells = len(cell_paths)
		start_generation_indices = np.searchsorted(
			cell_ids, np.arange(n_cells), side='left')
		# Inclusive last row of each generation; callers slice time[start:end + 1].
		end_generation_indices = np.append(
			start_generation_indices[1:], len(cell_ids)) - 1
		return (time.astype(int), start_generation_indices,
			end_generation_indices)

	def _plot_counts_dynamics_gen(self, plotOutDir, plotOutFileName, metadata,
			set_name, counts, names, molecule_type, time,
			start_generation_indices, end_generation_indices, seed):
		num_groups = len(names)
		cols = min(num_groups, 4)  # Number of columns for the grid
		rows = (num_groups + cols - 1) // cols  # Calculate the number of rows
		fig_width = cols * 10
		fig_height = fig_width / 2
		fig, axes = plt.subplots(
			nrows=rows, ncols=cols, figsize=(fig_width, fig_height),
			sharex=False)
		axes = np.atleast_1d(axes).flatten()

		for i, ax in enumerate(axes):
			if i >= num_groups:
				continue
			interest_counts = counts[:, i]
			gen_number = 0
			# Iterate over each generation and plot the segment separately, so
			# the lines stay disconnected across divisions.
			for start_idx, end_idx in zip(
					start_generation_indices, end_generation_indices):
				# end_idx + 1 to include the last point of the generation
				gen_time = time[start_idx:end_idx + 1]
				gen_counts = interest_counts[start_idx:end_idx + 1]
				if len(gen_time) < 2:
					gen_number += 1
					continue
				duration = gen_time[-1] - gen_time[0]
				normalized = (gen_time - gen_time[0]) / duration
				ax.plot(gen_number + normalized, gen_counts,
					color=COLOR_LINE, linewidth=6)
				gen_number += 1
			ax.set_xlabel('Generation number', fontsize=30)
			ax.set_ylabel('%s counts' % molecule_type, fontsize=30)
			ax.set_title(names[i], fontsize=30)
			ax.tick_params(axis='x', labelsize=30)
			ax.tick_params(axis='y', labelsize=30)
			max_gen = len(start_generation_indices)
			ax.set_xticks(range(max_gen + 1))
			# Set limit from 0 to the total number of generations
			ax.set_xlim(0, max_gen)
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)

			# Generation boundaries are already the integer x ticks above: x is
			# normalized to generation number, so every tick IS a boundary. The
			# old code drew axvlines at `end_generation_times / 60` -- minutes
			# on a generation-number axis, which put them far off the right of
			# the plot (or all at ~0 early on).

		# Remove any empty subplots
		for i in range(num_groups, rows * cols):
			fig.delaxes(axes[i])

		plt.tight_layout()
		exportFigure(plt, plotOutDir,
			plotOutFileName + '_%s_%s_seed%06d' % (
				molecule_type, set_name, seed),
			metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
