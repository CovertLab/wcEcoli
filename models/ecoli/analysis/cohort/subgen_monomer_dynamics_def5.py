"""
Definition-5 rewrite of subgen_monomer_dynamics.py.

Plots mRNA and monomer count dynamics (per generation) for subgenerational
genes, but selects those genes from the canonical Definition-5 def5_CI
classification (subgen iff the 95% CI is below 1 completed transcript/gen across
successful lineages) instead of the original hardcoded curated list.

Gene *selection* comes from the pre-computed raw extraction (run
subgen_raw_extract.py first); the time-series traces themselves are read from
simOut, since dynamics need the full within-generation trace.
"""

import os
import pickle

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, stacked_cell_identification)
from wholecell.io.tablereader import TableReader

# How many subgen genes to plot (highest Def-5 mean first, for clearest traces)
# and how many successful seeds to draw.
N_SUBGEN_TO_PLOT = 12
N_SEEDS_TO_PLOT = 3
COLOR_LINE = 'mediumseagreen'


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		clf = sc.canonical_def5_classification(plotOutDir)
		if clf['n_lineages'] == 0:
			print('No successful lineages found. Skipping.')
			return
		gene_ids, cistron_ids, monomer_ids = sc.load_raw_genes(plotOutDir)
		cat = clf['stats']['cat']
		mean = clf['stats']['mean']

		# Select the N subgen genes with the highest Def-5 mean (clearest traces).
		subgen_idx = np.where(cat == 'subgen')[0]
		if subgen_idx.size == 0:
			print('No subgenerational genes under Definition 5. Skipping.')
			return
		subgen_idx = subgen_idx[np.argsort(mean[subgen_idx])[::-1]]
		chosen = subgen_idx[:N_SUBGEN_TO_PLOT]
		chosen_cistrons = [cistron_ids[i] for i in chosen]
		chosen_monomers = [monomer_ids[i] for i in chosen]
		chosen_names = [gene_ids[i] for i in chosen]
		print('Plotting %d subgen genes: %s'
			% (len(chosen), ', '.join(chosen_names)))

		# Successful seeds to draw (first few, deterministic).
		seeds_to_plot = sorted(clf['successful_seeds'])[:N_SEEDS_TO_PLOT]
		if not seeds_to_plot:
			print('No successful seeds to plot. Skipping.')
			return

		# Map chosen genes into the listener subcolumns (from a real cell).
		probe = self.ap.get_cells(
			seed=[seeds_to_plot[0]],
			generation=np.arange(sc.IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(probe) == 0:
			print('No cells for probe seed. Skipping.')
			return
		rna_reader = TableReader(os.path.join(probe[0], 'simOut', 'RNACounts'))
		mRNA_ids = rna_reader.readAttribute('mRNA_cistron_ids')
		rna_reader.close()
		mRNA_id_to_index = {c: i for i, c in enumerate(mRNA_ids)}
		monomer_reader = TableReader(
			os.path.join(probe[0], 'simOut', 'MonomerCounts'))
		monomer_ids_table = monomer_reader.readAttribute('monomerIds')
		monomer_reader.close()
		monomer_id_to_index = {m: i for i, m in enumerate(monomer_ids_table)}

		mRNA_indices = np.array([mRNA_id_to_index[c] for c in chosen_cistrons])
		monomer_indices = np.array([
			monomer_id_to_index[m] for m in chosen_monomers])

		for seed in seeds_to_plot:
			cell_paths_per_seed = self.ap.get_cells(
				generation=np.arange(
					sc.IGNORE_FIRST_N_GENS, self.ap.n_generation),
				seed=[seed], only_successful=True)
			if len(cell_paths_per_seed) == 0:
				continue

			monomer_counts = read_stacked_columns(
				cell_paths_per_seed, 'MonomerCounts', 'monomerCounts'
				)[:, monomer_indices]
			cistron_counts = read_stacked_columns(
				cell_paths_per_seed, 'RNACounts', 'mRNA_cistron_counts'
				)[:, mRNA_indices]
			(time, _, end_generation_times, start_generation_indices,
				end_generation_indices) = self._extract_doubling_times(
				cell_paths_per_seed)

			self._plot_counts_dynamics_gen(
				plotOutDir, plotOutFileName, metadata, monomer_counts,
				chosen_names, 'monomer', time, start_generation_indices,
				end_generation_indices, end_generation_times, seed)
			self._plot_counts_dynamics_gen(
				plotOutDir, plotOutFileName, metadata, cistron_counts,
				chosen_names, 'mRNA', time, start_generation_indices,
				end_generation_indices, end_generation_times, seed)

	# ------------------------------------------------------------------ helpers

	def _extract_doubling_times(self, cell_paths):
		time = read_stacked_columns(cell_paths, 'Main', 'time').squeeze()
		doubling_times = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: (x[-1] - x[0])).squeeze().astype(int)
		end_generation_times = np.cumsum(doubling_times) + time[0]
		# Per-generation row boundaries from actual row counts. read_stacked_columns
		# stacks one cell (= one generation) per block, so labeling each row by its
		# source cell segments generations exactly, with no assumption that a row is
		# 1 s and no boundary-duplicate drift.
		cell_ids = stacked_cell_identification(
			cell_paths, 'Main', 'time').squeeze().astype(int)
		n_cells = len(cell_paths)
		start_generation_indices = np.searchsorted(
			cell_ids, np.arange(n_cells), side='left')
		# Inclusive last row of each generation; callers slice time[start:end + 1].
		end_generation_indices = np.append(
			start_generation_indices[1:], len(cell_ids)) - 1
		return (time.astype(int), doubling_times, end_generation_times,
			start_generation_indices, end_generation_indices)

	def _plot_counts_dynamics_gen(self, plotOutDir, plotOutFileName, metadata,
			counts, names, molecule_type, time, start_generation_indices,
			end_generation_indices, end_generation_times, seed):
		num_groups = len(names)
		cols = min(num_groups, 4)
		rows = (num_groups + cols - 1) // cols
		fig_width = cols * 10
		fig_height = fig_width / 2
		fig, axes = plt.subplots(
			nrows=rows, ncols=cols, figsize=(fig_width, fig_height), sharex=False)
		axes = np.atleast_1d(axes).flatten()

		for i, ax in enumerate(axes):
			if i >= num_groups:
				continue
			interest_counts = counts[:, i]
			gen_number = 0
			for start_idx, end_idx in zip(
					start_generation_indices, end_generation_indices):
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
			ax.set_xlim(0, max_gen)
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)

		for i in range(num_groups, rows * cols):
			fig.delaxes(axes[i])

		plt.tight_layout()
		exportFigure(plt, plotOutDir,
			plotOutFileName + '_%s_dynamics_gen_%s' % (molecule_type, seed),
			metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
