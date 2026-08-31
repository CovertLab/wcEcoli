"""
mRNA and monomer count traces, for three gene sets, on two x axes.

The same per-timestep traces are rendered under two x transforms:

  generation  x(t) = gen_number + (t - t_start) / (t_end - t_start), so every
              generation spans unit width. Shows per-generation behaviour: did
              the gene fire in this cell cycle at all, and did the burst carry
              over into the daughters.
  time        x(t) = (t - t_first) / 60, in minutes, rebased so the first
              plotted generation starts at 0. Shows kinetics: accumulation
              slope, degradation, and how long a protein persists after its
              mRNA stops.

Gene sets, all plotted in one run, with the set name in every output filename:

  anaerobic   the fixed anaerobic-respiration panel, one subunit per complex
              (sc.curated_panel('anaerobic')), over SEEDS_ANAEROBIC.
  curated10   the 10-gene panel shared with subgen_peak_counts.py and
              subgen_protein_distribution.py (sc.curated_panel('curated10')),
              over seeds chosen by gene coverage (see _coverage_ranked_seeds).
  subgen_top  the subgenerational genes with the highest
              Definition-5 mean, from the extraction cache, over the first
              N_SEEDS_TO_PLOT strict-successful seeds.

Writes, per gene set per plotted seed:
  <plotOutFileName>_{monomer,mRNA}_<set>_seed<SSSSSS>.pdf       generation axis
  <plotOutFileName>_{monomer,mRNA}_<set>_time_seed<SSSSSS>.pdf  minutes axis
  <plotOutFileName>_traces_<set>_seed<SSSSSS>.tsv               plotted traces
  <plotOutFileName>_run_metadata.json                           provenance

The generation-axis filenames are unchanged from before the minutes axis was
added, so existing outputs regenerate under exactly the same names.

Definition 5 always means the CI form (sc.classify_def5_ci). The `subgen_top`
set and curated10's coverage-based seed ranking both read subgen_extract.py's
cache, and each degrades on its own if it is missing: subgen_top is skipped with
a message, curated10 falls back to the first N strict-successful seeds.

Time-series traces always come from simOut, since dynamics need the full
within-generation trace; only gene SELECTION and seed RANKING are cached.
"""

import csv
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
SEEDS_ANAEROBIC = np.arange(20, 25)
N_SUBGEN_TO_PLOT = 12
N_SEEDS_TO_PLOT = 3
COLOR_LINE = 'mediumseagreen'  # 'skyblue'
PANEL = 'anaerobic'
PANEL_CURATED = 'curated10'
SEEDS_CURATED10 = None
N_SEEDS_CURATED10 = 3
CURATED10_COLS = 5
SET_COLS = {PANEL_CURATED: CURATED10_COLS}
X_AXES = ('generation', 'time')
SHADE_ZERO = True
ZERO_SHADE_ALPHA = 0.15
WRITE_TRACE_TSV = True


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

		gene_sets = [
			self._anaerobic_gene_set(sim_data, successful),
			self._curated_gene_set(sim_data, successful, plotOutDir),
			]
		subgen_set = self._subgen_gene_set(plotOutDir, successful)
		if subgen_set is not None:
			gene_sets.append(subgen_set)

		plotted = {}
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
			plotted[set_name] = {
				'n_genes': len(names),
				'genes': list(names),
				'seeds': [int(s) for s in seeds],
				}

		sim_metadata_path, sim_metadata = sc.load_sim_metadata(variantDir)
		sc.write_run_metadata(
			os.path.join(plotOutDir, plotOutFileName + '_run_metadata.json'),
			'subgen_monomer_dynamics.py',
			{
				'ignore_first_n_gens': IGNORE_FIRST_N_GENS,
				'x_axes': list(X_AXES),
				'shade_zero': SHADE_ZERO,
				'n_subgen_to_plot': N_SUBGEN_TO_PLOT,
				'n_seeds_curated10': N_SEEDS_CURATED10,
				'panels': [PANEL, PANEL_CURATED],
				},
			sim_metadata_path=sim_metadata_path,
			sim_metadata=sim_metadata,
			extra={'gene_sets': plotted,
				'n_successful_seeds': len(successful)},
			)

	# gene sets

	def _anaerobic_gene_set(self, sim_data, successful):
		"""The fixed anaerobic panel, ordered by monomer id."""
		monomers_of_interest, name_dict = sc.curated_panel(PANEL)
		cistron_ids = self._cistrons_for(sim_data, monomers_of_interest)
		names = [name_dict[m] for m in monomers_of_interest]
		seeds = [int(s) for s in SEEDS_ANAEROBIC if s in successful]
		if len(seeds) < len(SEEDS_ANAEROBIC):
			print('[%s] Plotting %d of %d requested seeds; the rest are not '
				'strict-successful seeds. This is the only reason the figure '
				'count differs between cohorts.'
				% (PANEL, len(seeds), len(SEEDS_ANAEROBIC)))
		return (PANEL, names, cistron_ids, list(monomers_of_interest), seeds)

	def _curated_gene_set(self, sim_data, successful, plotOutDir):
		monomers_of_interest, name_dict = sc.curated_panel(PANEL_CURATED)
		cistron_ids = self._cistrons_for(sim_data, monomers_of_interest)
		names = [name_dict[m] for m in monomers_of_interest]
		if SEEDS_CURATED10 is not None:
			seeds = [int(s) for s in SEEDS_CURATED10 if s in successful]
			if len(seeds) < len(SEEDS_CURATED10):
				print('[%s] Plotting %d of %d pinned seeds; the rest are not '
					'strict-successful.'
					% (PANEL_CURATED, len(seeds), len(SEEDS_CURATED10)))
		else:
			seeds = self._coverage_ranked_seeds(
				plotOutDir, monomers_of_interest, successful,
				N_SEEDS_CURATED10)
		return (PANEL_CURATED, names, cistron_ids,
			list(monomers_of_interest), seeds)

	def _cistrons_for(self, sim_data, monomer_ids):
		protein_id_to_cistron_id = {
			protein['id']: protein['cistron_id']
			for protein in sim_data.process.translation.monomer_data
			}
		return [protein_id_to_cistron_id[m] for m in monomer_ids]

	def _coverage_ranked_seeds(self, plotOutDir, monomer_ids, successful, n):
		"""Seeds ranked by how many of these genes are ever present.

		`sorted(successful)[:n]` is the right rule for subgen_top, which picks
		the HIGHEST-expressed subgen genes, but wrong for a fixed panel holding
		rarely expressed ones: on sim set 1 it returns seeds 1, 2, 3, and lacZ
		is zero in every generation of seeds 2 and 3, so two of three figures
		would carry a dead flat lacZ panel. Rank instead by how many of the
		panel's genes are ever nonzero across the generations actually plotted,
		breaking ties on the total number of gene-generations with signal.

		Ranking on the plotted window matters: a seed can otherwise rank highly
		on bursts that fall outside it. Reads the extraction cache (one small
		TSV); falls back to the default rule with a message when it is absent.
		"""
		fallback = sorted(int(s) for s in successful)[:n]
		try:
			seeds_arr, generations, is_succ, gene_ids, matrix = sc.load_raw_max(
				plotOutDir, 'protein')
			raw_gene_ids, _, raw_monomer_ids = sc.load_raw_genes(plotOutDir)
		except FileNotFoundError:
			print('[%s] No extraction in %s, so seeds fall back to the first %d '
				'strict-successful seeds (%s). Run subgen_extract.py to rank '
				'seeds by panel coverage instead.'
				% (PANEL_CURATED, plotOutDir, n,
					', '.join(str(s) for s in fallback)))
			return fallback

		# The extraction matrix is keyed by gene_id; the panel is monomer ids.
		monomer_to_gene = dict(zip(raw_monomer_ids, raw_gene_ids))
		gene_index = {g: i for i, g in enumerate(gene_ids)}
		cols, missing = [], []
		for m in monomer_ids:
			gene_id = monomer_to_gene.get(m)
			if gene_id is None or gene_id not in gene_index:
				missing.append(m)
				continue
			cols.append(gene_index[gene_id])
		if missing:
			print('[%s] %d panel monomers are absent from the extraction and do '
				'not count toward coverage: %s'
				% (PANEL_CURATED, len(missing), ', '.join(missing)))
		if not cols:
			return fallback

		keep = is_succ & (generations >= IGNORE_FIRST_N_GENS)
		if not keep.any():
			return fallback
		present = matrix[np.ix_(keep, np.array(cols))] > 0
		kept_seeds = seeds_arr[keep]

		eligible = {int(s) for s in successful}
		ranked = []
		for seed in sorted({int(s) for s in kept_seeds} & eligible):
			rows = present[kept_seeds == seed]
			if not rows.size:
				continue
			per_gene = rows.sum(axis=0)
			ranked.append(
				(int((per_gene > 0).sum()), int(per_gene.sum()), seed))
		if not ranked:
			return fallback
		# Most genes with signal first, then most gene-generations, then seed.
		ranked.sort(key=lambda t: (-t[0], -t[1], t[2]))
		chosen = [seed for _, _, seed in ranked[:n]]
		print('[%s] Seeds ranked by panel coverage over generations >= %d: '
			'chose %s (genes with signal / gene-generations -- %s).'
			% (PANEL_CURATED, IGNORE_FIRST_N_GENS,
				', '.join(str(s) for s in chosen),
				'; '.join('seed %d: %d of %d / %d' % (seed, n_genes, len(cols),
					total) for n_genes, total, seed in ranked[:n])))
		return chosen

	def _subgen_gene_set(self, plotOutDir, successful):
		"""The top-N def5_CI subgen genes, or None if the extraction is absent.
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
		cols = SET_COLS.get(set_name)

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
				for x_axis in X_AXES:
					self._plot_counts_dynamics(
						plotOutDir, plotOutFileName, metadata, set_name, counts,
						names, molecule_type, time, start_generation_indices,
						end_generation_indices, seed, x_axis, cols)

			if WRITE_TRACE_TSV:
				self._write_trace_tsv(
					plotOutDir, plotOutFileName, set_name, seed, names,
					monomer_counts, time, cell_paths_per_seed,
					start_generation_indices, end_generation_indices)

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
		return time, start_generation_indices, end_generation_indices

	def _shade_zero_runs(self, ax, x, counts):
		"""Shade contiguous stretches where the count is zero."""
		zero = np.asarray(counts) == 0
		if not zero.any():
			return
		# Run boundaries from the 0/1 mask, zero-padded so runs touching either
		# end are closed.
		edges = np.diff(np.concatenate(
			([0], zero.astype(np.int8), [0])))
		for start, end in zip(np.where(edges == 1)[0],
				np.where(edges == -1)[0] - 1):
			ax.axvspan(x[start], x[end], color=sc.MUTED,
				alpha=ZERO_SHADE_ALPHA, linewidth=0, zorder=0)

	def _plot_counts_dynamics(self, plotOutDir, plotOutFileName, metadata,
			set_name, counts, names, molecule_type, time,
			start_generation_indices, end_generation_indices, seed, x_axis,
			cols=None):
		num_groups = len(names)
		if cols is None:
			cols = min(num_groups, 4)  # Number of columns for the grid
		cols = min(cols, num_groups)
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
				if x_axis == 'generation':
					duration = gen_time[-1] - gen_time[0]
					x = gen_number + (gen_time - gen_time[0]) / duration
				else:
					x = (gen_time - time[0]) / 60.
				ax.plot(x, gen_counts, color=COLOR_LINE, linewidth=6)
				if SHADE_ZERO:
					self._shade_zero_runs(ax, x, gen_counts)
				gen_number += 1

			ax.set_ylabel('%s counts' % molecule_type, fontsize=30)
			ax.set_title(names[i], fontsize=30)
			ax.tick_params(axis='x', labelsize=30)
			ax.tick_params(axis='y', labelsize=30)
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)

			if x_axis == 'generation':
				ax.set_xlabel('Generation number', fontsize=30)
				max_gen = len(start_generation_indices)
				ax.set_xticks(range(max_gen + 1))
				# Set limit from 0 to the total number of generations
				ax.set_xlim(0, max_gen)
			else:
				ax.set_xlabel('Time (min)', fontsize=30)
				# Minutes on a MINUTES axis: this is the correct form of the
				# line removed from the generation axis above, not a regression.
				for end_idx in end_generation_indices[:-1]:
					ax.axvline((time[end_idx] - time[0]) / 60.,
						color=sc.GRID, linewidth=2, zorder=0)

		# Remove any empty subplots
		for i in range(num_groups, rows * cols):
			fig.delaxes(axes[i])

		plt.tight_layout()
		# The generation axis keeps its historical filename (no infix), so
		# existing outputs regenerate under the same names.
		infix = '' if x_axis == 'generation' else '_time'
		exportFigure(plt, plotOutDir,
			plotOutFileName + '_%s_%s%s_seed%06d' % (
				molecule_type, set_name, infix, seed),
			metadata)
		plt.close('all')

	def _write_trace_tsv(self, plotOutDir, plotOutFileName, set_name, seed,
			names, monomer_counts, time, cell_paths,
			start_generation_indices, end_generation_indices):
		"""Dump the plotted monomer traces, so both figures can be redrawn
		locally without simOut.

		One row per timestep: the absolute generation number, absolute time in
		seconds, minutes rebased to the first plotted generation, the
		within-generation fraction (the generation axis's x, minus the integer
		part), then one column per gene.
		"""
		path = os.path.join(plotOutDir, plotOutFileName
			+ '_traces_%s_seed%06d.tsv' % (set_name, seed))
		n_rows = len(time)
		generation = np.full(n_rows, -1, dtype=int)
		gen_fraction = np.zeros(n_rows)
		for cell_path, start_idx, end_idx in zip(
				cell_paths, start_generation_indices, end_generation_indices):
			generation[start_idx:end_idx + 1] = sc.parse_cell_id(cell_path)[1]
			span = time[end_idx] - time[start_idx]
			if span > 0:
				gen_fraction[start_idx:end_idx + 1] = (
					time[start_idx:end_idx + 1] - time[start_idx]) / span

		t0 = time[0]
		with open(path, 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(
				['generation', 'time_s', 'time_min', 'gen_fraction']
				+ list(names))
			for r in range(n_rows):
				writer.writerow([
					generation[r],
					'%.3f' % time[r],
					'%.5f' % ((time[r] - t0) / 60.),
					'%.6f' % gen_fraction[r],
					] + [int(v) for v in monomer_counts[r]])
		print('  wrote %s' % path)


if __name__ == '__main__':
	Plot().cli()
