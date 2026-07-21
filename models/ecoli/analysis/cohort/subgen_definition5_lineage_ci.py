"""
Per-lineage definition-5 subgenerational classification with confidence intervals.

Definition 5 is the mean number of completed transcripts per generation, per gene.
This analysis estimates it independently for each lineage (one seed = one
independent multi-generation simulation), then aggregates across lineages into a
mean, standard deviation, and 95% confidence interval (from the standard error of
the mean). Each gene is classified by where its CI sits relative to one completed
transcript per generation:

	never_expressed  -- mean == 0 across all lineages
	subgen           -- CI entirely below 1
	possibly_subgen  -- CI includes 1
	not_subgen       -- CI entirely above 1

The classification is computed twice: over all lineages, and over only
"successful" lineages (those that completed every generation and never had a cell
hit the 180-minute simulation-length cap), so we can see whether that filtering
changes the list.

A rate estimated across independent lineages scales well with the number of seeds:
the standard error shrinks as std / sqrt(n_lineages), so more seeds tighten the
confidence interval and let more genes be classified confidently -- unlike a
presence-frequency-at-1 definition, whose "always on" count degrades as the sample
grows.

This module is self-contained and imports no other subgenerational analysis.
"""

import pickle
import os
import csv
import json
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.io.tablereader import TableReader


# All shared constants and the Form B classifier now live in subgen_common, so
# this from-scratch reference implementation and the raw-extraction pipeline
# classify genes identically and cannot drift. The module-level aliases keep the
# local names used throughout this file.
IGNORE_FIRST_N_GENS = sc.IGNORE_FIRST_N_GENS
MAX_DOUBLING_MIN = sc.MAX_DOUBLING_MIN
DOUBLING_AT_MAX_TOL = sc.DOUBLING_AT_MAX_TOL
CI_Z = sc.CI_Z
CONFIDENCE = sc.CONFIDENCE
SYNTH_TABLE = sc.SYNTH_TABLE
SYNTH_COLUMN = sc.SYNTH_COLUMN
CATEGORIES = sc.CATEGORIES
PALETTE = sc.PALETTE
CAT_LABEL = sc.CAT_LABEL
INK, MUTED, GRID, SURF = sc.INK, sc.MUTED, sc.GRID, sc.SURF

# Thin aliases onto the canonical implementations in subgen_common.
_stats = sc.classify_def5_ci
_category_counts = sc.category_counts
_git_info = sc.git_info
_load_sim_metadata = sc.load_sim_metadata


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		analysis_run_time = datetime.now().isoformat(timespec='seconds')
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		n_generation = self.ap.n_generation
		if n_generation <= IGNORE_FIRST_N_GENS:
			print('Skipping analysis - not enough generations run.')
			return

		# --- Gene set: protein-coding mRNA cistrons, in cistron order ---
		mRNA_cistron_ids, monomer_ids, gene_ids = sc.get_mrna_gene_set(sim_data)
		n_genes = len(mRNA_cistron_ids)

		sim_metadata_path, sim_metadata = _load_sim_metadata(variantDir)
		total_init_sims = sim_metadata.get('total_init_sims')

		# Strict successful-lineage classification (the canonical implementation
		# shared with every other subgen analysis).
		success = sc.compute_lineage_success(
			self.ap, n_generation, total_init_sims=total_init_sims)
		seeds = success['seeds']
		all_seed_ids = success['all_seed_ids']
		doubling = success['doubling']
		successful_gens = success['successful_gens']
		n_at_180 = success['n_at_180']
		gens_at_180 = success['gens_at_180']
		completed_all = success['completed_all']
		in_successful = success['in_successful']
		print('Seeds present: %d ; generations: %d ; genes: %d'
			% (len(seeds), n_generation, n_genes))

		# --- Detect the def-5 listener column and build its index map ---
		has_synth, synth_indexes = self._synth_index_map(
			seeds, n_generation, mRNA_cistron_ids)
		if not has_synth:
			print('WARNING: %s/%s not found in this cohort. Definition-5 outputs '
				'and plots are skipped; doubling-time and lineage-filter tables '
				'are still produced. Re-run on simulations built after the '
				'listener was added to populate definition 5.'
				% (SYNTH_TABLE, SYNTH_COLUMN))

		# --- Per-lineage definition 5: mean completed transcripts per cell ---
		lineage_seeds = []
		lambda_rows = []
		if has_synth:
			for s in seeds:
				included = self.ap.get_cells(
					seed=[s],
					generation=np.arange(IGNORE_FIRST_N_GENS, n_generation),
					only_successful=True)
				accum = np.zeros(n_genes)
				n_cells = 0
				for cell_path in included:
					try:
						col = TableReader(os.path.join(cell_path, 'simOut',
							SYNTH_TABLE)).readColumn(SYNTH_COLUMN)
						accum += col.sum(axis=0)[synth_indexes]
						n_cells += 1
					except Exception as e:
						print('  Warning: could not read %s: %s' % (cell_path, e))
				if n_cells > 0:
					lineage_seeds.append(s)
					lambda_rows.append(accum / n_cells)

		# --- Aggregate to per-gene statistics (all and successful) ---
		if has_synth and lineage_seeds:
			lambda_matrix = np.array(lambda_rows)  # (n_lineages_all, n_genes)
			all_idx = np.arange(len(lineage_seeds))
			succ_idx = np.array([
				i for i, s in enumerate(lineage_seeds) if in_successful[s]],
				dtype=int)
			stats_all = _stats(lambda_matrix, all_idx, n_genes)
			stats_succ = _stats(lambda_matrix, succ_idx, n_genes) \
				if succ_idx.size else None
			print('Lineages: %d all, %d successful.'
				% (len(lineage_seeds), succ_idx.size))
		else:
			lambda_matrix = np.zeros((0, n_genes))
			stats_all = stats_succ = None

		# --- Write outputs ---
		prefix = os.path.join(plotOutDir, plotOutFileName)
		if stats_all is not None:
			self._write_pergene(prefix + '_pergene_all.tsv', stats_all,
				gene_ids, mRNA_cistron_ids, monomer_ids, n_genes)
			self._write_lineage_by_gene(prefix + '_lineage_by_gene.tsv',
				lambda_matrix, lineage_seeds, gene_ids, mRNA_cistron_ids, n_genes)
		if stats_succ is not None:
			self._write_pergene(prefix + '_pergene_successful.tsv', stats_succ,
				gene_ids, mRNA_cistron_ids, monomer_ids, n_genes)

		self._write_doubling(prefix + '_doubling_times.tsv',
			doubling, all_seed_ids, n_generation)
		self._write_filter_table(prefix + '_lineage_filter.tsv',
			all_seed_ids, seeds, successful_gens, doubling, n_at_180,
			gens_at_180, in_successful, completed_all, n_generation)

		# --- Plots ---
		if stats_all is not None:
			self._plot_histogram(prefix + '_hist_all.png', stats_all, 'all lineages')
			self._plot_all_vs_successful(
				prefix + '_all_vs_successful.png', stats_all, stats_succ)
		if stats_succ is not None:
			self._plot_histogram(
				prefix + '_hist_successful.png', stats_succ, 'successful lineages')
			self._plot_forest(prefix + '_forest_overlap1.png', stats_succ,
				gene_ids, mRNA_cistron_ids)
			self._plot_ranked(prefix + '_ranked_ci.png', stats_succ)

		# --- Provenance metadata ---
		self._write_metadata(prefix + '_run_metadata.json', variantDir,
			analysis_run_time, sim_metadata_path, sim_metadata, has_synth,
			n_genes, lineage_seeds, in_successful, stats_all, stats_succ)
		print('Done.')

	# ------------------------------------------------------------------ helpers

	def _synth_index_map(self, seeds, n_generation, mRNA_cistron_ids):
		"""Find a cell with the def-5 column and map genes to its cistron order.

		Tries cells until one reads successfully; a single unreadable cell no
		longer disables Definition 5 for the whole cohort.
		"""
		for s in seeds:
			cells = self.ap.get_cells(
				seed=[s],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_generation),
				only_successful=True)
			for cell_path in cells:
				ok, idx = sc.cistron_index_map(
					cell_path, SYNTH_TABLE, mRNA_cistron_ids)
				if ok:
					return True, idx
		return False, None

	def _write_pergene(self, path, st, gene_ids, cistron_ids, monomer_ids,
			n_genes):
		print('Writing %s' % path)
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['gene_id', 'cistron_id', 'protein_id', 'n_lineages',
				'mean', 'std', 'se', 'ci_lower', 'ci_upper', 'category'])
			for i in range(n_genes):
				w.writerow([gene_ids[i], cistron_ids[i], monomer_ids[i][:-3],
					st['n'], st['mean'][i], st['std'][i], st['se'][i],
					st['ci_low'][i], st['ci_high'][i], st['cat'][i]])

	def _write_lineage_by_gene(self, path, lambda_matrix, lineage_seeds,
			gene_ids, cistron_ids, n_genes):
		print('Writing %s' % path)
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['gene_id', 'cistron_id']
				+ ['seed_%06d' % s for s in lineage_seeds])
			for i in range(n_genes):
				w.writerow([gene_ids[i], cistron_ids[i]]
					+ ['%.6g' % lambda_matrix[j, i]
						for j in range(len(lineage_seeds))])

	def _write_doubling(self, path, doubling, all_seed_ids, n_generation):
		print('Writing %s' % path)
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['seed'] + ['gen_%d' % g for g in range(n_generation)])
			for s in all_seed_ids:
				w.writerow(['%06d' % s]
					+ ['%.4g' % x if x >= 0 else '-1' for x in doubling[s]])

	def _write_filter_table(self, path, all_seed_ids, seeds, successful_gens,
			doubling, n_at_180, gens_at_180, in_successful, completed_all,
			n_generation):
		print('Writing %s' % path)
		with open(path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['seed', 'n_gens_ran', 'reached_final_gen',
				'max_doubling_min', 'n_cells_at_180', 'completed_all_gens',
				'passed_180_filter', 'in_successful', 'reason'])
			for s in all_seed_ids:
				gens = successful_gens[s]
				valid = doubling[s][doubling[s] >= 0]
				max_dt = float(valid.max()) if valid.size else -1
				comp = completed_all[s]
				passed180 = n_at_180[s] == 0 and (s in seeds)
				if s not in seeds and not gens:
					reason = 'never ran'
				elif in_successful[s]:
					reason = 'ok'
				else:
					parts = []
					if not comp:
						parts.append('incomplete (%d/%d gens, last %d)' % (
							len(gens), n_generation, max(gens) if gens else -1))
					if n_at_180[s] > 0:
						parts.append('%d cell(s) at 180-min cap (gens %s)' % (
							n_at_180[s],
							','.join(str(g) for g in gens_at_180[s])))
					reason = '; '.join(parts) if parts else 'excluded'
				w.writerow(['%06d' % s, len(gens),
					(n_generation - 1) in gens, '%.4g' % max_dt, n_at_180[s],
					comp, passed180, in_successful[s], reason])

	def _write_metadata(self, path, variantDir, analysis_run_time,
			sim_metadata_path, sim_metadata, has_synth, n_genes, lineage_seeds,
			in_successful, stats_all, stats_succ):
		repo_dir = os.path.dirname(os.path.abspath(__file__))
		meta = {
			'analysis': {
				'script': os.path.basename(__file__),
				'run_time': analysis_run_time,
				'git': _git_info(repo_dir),
				'parameters': {
					'ignore_first_n_gens': IGNORE_FIRST_N_GENS,
					'ci_confidence': CONFIDENCE,
					'ci_z_multiplier': CI_Z,
					'max_doubling_min': MAX_DOUBLING_MIN,
					'definition5_available': has_synth,
					},
				},
			'simulation': {
				'metadata_source': sim_metadata_path,
				'git_hash': sim_metadata.get('git_hash'),
				'git_branch': sim_metadata.get('git_branch'),
				'run_time': sim_metadata.get('time'),
				'description': sim_metadata.get('description'),
				'variant': sim_metadata.get('variant'),
				'total_gens': sim_metadata.get('total_gens'),
				'total_init_sims': sim_metadata.get('total_init_sims'),
				},
			'genes': {'n_genes': n_genes},
			'lineages': {
				'n_all': len(lineage_seeds),
				'n_successful': int(sum(
					1 for s in lineage_seeds if in_successful[s])),
				'seeds_all': lineage_seeds,
				'seeds_successful': [
					s for s in lineage_seeds if in_successful[s]],
				},
			'category_counts': {
				'all': _category_counts(stats_all['cat'])
					if stats_all else None,
				'successful': _category_counts(stats_succ['cat'])
					if stats_succ else None,
				},
			}
		print('Writing %s' % path)
		with open(path, 'w') as f:
			json.dump(meta, f, indent=2)

	# -------------------------------------------------------------------- plots

	def _base_axes(self):
		plt.rcParams.update({
			'font.family': 'DejaVu Sans', 'font.size': 11,
			'axes.edgecolor': MUTED, 'text.color': INK, 'axes.labelcolor': INK,
			'xtick.color': MUTED, 'ytick.color': MUTED})

	def _plot_histogram(self, path, st, label):
		self._base_axes()
		mean, cat = st['mean'], st['cat']
		pos = mean[mean > 0]
		never = int(np.sum(mean == 0))
		counts = _category_counts(cat)
		total = len(mean)
		if pos.size == 0:
			return
		lo, hi = np.floor(np.log10(pos.min())), np.ceil(np.log10(pos.max()))
		edges = np.logspace(lo, hi, int((hi - lo) * 10) + 1)
		hist, _ = np.histogram(pos, bins=edges)
		centers = np.sqrt(edges[:-1] * edges[1:])
		colors = [PALETTE['subgen'] if c < 1 else PALETTE['not_subgen']
			for c in centers]
		ymax = max(hist.max(), never, 1) * 1.12

		fig = plt.figure(figsize=(9.4, 5.4), dpi=160)
		fig.patch.set_facecolor(SURF)
		gs = GridSpec(1, 2, width_ratios=[1, 16], wspace=0.05,
			left=0.09, right=0.975, top=0.86, bottom=0.13)
		axm = fig.add_subplot(gs[1])
		axn = fig.add_subplot(gs[0], sharey=axm)
		for ax in (axm, axn):
			ax.set_facecolor(SURF)
			ax.spines['top'].set_visible(False)
			ax.set_axisbelow(True)

		axn.bar([0], [never], width=0.72, color=PALETTE['never_expressed'],
			zorder=3)
		axn.set_xlim(-0.55, 0.55)
		axn.set_xticks([0])
		axn.set_xticklabels(['0'])
		axn.spines['right'].set_visible(False)
		axn.yaxis.grid(True, color=GRID, lw=0.8)
		axn.set_ylabel('Number of genes')
		axn.set_ylim(0, ymax)
		if never:
			axn.annotate('%d' % never, (0, never), textcoords='offset points',
				xytext=(0, 4), ha='center', fontsize=9.5,
				color=PALETTE['never_expressed'], fontweight='bold')

		axm.bar(edges[:-1], hist, width=np.diff(edges), align='edge',
			color=colors, edgecolor=SURF, linewidth=0.5, zorder=3)
		axm.set_xscale('log')
		axm.set_xlim(edges[0], edges[-1])
		axm.spines['left'].set_visible(False)
		axm.tick_params(axis='y', length=0)
		plt.setp(axm.get_yticklabels(), visible=False)
		axm.yaxis.grid(True, color=GRID, lw=0.8)
		axm.set_xlabel('Mean completed transcripts per generation  (log scale)')
		axm.axvline(1, color=INK, lw=1.4, ls=(0, (5, 3)), zorder=4)
		axm.annotate('subgenerational  <  1 / gen', (1, ymax * 0.52),
			xytext=(-7, 0), textcoords='offset points', ha='right',
			va='center', rotation=90, fontsize=10, color=INK, fontweight='bold')

		fig.text(0.09, 0.945, 'Per-gene mean rate (definition 5) - %s' % label,
			fontsize=15, fontweight='bold', color=INK)
		fig.text(0.09, 0.90, '%d genes  ;  bars colored by threshold; category '
			'counts (CI-based) below' % total, fontsize=9.5, color=MUTED)
		legend = [Patch(fc=PALETTE[c], ec='none',
			label='%s: %d' % (CAT_LABEL[c], counts[c])) for c in CATEGORIES]
		axm.legend(handles=legend, loc='upper right', frameon=False,
			fontsize=9.5, handlelength=1.1, handleheight=1.1, borderpad=0.2)
		fig.savefig(path, facecolor=SURF)
		plt.close(fig)
		print('Writing %s' % path)

	def _plot_forest(self, path, st, gene_ids, cistron_ids):
		self._base_axes()
		idx = np.where(st['cat'] == 'possibly_subgen')[0]
		fig_dpi = 150
		if idx.size == 0:
			fig, ax = plt.subplots(figsize=(7, 2), dpi=fig_dpi)
			fig.patch.set_facecolor(SURF)
			ax.text(0.5, 0.5, 'No genes have a CI that includes 1.',
				ha='center', va='center', color=MUTED)
			ax.axis('off')
			fig.savefig(path, facecolor=SURF)
			plt.close(fig)
			print('Writing %s' % path)
			return
		order = idx[np.argsort(st['mean'][idx])]
		n = order.size
		label_genes = n <= 60
		height = min(42, max(3.5, n * (0.22 if label_genes else 0.1)))
		fig, ax = plt.subplots(figsize=(8.2, height), dpi=fig_dpi)
		fig.patch.set_facecolor(SURF)
		ax.set_facecolor(SURF)
		y = np.arange(n)
		mean = st['mean'][order]
		err_low = mean - st['ci_low'][order]
		err_high = st['ci_high'][order] - mean
		ax.errorbar(mean, y, xerr=[err_low, err_high], fmt='o',
			ms=4 if label_genes else 2.5, color=PALETTE['possibly_subgen'],
			ecolor=PALETTE['possibly_subgen'], elinewidth=1.2, capsize=2,
			alpha=0.9, zorder=3)
		ax.axvline(1, color=INK, lw=1.4, ls=(0, (5, 3)), zorder=2)
		ax.set_ylim(-1, n)
		ax.set_xlabel('Mean completed transcripts per generation  (95% CI)')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.grid(True, color=GRID, lw=0.8)
		ax.set_axisbelow(True)
		if label_genes:
			ax.set_yticks(y)
			ax.set_yticklabels(['%s (%s)' % (gene_ids[i], cistron_ids[i])
				for i in order], fontsize=7.5)
		else:
			ax.set_yticks([])
			ax.set_ylabel('%d genes whose CI includes 1 (sorted by mean)' % n)
		ax.set_title('Boundary genes: 95%% CI includes 1  (%d genes)' % n,
			fontsize=13, fontweight='bold', color=INK, loc='left')
		fig.tight_layout()
		fig.savefig(path, facecolor=SURF)
		plt.close(fig)
		print('Writing %s' % path)

	def _plot_ranked(self, path, st):
		self._base_axes()
		mean = st['mean']
		expressed = np.where(mean > 0)[0]
		order = expressed[np.argsort(mean[expressed])]
		x = np.arange(order.size)
		m = mean[order]
		lo = np.maximum(st['ci_low'][order], m * 1e-3 + 1e-6)
		hi = st['ci_high'][order]
		cats = st['cat'][order]
		colors = np.array([PALETTE[c] for c in cats])

		fig, ax = plt.subplots(figsize=(10, 5.6), dpi=160)
		fig.patch.set_facecolor(SURF)
		ax.set_facecolor(SURF)
		ax.fill_between(x, lo, hi, color=MUTED, alpha=0.18, lw=0, zorder=1)
		ax.scatter(x, m, s=4, c=colors, zorder=3, linewidths=0)
		ax.set_yscale('log')
		ax.axhline(1, color=INK, lw=1.4, ls=(0, (5, 3)), zorder=2)
		ax.set_xlim(0, order.size)
		ax.set_xlabel('Genes ordered by mean rate')
		ax.set_ylabel('Mean completed transcripts / gen  (95% CI band)')
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.yaxis.grid(True, color=GRID, lw=0.8)
		ax.set_axisbelow(True)
		counts = _category_counts(st['cat'])
		legend = [Patch(fc=PALETTE[c], ec='none',
			label='%s: %d' % (CAT_LABEL[c], counts[c]))
			for c in ['subgen', 'possibly_subgen', 'not_subgen']]
		ax.legend(handles=legend, loc='upper left', frameon=False, fontsize=9.5,
			handlelength=1.1, handleheight=1.1)
		ax.set_title('Definition-5 rate with 95% CI, all genes ranked '
			'(successful lineages)', fontsize=13, fontweight='bold',
			color=INK, loc='left')
		fig.tight_layout()
		fig.savefig(path, facecolor=SURF)
		plt.close(fig)
		print('Writing %s' % path)

	def _plot_all_vs_successful(self, path, stats_all, stats_succ):
		self._base_axes()
		fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 5.2), dpi=160,
			gridspec_kw={'width_ratios': [1.15, 1]})
		fig.patch.set_facecolor(SURF)
		for ax in (ax0, ax1):
			ax.set_facecolor(SURF)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		if stats_succ is not None:
			ma, ms = stats_all['mean'], stats_succ['mean']
			both = (ma > 0) & (ms > 0)
			ax0.scatter(ma[both], ms[both], s=5, c=MUTED, alpha=0.4, linewidths=0)
			lim_lo = min(ma[both].min(), ms[both].min())
			lim_hi = max(ma[both].max(), ms[both].max())
			ax0.plot([lim_lo, lim_hi], [lim_lo, lim_hi], color=INK, lw=1,
				ls=(0, (4, 3)))
			ax0.axvline(1, color=PALETTE['not_subgen'], lw=1, alpha=0.6)
			ax0.axhline(1, color=PALETTE['not_subgen'], lw=1, alpha=0.6)
			ax0.set_xscale('log')
			ax0.set_yscale('log')
			ax0.set_xlabel('mean rate - all lineages')
			ax0.set_ylabel('mean rate - successful lineages')
			ax0.set_title('Per-gene mean: all vs successful', fontsize=12,
				fontweight='bold', color=INK, loc='left')
		else:
			ax0.text(0.5, 0.5, 'No successful lineages', ha='center',
				va='center', color=MUTED)
			ax0.axis('off')

		# Category-count bars
		ca = _category_counts(stats_all['cat'])
		cs = _category_counts(stats_succ['cat']) if stats_succ else \
			{c: 0 for c in CATEGORIES}
		xpos = np.arange(len(CATEGORIES))
		width = 0.4
		ax1.bar(xpos - width / 2, [ca[c] for c in CATEGORIES], width,
			color=[PALETTE[c] for c in CATEGORIES], edgecolor=SURF, label='all')
		ax1.bar(xpos + width / 2, [cs[c] for c in CATEGORIES], width,
			color=[PALETTE[c] for c in CATEGORIES], edgecolor=SURF, alpha=0.55,
			hatch='//', label='successful')
		for i, c in enumerate(CATEGORIES):
			ax1.annotate(str(ca[c]), (i - width / 2, ca[c]), ha='center',
				va='bottom', fontsize=8, color=INK)
			ax1.annotate(str(cs[c]), (i + width / 2, cs[c]), ha='center',
				va='bottom', fontsize=8, color=INK)
		ax1.set_xticks(xpos)
		ax1.set_xticklabels(['subgen', 'possibly', 'not', 'never'], fontsize=9)
		ax1.set_ylabel('Number of genes')
		ax1.yaxis.grid(True, color=GRID, lw=0.8)
		ax1.set_axisbelow(True)
		ax1.legend(frameon=False, fontsize=9)
		ax1.set_title('Genes per category', fontsize=12, fontweight='bold',
			color=INK, loc='left')
		fig.tight_layout()
		fig.savefig(path, facecolor=SURF)
		plt.close(fig)
		print('Writing %s' % path)


if __name__ == '__main__':
	Plot().cli()
