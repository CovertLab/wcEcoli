"""
Quantify whether the SET of subgenerational genes expressed per cell is unique,
while the NUMBER expressed is nearly constant.

Inputs
------
  binary_matrix.tsv : first column `cistron_id`, remaining columns one per cell,
      values 0/1 (1 = the cell had >=1 completed transcript of that gene, i.e. the
      same definition-5 "expressed" criterion). Header row carries the cell ids.
  pergene_successful.tsv : the definition-5 per-gene table; genes with
      category == 'subgen' (matched on cistron_id) define the gene set.

What it does
------------
  A. NUMBER per cell    -- distribution of how many subgen genes each cell
     expresses, compared to the Poisson-binomial expectation (are counts as tight
     as / tighter than independent expression predicts?).
  B. UNIQUENESS of set  -- raw distinct-set count (with the combinatorial context
     that ~100% unique is expected by chance), plus the pairwise Jaccard overlap
     distribution versus two nulls: an analytic independence null and a
     fixed-margins (curveball) null that holds both gene frequencies and per-cell
     counts fixed.
  C. STRUCTURE          -- PCA of the cell x gene matrix (discrete cell types vs a
     structureless stochastic continuum).
  D. POPULATION COVERAGE-- accumulation curve: distinct subgen genes expressed by
     n cells, showing the population collectively covers the pool (bet-hedging)
     while each cell samples a small slice.

Usage
-----
  python subgen_set_uniqueness.py <binary_matrix.tsv> <pergene_successful.tsv> \
      <output_dir> [--category subgen] [--n-null 100]

Self-contained: numpy + matplotlib only.
"""

import os
import csv
import sys
import json
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


INK = '#1b2530'
MUTED = '#5b6672'
GRID = '#e4e8ec'
SURF = '#fcfcfb'
ACCENT = '#1667B8'      # observed
NULLCLR = '#C7362F'     # null / reference
TEAL = '#109C9C'


# ------------------------------------------------------------------ loading

def _to01(values):
	return [1 if v not in ('', '0', '0.0', 'False') else 0 for v in values]


def load_binary_matrix(path):
	"""Load a binary expression matrix in either supported layout, returning
	(row_ids, cell_ids, M, seeds) with M an int8 array [n_genes, n_cells].

	- genes-as-rows: first column is the gene id (e.g. `cistron_id`), remaining
	  columns are one per cell.
	- cells-as-rows (Adrian's `_expressed.tsv`): columns are `seed`, `generation`,
	  then one per gene (header = gene id); each row is one cell.

	`seeds` is an int array of the lineage seed for each cell (parsed from the
	`seed` column or from `seed<N>_gen<N>` labels), or None if unavailable.
	"""
	with open(path) as f:
		reader = csv.reader(f, delimiter='\t')
		header = next(reader)
		if header[0] == 'seed':                       # cells-as-rows
			row_ids = header[2:]
			cols, seeds, cell_ids = [], [], []
			for rec in reader:
				if not rec:
					continue
				seeds.append(int(rec[0]))
				cell_ids.append('seed%06d_gen%03d' % (int(rec[0]), int(rec[1])))
				cols.append(_to01(rec[2:]))
			M = np.array(cols, dtype=np.int8).T       # genes x cells
			return row_ids, cell_ids, M, np.array(seeds)
		cell_ids = header[1:]                          # genes-as-rows
		row_ids, rows = [], []
		for rec in reader:
			if not rec:
				continue
			row_ids.append(rec[0])
			rows.append(_to01(rec[1:]))
	M = np.array(rows, dtype=np.int8)
	return row_ids, cell_ids, M, _seeds_from_labels(cell_ids)


def _seeds_from_labels(cell_ids):
	"""Parse lineage seed from `seed<N>_gen<N>` cell labels, else None."""
	seeds = []
	for c in cell_ids:
		if c.startswith('seed') and '_gen' in c:
			try:
				seeds.append(int(c[4:c.index('_gen')]))
				continue
			except ValueError:
				pass
		return None
	return np.array(seeds)


def load_subgen_ids(path, category='subgen'):
	"""Return (gene_id_set, cistron_id_set) for rows with the given category, so
	the matrix can be matched on whichever id it is keyed by."""
	gene_ids, cistron_ids = set(), set()
	with open(path) as f:
		for r in csv.DictReader(f, delimiter='\t'):
			if r['category'] == category:
				gene_ids.add(r['gene_id'])
				cistron_ids.add(r['cistron_id'])
	return gene_ids, cistron_ids


# ------------------------------------------------------------------ analysis

def subset_to_subgen(M, row_ids, subgen_gene_ids, subgen_cistron_ids):
	"""Keep only the rows (genes) in the subgen set, auto-matching whichever id
	(gene_id or cistron_id) the matrix rows are keyed by."""
	n_gene = sum(1 for r in row_ids if r in subgen_gene_ids)
	n_cistron = sum(1 for r in row_ids if r in subgen_cistron_ids)
	target = subgen_gene_ids if n_gene >= n_cistron else subgen_cistron_ids
	index = {rid: i for i, rid in enumerate(row_ids)}
	kept = [rid for rid in row_ids if rid in target]
	missing = sorted(target - set(row_ids))
	return M[[index[rid] for rid in kept], :], kept, missing


def pairwise_jaccard_mean(B):
	"""Mean pairwise Jaccard over cell pairs (upper triangle), plus the matrix
	diagonal-excluded intersection counts for reuse."""
	Bi = B.astype(np.int32)
	inter = Bi.T @ Bi                      # cells x cells, shared-gene counts
	k = np.diag(inter).astype(np.float64)  # per-cell set sizes
	union = k[:, None] + k[None, :] - inter
	iu = np.triu_indices(inter.shape[0], k=1)
	u = union[iu]
	nz = u > 0
	jac = np.zeros(u.shape)
	jac[nz] = inter[iu][nz] / u[nz]
	return float(jac.mean()), jac, k


def analyze(B):
	n_genes, n_cells = B.shape
	p = B.mean(axis=1)                     # per-gene expression frequency
	k = B.sum(axis=0).astype(np.float64)   # per-cell subgen count

	# A. count uniformity vs Poisson-binomial null
	mu = float(p.sum())
	var_null = float((p * (1 - p)).sum())
	count_stats = {
		'mean': float(k.mean()), 'sd': float(k.std(ddof=1)),
		'cv': float(k.std(ddof=1) / k.mean()) if k.mean() else float('nan'),
		'min': int(k.min()), 'max': int(k.max()),
		'null_mean': mu, 'null_sd': var_null ** 0.5,
		'dispersion_ratio': float(k.var() / var_null) if var_null else float('nan'),
		}

	# B. raw uniqueness (with combinatorial context)
	cols = [tuple(np.nonzero(B[:, c])[0]) for c in range(n_cells)]
	from collections import Counter
	set_counts = Counter(cols)
	n_distinct = len(set_counts)
	largest_group = max(set_counts.values())
	n_unique_cells = sum(1 for v in set_counts.values() if v == 1)

	# B. pairwise overlap observed vs independence null
	obs_jac, jac_vec, _ = pairwise_jaccard_mean(B)
	exp_shared = float((p ** 2).sum())
	exp_size = mu
	exp_jac = exp_shared / (2 * exp_size - exp_shared) if exp_size else float('nan')

	# D. accumulation curve (analytic expectation)
	grid = np.unique(np.clip(
		np.round(np.geomspace(1, n_cells, 40)).astype(int), 1, n_cells))
	cover = np.array([(1 - (1 - p) ** int(n)).sum() for n in grid])

	return {
		'n_genes': n_genes, 'n_cells': n_cells,
		'gene_freq': p, 'per_cell_count': k, 'count_stats': count_stats,
		'n_distinct_sets': n_distinct, 'largest_identical_group': largest_group,
		'n_cells_with_unique_set': n_unique_cells,
		'obs_mean_jaccard': obs_jac, 'null_mean_jaccard_independence': exp_jac,
		'jaccard_vec': jac_vec,
		'accum_grid': grid, 'accum_cover': cover,
		'total_covered': int((p > 0).sum()),
		}


# --------------------------------------------------- fixed-margins (curveball)

def curveball(B, n_iter, rng):
	"""Randomize a binary matrix preserving both row and column sums."""
	rows = [set(np.nonzero(r)[0].tolist()) for r in B]
	R = len(rows)
	for _ in range(n_iter):
		i, j = rng.integers(0, R), rng.integers(0, R)
		if i == j:
			continue
		A, Bs = rows[i], rows[j]
		inter = A & Bs
		excl = list((A - Bs) | (Bs - A))
		if not excl:
			continue
		nA = len(A) - len(inter)
		rng.shuffle(excl)
		rows[i] = inter | set(excl[:nA])
		rows[j] = inter | set(excl[nA:])
	out = np.zeros_like(B)
	for r, s in enumerate(rows):
		if s:
			out[r, list(s)] = 1
	return out


def curveball_null(B, n_null, rng):
	"""Null distribution of mean pairwise Jaccard under fixed margins."""
	n_ones = int(B.sum())
	n_iter = max(2000, 5 * n_ones)
	vals = []
	cur = B.copy()
	for _ in range(n_null):
		cur = curveball(cur, n_iter, rng)
		vals.append(pairwise_jaccard_mean(cur)[0])
	return np.array(vals)


def pca_cells(B, n_components=10):
	"""PCA on the cell x gene matrix; return explained-variance ratio + top-2 scores."""
	X = B.T.astype(np.float64)
	X = X - X.mean(axis=0, keepdims=True)
	U, S, Vt = np.linalg.svd(X, full_matrices=False)
	var = (S ** 2)
	evr = var / var.sum()
	scores = U[:, :2] * S[:2]
	return evr[:n_components], scores


# ------------------------------------------------------------------- plots

def _ax():
	plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
		'axes.edgecolor': MUTED, 'text.color': INK, 'axes.labelcolor': INK,
		'xtick.color': MUTED, 'ytick.color': MUTED})


def plot_count_hist(res, path):
	_ax()
	k = res['per_cell_count']
	cs = res['count_stats']
	fig, ax = plt.subplots(figsize=(8, 4.8), dpi=160)
	fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
	bins = np.arange(k.min(), k.max() + 2) - 0.5
	ax.hist(k, bins=bins, color=ACCENT, edgecolor=SURF, alpha=0.9)
	xs = np.linspace(k.min(), k.max(), 200)
	from math import pi
	pdf = (np.exp(-(xs - cs['null_mean']) ** 2 / (2 * cs['null_sd'] ** 2))
		/ (cs['null_sd'] * (2 * pi) ** 0.5)) * len(k) * (bins[1] - bins[0])
	ax.plot(xs, pdf, color=NULLCLR, lw=2,
		label='independence null (Poisson-binomial)')
	ax.axvline(cs['mean'], color=INK, lw=1.2, ls=(0, (5, 3)))
	ax.set_xlabel('subgen genes expressed per cell')
	ax.set_ylabel('number of cells')
	ax.set_title('How many subgen genes each cell expresses\n'
		'mean %.1f, CV %.2f  (null CV %.2f, dispersion %.2f)'
		% (cs['mean'], cs['cv'], cs['null_sd'] / cs['null_mean'],
			cs['dispersion_ratio']), fontsize=12, color=INK, loc='left')
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
	ax.legend(frameon=False, fontsize=9)
	fig.tight_layout(); fig.savefig(path, facecolor=SURF); plt.close(fig)
	print('wrote', path)


def plot_jaccard(res, path, null_vals=None):
	_ax()
	jac = res['jaccard_vec']
	fig, ax = plt.subplots(figsize=(8, 4.8), dpi=160)
	fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
	ax.hist(jac, bins=60, color=ACCENT, edgecolor=SURF, alpha=0.85,
		label='observed cell pairs')
	ax.axvline(res['obs_mean_jaccard'], color=INK, lw=1.6, ls=(0, (5, 3)),
		label='observed mean %.3f' % res['obs_mean_jaccard'])
	ax.axvline(res['null_mean_jaccard_independence'], color=NULLCLR, lw=1.6,
		label='independence null %.3f' % res['null_mean_jaccard_independence'])
	if null_vals is not None and len(null_vals):
		ax.axvspan(null_vals.min(), null_vals.max(), color=TEAL, alpha=0.15,
			label='fixed-margins null range')
	ax.set_xlabel('pairwise Jaccard overlap of subgen sets')
	ax.set_ylabel('number of cell pairs')
	ax.set_title('How much two cells share their subgen set', fontsize=12,
		color=INK, loc='left')
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
	ax.legend(frameon=False, fontsize=9)
	fig.tight_layout(); fig.savefig(path, facecolor=SURF); plt.close(fig)
	print('wrote', path)


def plot_gene_freq(res, path):
	_ax()
	p = res['gene_freq']
	fig, ax = plt.subplots(figsize=(8, 4.8), dpi=160)
	fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
	ax.hist(p, bins=50, color=TEAL, edgecolor=SURF, alpha=0.9)
	ax.set_xlabel('fraction of cells expressing the gene')
	ax.set_ylabel('number of subgen genes')
	ax.set_title('The sampling pool: how often each subgen gene fires\n'
		'%d genes; %d expressed in >=1 cell' % (len(p), int((p > 0).sum())),
		fontsize=12, color=INK, loc='left')
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
	fig.tight_layout(); fig.savefig(path, facecolor=SURF); plt.close(fig)
	print('wrote', path)


def plot_accumulation(res, path):
	_ax()
	g, cov = res['accum_grid'], res['accum_cover']
	fig, ax = plt.subplots(figsize=(8, 4.8), dpi=160)
	fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
	ax.plot(g, cov, color=ACCENT, lw=2.4)
	ax.axhline(res['n_genes'], color=MUTED, lw=1, ls=(0, (4, 3)),
		label='all %d subgen genes' % res['n_genes'])
	ax.set_xscale('log')
	ax.set_xlabel('number of cells sampled (log)')
	ax.set_ylabel('distinct subgen genes expressed')
	ax.set_title('Population coverage of the subgen pool', fontsize=12,
		color=INK, loc='left')
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
	ax.legend(frameon=False, fontsize=9)
	fig.tight_layout(); fig.savefig(path, facecolor=SURF); plt.close(fig)
	print('wrote', path)


def plot_pca(scores, evr, k, path):
	_ax()
	fig, ax = plt.subplots(figsize=(7.4, 6), dpi=160)
	fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
	sc = ax.scatter(scores[:, 0], scores[:, 1], c=k, s=10, cmap='viridis',
		alpha=0.8, linewidths=0)
	fig.colorbar(sc, ax=ax, label='subgen genes per cell')
	ax.set_xlabel('PC1 (%.1f%%)' % (100 * evr[0]))
	ax.set_ylabel('PC2 (%.1f%%)' % (100 * evr[1]))
	ax.set_title('Cells in subgen-expression space\n'
		'(discrete clusters = programs; a blob = stochastic continuum)',
		fontsize=12, color=INK, loc='left')
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
	fig.tight_layout(); fig.savefig(path, facecolor=SURF); plt.close(fig)
	print('wrote', path)


def within_between_lineage(B, seeds):
	"""Mean pairwise Jaccard for cell pairs in the SAME lineage vs DIFFERENT
	lineages. Higher within-than-between overlap = inherited (lineage) memory."""
	Bi = B.astype(np.int32)
	inter = Bi.T @ Bi
	k = np.diag(inter).astype(np.float64)
	union = k[:, None] + k[None, :] - inter
	iu = np.triu_indices(inter.shape[0], k=1)
	u = union[iu]
	jac = np.zeros(u.shape)
	nz = u > 0
	jac[nz] = inter[iu][nz] / u[nz]
	same = seeds[iu[0]] == seeds[iu[1]]
	return {
		'within_mean': float(jac[same].mean()) if same.any() else float('nan'),
		'between_mean': float(jac[~same].mean()) if (~same).any() else float('nan'),
		'within_vals': jac[same], 'between_vals': jac[~same]}


def plot_within_between(wb, path):
	_ax()
	fig, ax = plt.subplots(figsize=(8, 4.8), dpi=160)
	fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
	bins = np.linspace(0, max(wb['within_vals'].max(), wb['between_vals'].max()), 60)
	ax.hist(wb['between_vals'], bins=bins, color=MUTED, alpha=0.55,
		density=True, label='different lineages')
	ax.hist(wb['within_vals'], bins=bins, color=ACCENT, alpha=0.7,
		density=True, label='same lineage')
	ax.axvline(wb['between_mean'], color=MUTED, lw=1.6, ls=(0, (5, 3)),
		label='between mean %.3f' % wb['between_mean'])
	ax.axvline(wb['within_mean'], color=ACCENT, lw=1.6, ls=(0, (5, 3)),
		label='within mean %.3f' % wb['within_mean'])
	ax.set_xlabel('pairwise Jaccard overlap of subgen sets')
	ax.set_ylabel('density of cell pairs')
	ax.set_title('Does inheritance leave a mark?\n'
		'sister/descendant cells (same lineage) vs unrelated cells',
		fontsize=12, color=INK, loc='left')
	ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
	ax.legend(frameon=False, fontsize=9)
	fig.tight_layout(); fig.savefig(path, facecolor=SURF); plt.close(fig)
	print('wrote', path)


# -------------------------------------------------------------------- main

def run(matrix_path, pergene_path, out_dir, category='subgen', n_null=100,
		seed=0):
	os.makedirs(out_dir, exist_ok=True)
	row_ids, cell_ids, M, seeds = load_binary_matrix(matrix_path)
	subgen_gene_ids, subgen_cistron_ids = load_subgen_ids(pergene_path, category)
	B, kept, missing = subset_to_subgen(
		M, row_ids, subgen_gene_ids, subgen_cistron_ids)
	n_subgen = len(subgen_gene_ids)
	print('Matrix: %d rows x %d cells. Subgen genes matched: %d/%d (%d missing).'
		% (M.shape[0], M.shape[1], len(kept), n_subgen, len(missing)))

	res = analyze(B)
	rng = np.random.default_rng(seed)
	null_vals = curveball_null(B, n_null, rng) if n_null > 0 else np.array([])
	evr, scores = pca_cells(B)
	wb = within_between_lineage(B, seeds) if seeds is not None else None

	prefix = os.path.join(out_dir, 'subgen_set_uniqueness')
	plot_count_hist(res, prefix + '_count_hist.png')
	plot_jaccard(res, prefix + '_jaccard.png', null_vals)
	plot_gene_freq(res, prefix + '_gene_freq.png')
	plot_accumulation(res, prefix + '_accumulation.png')
	plot_pca(scores, evr, res['per_cell_count'], prefix + '_pca.png')
	if wb is not None:
		plot_within_between(wb, prefix + '_within_between_lineage.png')

	cs = res['count_stats']
	null_z = (float('nan') if not len(null_vals)
		else (res['obs_mean_jaccard'] - null_vals.mean()) / (null_vals.std() + 1e-12))
	summary = {
		'inputs': {'matrix': matrix_path, 'pergene': pergene_path,
			'category': category, 'n_cells': res['n_cells'],
			'n_subgen_genes': res['n_genes'], 'n_missing_from_matrix': len(missing)},
		'number_per_cell': cs,
		'uniqueness': {
			'n_distinct_sets': res['n_distinct_sets'],
			'n_cells': res['n_cells'],
			'largest_identical_group': res['largest_identical_group'],
			'n_cells_with_unique_set': res['n_cells_with_unique_set']},
		'set_overlap': {
			'observed_mean_jaccard': res['obs_mean_jaccard'],
			'independence_null_mean_jaccard': res['null_mean_jaccard_independence'],
			'fixed_margins_null_mean': (None if not len(null_vals)
				else float(null_vals.mean())),
			'fixed_margins_null_sd': (None if not len(null_vals)
				else float(null_vals.std())),
			'fixed_margins_z': None if not len(null_vals) else float(null_z)},
		'coverage': {
			'total_subgen_genes': res['n_genes'],
			'covered_by_population': res['total_covered'],
			'covered_fraction': res['total_covered'] / res['n_genes']},
		'lineage_overlap': None if wb is None else {
			'within_lineage_mean_jaccard': wb['within_mean'],
			'between_lineage_mean_jaccard': wb['between_mean']},
		'pca_top_explained_variance': [float(x) for x in evr[:5]],
		}
	with open(prefix + '_summary.json', 'w') as f:
		json.dump(summary, f, indent=2)
	print('wrote', prefix + '_summary.json')

	print('\n===== SUMMARY =====')
	print('Cells: %d   subgen genes: %d' % (res['n_cells'], res['n_genes']))
	print('Number per cell: %.1f +/- %.1f  (CV %.2f; null CV %.2f; dispersion %.2f)'
		% (cs['mean'], cs['sd'], cs['cv'], cs['null_sd'] / cs['null_mean'],
			cs['dispersion_ratio']))
	print('Distinct sets: %d/%d  (largest identical group: %d)'
		% (res['n_distinct_sets'], res['n_cells'], res['largest_identical_group']))
	print('Mean pairwise Jaccard: %.3f  (independence null %.3f; fixed-margins %s)'
		% (res['obs_mean_jaccard'], res['null_mean_jaccard_independence'],
			'NA' if not len(null_vals) else '%.3f +/- %.3f, z=%.1f'
			% (null_vals.mean(), null_vals.std(), null_z)))
	print('Population covers %d/%d subgen genes (%.0f%%); PC1/PC2 explain %.1f%%/%.1f%%'
		% (res['total_covered'], res['n_genes'],
			100 * res['total_covered'] / res['n_genes'], 100 * evr[0], 100 * evr[1]))
	if wb is not None:
		print('Lineage overlap: within %.3f vs between %.3f'
			% (wb['within_mean'], wb['between_mean']))
	return summary


def main():
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument('matrix')
	ap.add_argument('pergene')
	ap.add_argument('out_dir')
	ap.add_argument('--category', default='subgen')
	ap.add_argument('--n-null', type=int, default=100)
	ap.add_argument('--seed', type=int, default=0)
	a = ap.parse_args()
	run(a.matrix, a.pergene, a.out_dir, a.category, a.n_null, a.seed)


if __name__ == '__main__':
	main()
