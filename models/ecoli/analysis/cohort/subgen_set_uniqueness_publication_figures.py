"""
Publication-ready figures for the subgenerational per-cell expression story.

Produces the "clean" panels of the Figure-3 series, all sharing one visual style
(DejaVu Sans, a single blue accent, an amber mean line, no top titles, and
editable-text PDF/PNG/SVG so they resize in Illustrator):

  count_distribution -- distribution over cells of the NUMBER of subgenerational
      genes expressed per cell, with a dashed mean line (mean, SD in the legend).
  jaccard            -- distribution of pairwise Jaccard similarity over all cell
      pairs, with a dashed mean line (mean, SD in the legend).
  fig3B_heatmap      -- genes x cells binary heatmap (rows sorted by expression
      frequency) topped by a per-cell "% of subgen genes expressed" bar carrying a
      dashed mean line (mean, SD in the legend). Drawn in black-and-white and in
      color, for each requested number of sampled cells.

Inputs match subgen_set_uniqueness.py: the binary "expressed" matrix (the
definition-5 0/1 criterion) and the definition-5 per-gene table (rows with
category == subgen define the gene set).

The spread statistic is SD (cell-to-cell variability -- the width the figures
actually show), not SE: these panels make a heterogeneity claim, so they report
how much cells differ, not how precisely the mean is known.

Heatmap cell sampling is reproducible across machines: numpy's PCG64
(default_rng(seed)) draws from cells first sorted into a canonical order, so the
selection does not depend on input row order.

Usage
-----
  python subgen_set_uniqueness_publication_figures.py <binary_matrix.tsv> \
      <pergene.tsv> <output_dir> [--category subgen] [--n-cells 120 100 50] \
      [--seed 1]

Self-contained: numpy + matplotlib only.
"""

import os
import csv
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter


ACCENT = '#1667B8'      # single blue accent (matches panels 3C/D/E)
BLACK = '#2b2b2b'       # black-and-white "expressed" fill
MEAN = '#E8A33D'        # amber mean line (colorblind-safe against blue and black)
MUTED = '#5b6672'       # reference lines / secondary text
COMMA = FuncFormatter(lambda x, _: '{:,.0f}'.format(x))


# ------------------------------------------------------------------ loading

def _to01(values):
	return [1 if v not in ('', '0', '0.0', 'False') else 0 for v in values]


def load_binary_matrix(path):
	"""Load a binary expression matrix in either supported layout, returning
	(row_ids, cell_ids, M) with M an int8 array [n_genes, n_cells].

	- genes-as-rows: first column is the gene id, remaining columns one per cell.
	- cells-as-rows (the `_expressed.tsv` layout): columns are `seed`,
	  `generation`, then one per gene (header = gene id); each row is one cell.
	"""
	with open(path) as f:
		reader = csv.reader(f, delimiter='\t')
		header = next(reader)
		if header[0] == 'seed':                       # cells-as-rows
			row_ids = header[2:]
			cols, cell_ids = [], []
			for rec in reader:
				if not rec:
					continue
				cell_ids.append('seed%06d_gen%03d' % (int(rec[0]), int(rec[1])))
				cols.append(_to01(rec[2:]))
			M = np.array(cols, dtype=np.int8).T       # genes x cells
			return row_ids, cell_ids, M
		cell_ids = header[1:]                          # genes-as-rows
		row_ids, rows = [], []
		for rec in reader:
			if not rec:
				continue
			row_ids.append(rec[0])
			rows.append(_to01(rec[1:]))
	return row_ids, cell_ids, np.array(rows, dtype=np.int8)


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


def subset_to_subgen(M, row_ids, subgen_gene_ids, subgen_cistron_ids):
	"""Keep only the rows (genes) in the subgen set, auto-matching whichever id
	(gene_id or cistron_id) the matrix rows are keyed by. Returns (B, kept)."""
	n_gene = sum(1 for r in row_ids if r in subgen_gene_ids)
	n_cistron = sum(1 for r in row_ids if r in subgen_cistron_ids)
	target = subgen_gene_ids if n_gene >= n_cistron else subgen_cistron_ids
	index = {rid: i for i, rid in enumerate(row_ids)}
	kept = [rid for rid in row_ids if rid in target]
	return M[[index[rid] for rid in kept], :], kept


# ------------------------------------------------------------------ style / io

def _apply_style():
	plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13,
		'axes.labelsize': 14.5, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
		'axes.linewidth': 1.0, 'figure.facecolor': 'white',
		'axes.facecolor': 'white', 'savefig.facecolor': 'white',
		'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})


def _despine(ax):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)


def _save(fig, out_dir, name, dpi):
	for ext in ('pdf', 'png', 'svg'):
		fig.savefig(os.path.join(out_dir, name + '.' + ext), dpi=dpi,
			bbox_inches='tight')
	plt.close(fig)


# ------------------------------------------------------------------ figures

def figure_count_distribution(B, n_total, out_dir):
	"""Distribution over cells of subgen genes expressed per cell."""
	k = B.sum(axis=0)
	mean_k, sd_k = k.mean(), k.std()
	fig, ax = plt.subplots(figsize=(6.0, 4.2))
	ax.hist(k, bins=35, color=ACCENT, edgecolor='white', linewidth=0.4)
	ax.axvline(mean_k, color=MEAN, lw=2.0, ls='--',
		label='mean = %.0f\nSD = %.0f' % (mean_k, sd_k))
	ax.set_xlabel('Subgenerational genes expressed per cell')
	ax.set_ylabel('Number of cells')
	ax.legend(frameon=False, fontsize=11, loc='upper right', handlelength=1.4,
		borderaxespad=0.3)
	_despine(ax)
	ax.yaxis.set_major_formatter(COMMA)
	_save(fig, out_dir, 'count_distribution', dpi=300)
	print('count_distribution: n=%d cells, range %d-%d, mean %.1f SD %.1f'
		% (B.shape[1], k.min(), k.max(), mean_k, sd_k))


def figure_jaccard(B, out_dir):
	"""Distribution of pairwise Jaccard similarity over all cell pairs."""
	Bi = B.astype(np.int32)
	inter = Bi.T @ Bi
	kk = np.diag(inter).astype(np.float64)
	iu = np.triu_indices(B.shape[1], 1)
	union = kk[iu[0]] + kk[iu[1]] - inter[iu]
	jac = inter[iu] / np.where(union > 0, union, 1)
	mean_j, sd_j = jac.mean(), jac.std()
	fig, ax = plt.subplots(figsize=(5.4, 4.0))
	ax.hist(jac, bins=60, color=ACCENT, edgecolor='white', linewidth=0.2)
	ax.axvline(mean_j, color=MEAN, lw=2.0, ls='--',
		label='mean = %.3f\nSD = %.3f' % (mean_j, sd_j))
	ax.set_xlabel('Pairwise Jaccard similarity of cell pairs')
	ax.set_ylabel('Number of cell pairs')
	ax.legend(frameon=False, fontsize=11, loc='upper right', handlelength=1.4,
		borderaxespad=0.3)
	_despine(ax)
	ax.yaxis.set_major_formatter(COMMA)
	_save(fig, out_dir, 'jaccard', dpi=300)
	print('jaccard: mean %.4f SD %.4f over %d cell pairs'
		% (mean_j, sd_j, jac.size))


def figure_heatmap(G, n_total, n_cells, seed, on_color, name, out_dir):
	"""Genes x cells binary heatmap with a per-cell %-expressed bar + mean line.

	`G` is the subgen matrix with rows already sorted by expression frequency and
	columns already in canonical (reproducible) order.
	"""
	n_sub, n_avail = G.shape
	rng = np.random.default_rng(seed)
	sel = np.sort(rng.choice(n_avail, size=n_cells, replace=False))
	S = G[:, sel]
	pct = 100 * S.sum(axis=0) / n_sub
	fw = 4.6 if n_cells == 120 else 1.2 + n_cells * 0.029

	fig = plt.figure(figsize=(fw, 5.6))
	gs = GridSpec(2, 1, height_ratios=[1, 6], hspace=0.06)
	axb = fig.add_subplot(gs[0])
	axh = fig.add_subplot(gs[1], sharex=axb)

	axb.bar(np.arange(n_cells), pct, width=1.0, color=on_color, linewidth=0)
	axb.axhline(100, color=MUTED, lw=1.2, ls=(0, (4, 3)))
	axb.text(n_cells * 0.99, 100, 'all {:,}'.format(n_sub), ha='right',
		va='bottom', fontsize=10, color=MUTED)
	axb.axhline(pct.mean(), color=MEAN, lw=1.8, ls='--',
		label='mean %.1f%%\nSD %.1f%%' % (pct.mean(), pct.std()))
	axb.set_ylim(0, 108)
	axb.set_yticks([0, 50, 100])
	axb.set_yticklabels(['0', '50', '100'])
	axb.set_ylabel('% of\nsubgen\ngenes', fontsize=11, rotation=0, ha='right',
		va='center', labelpad=6)
	axb.set_xlim(-0.5, n_cells - 0.5)
	axb.tick_params(labelbottom=False, bottom=False)
	for sp in ('top', 'right', 'bottom'):
		axb.spines[sp].set_visible(False)
	axb.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False,
		fontsize=9, handlelength=1.4, borderaxespad=0.0)

	axh.imshow(S, aspect='auto', interpolation='nearest',
		cmap=ListedColormap(['white', on_color]), vmin=0, vmax=1)
	axh.set_xticks([])
	axh.set_yticks([])
	for s in axh.spines.values():
		s.set_color('#999')
		s.set_linewidth(0.8)
	axh.set_xlabel('Cells (n = %d)' % n_cells)
	axh.set_ylabel('Subgenerational genes (n = {:,})'.format(n_sub))
	leg = [Patch(facecolor=on_color, label='Expressed'),
		Patch(facecolor='white', edgecolor='#999', label='Not expressed')]
	axh.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.05),
		ncol=2, frameon=False, fontsize=11, handlelength=1.1, handleheight=1.1,
		columnspacing=1.4)

	_save(fig, out_dir, name, dpi=600)
	print('%-30s N=%3d  pct min %.1f mean %.1f SD %.1f max %.1f'
		% (name, n_cells, pct.min(), pct.mean(), pct.std(), pct.max()))


# ------------------------------------------------------------------ main

def main():
	ap = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument('binary_matrix', help='binary "expressed" matrix TSV')
	ap.add_argument('pergene', help='definition-5 per-gene table TSV')
	ap.add_argument('output_dir', help='directory for the figures')
	ap.add_argument('--category', default='subgen',
		help='per-gene category defining the gene set (default: subgen)')
	ap.add_argument('--n-cells', type=int, nargs='+', default=[120],
		help='cell counts to draw the heatmap for (default: 120)')
	ap.add_argument('--seed', type=int, default=1,
		help='PCG64 seed for reproducible cell sampling (default: 1)')
	args = ap.parse_args()

	os.makedirs(args.output_dir, exist_ok=True)
	_apply_style()

	row_ids, cell_ids, M = load_binary_matrix(args.binary_matrix)
	gene_ids, cistron_ids = load_subgen_ids(args.pergene, args.category)
	B, kept = subset_to_subgen(M, row_ids, gene_ids, cistron_ids)
	n_sub, n_cells_total = B.shape
	print('matrix: %d subgen genes x %d cells (category=%s)'
		% (n_sub, n_cells_total, args.category))

	figure_count_distribution(B, n_cells_total, args.output_dir)
	figure_jaccard(B, args.output_dir)

	# Canonical, order-independent column order for reproducible sampling: the
	# zero-padded `seed%06d_gen%03d` labels sort lexicographically as (seed, gen).
	order = np.argsort(cell_ids)
	G = B[:, order]
	G = G[np.argsort(-G.sum(axis=1))]        # rows by expression frequency
	for n in args.n_cells:
		for on_color, tag in ((BLACK, 'bw'), (ACCENT, 'color')):
			figure_heatmap(G, n_cells_total, n, args.seed, on_color,
				'fig3B_heatmap_%s_%dcells' % (tag, n), args.output_dir)


if __name__ == '__main__':
	main()
