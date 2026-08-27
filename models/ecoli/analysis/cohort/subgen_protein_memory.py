"""
Subgenerational transcription vs. protein absence ("protein memory") figures.

Tests the hypothesis that frequently-transcribed subgenerational genes rarely run
out of protein, because long protein half-lives let a cell carry inherited protein
through a generation that fails to transcribe. For each subgen gene it plots how
sub-generationally it is transcribed against how often its protein is absent.

Both axes are per-lineage *rates* aggregated to a 95% CI (the def-5 estimator), so
they are directly comparable and scale honestly with more seeds:

  y  protein-absence rate  = expected fraction of cell-cycle TIME with zero protein
      copies. Estimated per cell cycle (time-weighted), averaged within a lineage,
      then averaged across lineages. Read as: freeze a random cell at a random
      instant -- chance of seeing zero copies. Source: the extractor's
      `_frac_protein_zero_per_cell.tsv` (see subgen_raw_extract.py).
  x1 transcript-off frequency = fraction of generations with zero completed
      transcripts (same units as y, so the y = x diagonal reads as "memory").
  x2 def-5 rate = mean completed transcripts per generation (log x-axis).
      Both x from `_synth_per_cell.tsv`.

The independent replicate is the lineage (cells within a lineage share inherited
protein), so every rate is estimated per lineage and CIs come from the spread
across lineages -- never by pooling cells and never by counting "ever/never absent"
genes (that count would drift with sample size like presence-frequency-at-1).

Inputs
------
  frac_zero.tsv : `_frac_protein_zero_per_cell.tsv` -- wide, `seed`, `generation`,
      `is_successful_lineage` + one column per gene (fraction of cycle at zero).
  synth.tsv     : `_synth_per_cell.tsv` -- same layout, completed transcripts per
      generation.
  genes.tsv     : `_genes.tsv` -- gene_id, cistron_id, monomer_id (column key).
  pergene.tsv   : def-5 per-gene table; genes with category == subgen define the
      gene set (matched on gene_id).

Usage
-----
  python subgen_protein_memory.py <frac_zero.tsv> <synth.tsv> <genes.tsv> \
      <pergene.tsv> <output_dir> [--category subgen] [--error-bars]

Self-contained: numpy + matplotlib only.
"""

import os
import sys
import csv
import json
import argparse
import subprocess

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Standalone CLI: make `models.ecoli.analysis.cohort.subgen_common` importable even
# when this file is run by path from an arbitrary cwd with no PYTHONPATH set. Only
# subgen_common is imported, and it keeps its `wholecell` imports lazy, so this
# script still needs nothing beyond numpy + matplotlib.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
	os.path.dirname(os.path.abspath(__file__))))))
if _REPO_ROOT not in sys.path:
	sys.path.insert(0, _REPO_ROOT)
from models.ecoli.analysis.cohort import subgen_common as sc


ACCENT = '#1667B8'      # single blue accent (matches the other subgen figures)
MEAN = '#E8A33D'        # amber reference line
INK = '#1b2530'
MUTED = '#5b6672'
CI_Z = 1.96             # 95% CI from the SE (normal approx; matches subgen_common)
N_META = 3              # seed, generation, is_successful_lineage
COMMA = FuncFormatter(lambda x, _: '{:,.0f}'.format(x))


# ------------------------------------------------------------------ loading

def load_per_cell(path):
	"""Load a wide per-cell matrix (seed, generation, is_successful_lineage +
	one column per gene). Returns (seeds, gens, is_succ, gene_ids, M) with M a
	float array [n_cells, n_genes]."""
	with open(path) as f:
		reader = csv.reader(f, delimiter='\t')
		header = next(reader)
		gene_ids = header[N_META:]
		seeds, gens, is_succ, rows = [], [], [], []
		for rec in reader:
			if not rec:
				continue
			seeds.append(int(rec[0]))
			gens.append(int(rec[1]))
			is_succ.append(int(rec[2]))
			rows.append([float(v) for v in rec[N_META:]])
	return (np.array(seeds), np.array(gens), np.array(is_succ), gene_ids,
		np.array(rows, dtype=np.float64))


def load_subgen_gene_ids(path, category='subgen'):
	"""Return the set of gene_ids whose category matches (from the def-5 table).

	Delegates to sc.load_def5_categories, which accepts either canonical per-gene
	table's column spelling. See subgen_set_uniqueness.py for the details.
	"""
	return set(sc.load_def5_categories(path, category=category)['gene_ids'])


# ------------------------------------------------------------------ statistics

def per_lineage_means(values, seeds):
	"""Mean of `values` (n_cells x n_genes) within each lineage (seed).

	Returns (lineage_ids, L) with L an [n_lineages, n_genes] array of per-lineage
	means -- each lineage contributes one unweighted data point per gene."""
	lineage_ids = np.unique(seeds)
	L = np.empty((lineage_ids.size, values.shape[1]), dtype=np.float64)
	for i, s in enumerate(lineage_ids):
		L[i] = values[seeds == s].mean(axis=0)
	return lineage_ids, L


def summarize(L):
	"""Across-lineage mean, std, SE and 95% CI for an [n_lineages, n_genes]
	matrix of per-lineage estimates."""
	n = L.shape[0]
	mean = L.mean(axis=0)
	std = L.std(axis=0, ddof=1) if n > 1 else np.zeros(L.shape[1])
	se = std / np.sqrt(n)
	return {'mean': mean, 'std': std, 'se': se,
		'ci_lower': np.clip(mean - CI_Z * se, 0, None),
		'ci_upper': mean + CI_Z * se}


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


def _save(fig, out_dir, name, dpi=300):
	for ext in ('pdf', 'png', 'svg'):
		fig.savefig(os.path.join(out_dir, name + '.' + ext), dpi=dpi,
			bbox_inches='tight')
	plt.close(fig)


def _git_info():
	repo = os.path.dirname(os.path.abspath(__file__))
	try:
		out = subprocess.check_output(
			['git', '-C', repo, 'rev-parse', 'HEAD'],
			stderr=subprocess.DEVNULL).decode().strip()
		branch = subprocess.check_output(
			['git', '-C', repo, 'rev-parse', '--abbrev-ref', 'HEAD'],
			stderr=subprocess.DEVNULL).decode().strip()
		return {'git_hash': out, 'git_branch': branch}
	except Exception as e:
		return {'git_hash': None, 'git_branch': None, 'error': str(e)}


# ------------------------------------------------------------------ figures

def figure_vs_transcript_off(x, y, out_dir, error_bars):
	"""Scatter: transcript-off frequency (x) vs protein-absence rate (y)."""
	fig, ax = plt.subplots(figsize=(5.6, 5.4))
	# Reference diagonal: genes below it keep protein even when transcription is
	# off (the memory signal).
	ax.plot([0, 1], [0, 1], color=MUTED, lw=1.4, ls='--', zorder=1)
	if error_bars:
		ax.errorbar(x['mean'], y['mean'],
			xerr=CI_Z * x['se'], yerr=CI_Z * y['se'], fmt='none',
			ecolor=ACCENT, elinewidth=0.4, alpha=0.25, zorder=2)
	ax.scatter(x['mean'], y['mean'], s=14, color=ACCENT, alpha=0.4,
		edgecolor='none', zorder=3)
	ax.set_xlim(-0.02, 1.02)
	ax.set_ylim(-0.02, 1.02)
	ax.set_aspect('equal')
	ax.set_xlabel('Fraction of cells with\nno transcription event')
	ax.set_ylabel('Fraction of cell-cycle time\nwith absent protein')
	_despine(ax)
	_save(fig, out_dir, 'subgen_protein_memory_vs_transcript_off')


def figure_vs_def5_rate(x, y, out_dir, error_bars):
	"""Scatter: def-5 rate on a log x-axis vs protein-absence rate (y)."""
	rate = x['mean']
	pos = rate > 0
	fig, ax = plt.subplots(figsize=(6.0, 4.6))
	if error_bars:
		ax.errorbar(rate[pos], y['mean'][pos],
			yerr=CI_Z * y['se'][pos], fmt='none', ecolor=ACCENT,
			elinewidth=0.4, alpha=0.25, zorder=2)
	ax.scatter(rate[pos], y['mean'][pos], s=14, color=ACCENT, alpha=0.4,
		edgecolor='none', zorder=3)
	ax.set_xscale('log')
	ax.set_ylim(-0.02, 1.02)
	ax.set_xlabel('Def-5 transcription rate\n(mean completed transcripts per generation)')
	ax.set_ylabel('Protein-absence rate\n(fraction of cell-cycle time at 0 copies)')
	_despine(ax)
	_save(fig, out_dir, 'subgen_protein_memory_vs_def5_rate')


def figure_binned(x_off, y, out_dir, n_bins=10):
	"""Binned trend: mean protein-absence rate per transcript-off-frequency
	decile, with a 95% CI band over the genes in each bin."""
	xv = x_off['mean']
	yv = y['mean']
	edges = np.quantile(xv, np.linspace(0, 1, n_bins + 1))
	edges[-1] += 1e-9
	centers, means, los, his = [], [], [], []
	for lo, hi in zip(edges[:-1], edges[1:]):
		m = (xv >= lo) & (xv < hi)
		if m.sum() < 2:
			continue
		vals = yv[m]
		se = vals.std(ddof=1) / np.sqrt(vals.size)
		centers.append(xv[m].mean())
		means.append(vals.mean())
		los.append(vals.mean() - CI_Z * se)
		his.append(vals.mean() + CI_Z * se)
	centers, means = np.array(centers), np.array(means)
	fig, ax = plt.subplots(figsize=(6.0, 4.4))
	ax.plot([0, 1], [0, 1], color=MUTED, lw=1.2, ls='--', zorder=1)
	ax.fill_between(centers, los, his, color=ACCENT, alpha=0.2, zorder=2)
	ax.plot(centers, means, color=ACCENT, lw=2.2, marker='o', ms=5, zorder=3)
	ax.set_xlim(-0.02, 1.02)
	ax.set_ylim(-0.02, 1.02)
	ax.set_xlabel('Transcript-off frequency (decile bins)')
	ax.set_ylabel('Mean protein-absence rate')
	_despine(ax)
	_save(fig, out_dir, 'subgen_protein_memory_binned')


# ------------------------------------------------------------------ output

def write_pergene(path, gene_ids, cistron_ids, n_lineages, y, x_off, x_rate):
	with open(path, 'w') as f:
		w = csv.writer(f, delimiter='\t')
		cols = ['gene_id', 'cistron_id', 'n_lineages']
		for tag in ('protein_zero_frac', 'transcript_off_freq', 'def5_rate'):
			cols += ['%s_mean' % tag, '%s_se' % tag,
				'%s_ci_lower' % tag, '%s_ci_upper' % tag]
		w.writerow(cols)
		for i, g in enumerate(gene_ids):
			row = [g, cistron_ids.get(g, ''), n_lineages]
			for s in (y, x_off, x_rate):
				row += ['%.6g' % s['mean'][i], '%.6g' % s['se'][i],
					'%.6g' % s['ci_lower'][i], '%.6g' % s['ci_upper'][i]]
			w.writerow(row)


# ------------------------------------------------------------------ main

def main():
	ap = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument('frac_zero', help='_frac_protein_zero_per_cell.tsv')
	ap.add_argument('synth', help='_synth_per_cell.tsv')
	ap.add_argument('genes', help='_genes.tsv (gene_id, cistron_id, monomer_id)')
	ap.add_argument('pergene', help='def-5 per-gene table (defines the gene set)')
	ap.add_argument('output_dir', help='directory for the figures + tables')
	ap.add_argument('--category', default='subgen',
		help='per-gene category defining the gene set (default: subgen)')
	ap.add_argument('--error-bars', action='store_true',
		help='draw faint 95%% CI whiskers on the scatter plots')
	args = ap.parse_args()

	os.makedirs(args.output_dir, exist_ok=True)
	_apply_style()

	seeds_z, gens_z, succ_z, gene_ids_z, fz = load_per_cell(args.frac_zero)
	seeds_s, gens_s, succ_s, gene_ids_s, synth = load_per_cell(args.synth)
	if gene_ids_z != gene_ids_s:
		raise ValueError('frac_zero and synth have different gene columns')
	if not (np.array_equal(seeds_z, seeds_s) and np.array_equal(gens_z, gens_s)):
		raise ValueError('frac_zero and synth are not row-aligned '
			'(same extraction run required)')

	# Restrict to successful lineages, then to the subgen gene set.
	keep = succ_z == 1
	seeds = seeds_z[keep]
	fz = fz[keep]
	synth = synth[keep]

	cistron_by_gene = {}
	with open(args.genes) as f:
		for r in csv.DictReader(f, delimiter='\t'):
			cistron_by_gene[r['gene_id']] = r['cistron_id']
	subgen = load_subgen_gene_ids(args.pergene, args.category)
	sub_idx = [i for i, g in enumerate(gene_ids_z) if g in subgen]
	sub_gene_ids = [gene_ids_z[i] for i in sub_idx]
	n_missing = len(subgen) - len(sub_idx)
	fz = fz[:, sub_idx]
	synth = synth[:, sub_idx]
	print('cells (successful) = %d, subgen genes = %d (%d in set not in matrix)'
		% (fz.shape[0], len(sub_idx), n_missing))

	# Per lineage, then across lineages, for each axis.
	lineage_ids, L_y = per_lineage_means(fz, seeds)
	_, L_xoff = per_lineage_means((synth == 0).astype(np.float64), seeds)
	_, L_xrate = per_lineage_means(synth, seeds)
	y = summarize(L_y)
	x_off = summarize(L_xoff)
	x_rate = summarize(L_xrate)
	n_lineages = lineage_ids.size
	print('lineages = %d; protein-absence rate mean %.3f (range %.3f-%.3f)'
		% (n_lineages, y['mean'].mean(), y['mean'].min(), y['mean'].max()))
	below = np.mean(y['mean'] < x_off['mean'])
	print('fraction of subgen genes below the y=x diagonal (memory): %.1f%%'
		% (100 * below))

	write_pergene(
		os.path.join(args.output_dir, 'subgen_protein_memory_pergene.tsv'),
		sub_gene_ids, cistron_by_gene, n_lineages, y, x_off, x_rate)
	figure_vs_transcript_off(x_off, y, args.output_dir, args.error_bars)
	figure_vs_def5_rate(x_rate, y, args.output_dir, args.error_bars)
	figure_binned(x_off, y, args.output_dir)

	with open(os.path.join(args.output_dir,
			'subgen_protein_memory_run_metadata.json'), 'w') as f:
		json.dump({
			'inputs': {'frac_zero': args.frac_zero, 'synth': args.synth,
				'genes': args.genes, 'pergene': args.pergene,
				'category': args.category},
			'n_lineages': int(n_lineages),
			'n_subgen_genes': len(sub_idx),
			'n_missing_from_matrix': int(n_missing),
			'n_successful_cells': int(fz.shape[0]),
			'ci_z': CI_Z, 'confidence': 0.95,
			'frac_below_diagonal': float(below),
			'git': _git_info(),
			}, f, indent=2)
	print('Wrote outputs to %s' % args.output_dir)


if __name__ == '__main__':
	main()
