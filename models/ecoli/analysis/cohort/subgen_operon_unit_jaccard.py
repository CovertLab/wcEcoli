"""
Operon-unit robustness check for the subgen per-cell uniqueness / Jaccard result.

Recomputes the mean pairwise cell-to-cell Jaccard (and its independence and
fixed-margins nulls) at TWO levels and prints them side by side:

  cistron level -- the original analysis: one row per subgen gene (cistron).
  operon-unit level -- co-transcribed subgen cistrons collapsed into one unit
      (a unit is "expressed" in a cell if ANY of its subgen cistrons fired).
      Units are the connected components of the graph that links two subgen
      cistrons whenever they share a transcription unit (handles overlapping TUs).

If operons were inflating cell-to-cell similarity, collapsing them would pull the
observed Jaccard down toward its null. The expected result (the mean Jaccard is
governed by per-unit population frequency and is invariant to within-cell
co-transcription) is that observed stays ~equal to the null at both levels.

The Jaccard, independence-null, and curveball-null implementations are copied
verbatim from subgen_set_uniqueness.py so the two levels are directly comparable.

Inputs
------
  binary_matrix.tsv : the def-5 "expressed" 0/1 matrix (cells-as-rows: seed,
      generation, then one column per gene).
  pergene.tsv       : def-5 per-gene table; category == subgen defines the set.
  sim_data.cPickle  : SimulationDataEcoli pickle (operons ON) for the cistron->TU
      mapping.

Usage
-----
  python subgen_operon_unit_jaccard.py <binary_matrix.tsv> <pergene.tsv> \
      <sim_data.cPickle> <output_dir> [--category subgen] [--n-null 10] [--seed 0]

Self-contained: numpy + scipy + the sim_data pickle.
"""

import os
import sys
import csv
import json
import pickle
import argparse

import numpy as np
import scipy.sparse as sp

# Standalone CLI: make `models.ecoli.analysis.cohort.subgen_common` importable even
# when this file is run by path from an arbitrary cwd with no PYTHONPATH set. Only
# subgen_common is imported, and it keeps its `wholecell` imports lazy, so this
# script still needs nothing beyond numpy (+ matplotlib/scipy where used).
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
	os.path.dirname(os.path.abspath(__file__))))))
if _REPO_ROOT not in sys.path:
	sys.path.insert(0, _REPO_ROOT)
from models.ecoli.analysis.cohort import subgen_common as sc


# ------------------------------------------------------------------ loading

def _to01(values):
	return [1 if v not in ('', '0', '0.0', 'False') else 0 for v in values]


def load_binary_matrix(path):
	"""Load the cells-as-rows expressed matrix. Returns (gene_ids, M) with M an
	int8 array [n_genes, n_cells]."""
	with open(path) as f:
		reader = csv.reader(f, delimiter='\t')
		header = next(reader)
		gene_ids = header[2:]
		cols = []
		for rec in reader:
			if not rec:
				continue
			cols.append(_to01(rec[2:]))
	return gene_ids, np.array(cols, dtype=np.int8).T


def load_subgen(path, category='subgen'):
	"""Return dict gene_id -> cistron_id for rows with the given category.

	Delegates to sc.load_def5_categories, which accepts either canonical per-gene
	table's column spelling. See subgen_set_uniqueness.py for the details.
	"""
	sel = sc.load_def5_categories(path, category=category)
	return dict(zip(sel['gene_ids'], sel['cistron_ids']))


# --------------------------------------------- jaccard + nulls (from uniqueness)

def pairwise_jaccard_mean(B):
	Bi = B.astype(np.int32)
	inter = Bi.T @ Bi
	k = np.diag(inter).astype(np.float64)
	union = k[:, None] + k[None, :] - inter
	iu = np.triu_indices(inter.shape[0], k=1)
	u = union[iu]
	nz = u > 0
	jac = np.zeros(u.shape)
	jac[nz] = inter[iu][nz] / u[nz]
	return float(jac.mean())


def independence_null_jaccard(B):
	"""Analytic independence null: Sp^2 / (2*Sp - Sp^2) with p the per-row freq."""
	p = B.mean(axis=1)
	exp_shared = float((p ** 2).sum())
	exp_size = float(p.sum())
	return exp_shared / (2 * exp_size - exp_shared) if exp_size else float('nan')


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
	n_ones = int(B.sum())
	n_iter = max(2000, 5 * n_ones)
	vals = []
	cur = B.copy()
	for _ in range(n_null):
		cur = curveball(cur, n_iter, rng)
		vals.append(pairwise_jaccard_mean(cur))
	return np.array(vals)


# ------------------------------------------------------------------ operon units

class _UF:
	def __init__(self, n):
		self.p = list(range(n))

	def find(self, x):
		while self.p[x] != x:
			self.p[x] = self.p[self.p[x]]
			x = self.p[x]
		return x

	def union(self, a, b):
		ra, rb = self.find(a), self.find(b)
		if ra != rb:
			self.p[rb] = ra


def operon_units(sub_cistron_rows, M):
	"""Group subgen cistron rows (indices into sim_data) into operon units:
	connected components over "share a TU". Returns a list of component labels,
	one per entry of sub_cistron_rows."""
	pos = {c: i for i, c in enumerate(sub_cistron_rows)}
	uf = _UF(len(sub_cistron_rows))
	for t in range(M.shape[1]):
		members = [pos[c] for c in np.where(M[:, t])[0] if c in pos]
		for m in members[1:]:
			uf.union(members[0], m)
	roots = [uf.find(i) for i in range(len(sub_cistron_rows))]
	relabel = {r: i for i, r in enumerate(sorted(set(roots)))}
	return np.array([relabel[r] for r in roots])


def collapse(B, labels):
	"""OR-collapse rows of B sharing a label into one unit row."""
	n_units = labels.max() + 1
	U = np.zeros((n_units, B.shape[1]), dtype=np.int8)
	for u in range(n_units):
		U[u] = (B[labels == u].sum(axis=0) > 0).astype(np.int8)
	return U


# ------------------------------------------------------------------ reporting

def level_stats(B, n_null, rng):
	k = B.sum(axis=0).astype(np.float64)
	p = B.mean(axis=1)
	var_null = float((p * (1 - p)).sum())
	cb = curveball_null(B, n_null, rng) if n_null > 0 else np.array([])
	return {
		'n_rows': int(B.shape[0]),
		'n_ones': int(B.sum()),
		'count_mean': float(k.mean()),
		'count_cv': float(k.std(ddof=1) / k.mean()),
		'dispersion_ratio': float(k.var() / var_null) if var_null else float('nan'),
		'obs_jaccard': pairwise_jaccard_mean(B),
		'indep_null_jaccard': independence_null_jaccard(B),
		'curveball_mean': float(cb.mean()) if cb.size else None,
		'curveball_sd': float(cb.std()) if cb.size else None,
		}


def main():
	ap = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument('binary_matrix')
	ap.add_argument('pergene')
	ap.add_argument('sim_data')
	ap.add_argument('output_dir')
	ap.add_argument('--category', default='subgen')
	ap.add_argument('--n-null', type=int, default=10,
		help='curveball shuffles per level (0 to skip; default 10)')
	ap.add_argument('--seed', type=int, default=0)
	args = ap.parse_args()
	os.makedirs(args.output_dir, exist_ok=True)
	rng = np.random.default_rng(args.seed)

	gene_ids, M_all = load_binary_matrix(args.binary_matrix)
	subgen = load_subgen(args.pergene, args.category)          # gene_id -> cistron_id
	sub_cols = [i for i, g in enumerate(gene_ids) if g in subgen]
	sub_gene_ids = [gene_ids[i] for i in sub_cols]
	B_cistron = M_all[sub_cols, :]                            # subgen cistron x cell

	with open(args.sim_data, 'rb') as f:
		sim_data = pickle.load(f)
	tx = sim_data.process.transcription
	cis_ids = np.asarray(tx.cistron_data['id'])
	cis_idx = {c: i for i, c in enumerate(cis_ids)}
	Msp = tx.cistron_tu_mapping_matrix
	Msp = (Msp.toarray() if sp.issparse(Msp) else np.asarray(Msp)) > 0
	# sim_data cistron rows for each subgen gene, in B_cistron row order.
	sub_cistron_rows = np.array([cis_idx[subgen[g]] for g in sub_gene_ids])

	labels = operon_units(sub_cistron_rows, Msp)
	B_unit = collapse(B_cistron, labels)
	print('subgen cistrons: %d  ->  operon-units: %d'
		% (B_cistron.shape[0], B_unit.shape[0]))

	cis_stats = level_stats(B_cistron, args.n_null, rng)
	unit_stats = level_stats(B_unit, args.n_null, rng)

	def fmt(v):
		return '%.4f' % v if isinstance(v, float) else str(v)
	keys = ['n_rows', 'n_ones', 'count_mean', 'count_cv', 'dispersion_ratio',
		'obs_jaccard', 'indep_null_jaccard', 'curveball_mean', 'curveball_sd']
	print('\n%-20s %14s %14s' % ('metric', 'cistron', 'operon-unit'))
	for key in keys:
		print('%-20s %14s %14s'
			% (key, fmt(cis_stats[key]), fmt(unit_stats[key])))

	with open(os.path.join(args.output_dir,
			'subgen_operon_unit_jaccard_summary.tsv'), 'w') as f:
		w = csv.writer(f, delimiter='\t')
		w.writerow(['metric', 'cistron', 'operon_unit'])
		for key in keys:
			w.writerow([key, cis_stats[key], unit_stats[key]])
	with open(os.path.join(args.output_dir,
			'subgen_operon_unit_jaccard_run_metadata.json'), 'w') as f:
		json.dump({'inputs': {'binary_matrix': args.binary_matrix,
			'pergene': args.pergene, 'sim_data': args.sim_data,
			'category': args.category, 'n_null': args.n_null, 'seed': args.seed},
			'cistron': cis_stats, 'operon_unit': unit_stats}, f, indent=2)
	print('\nWrote summary + metadata to %s' % args.output_dir)


if __name__ == '__main__':
	main()
