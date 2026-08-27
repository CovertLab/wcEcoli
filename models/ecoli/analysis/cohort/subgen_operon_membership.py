"""
Operon membership of the subgenerational gene set.

Quantifies how many subgenerational genes are co-transcribed with other genes,
to assess whether operon structure could confound the per-cell uniqueness /
Jaccard-vs-chance analysis (co-transcribed cistrons are not independently
expressed). Uses the cistron -> transcription-unit (TU) 
mapping from sim_data;
a "polycistronic" TU (>= 2 cistrons) is an operon.

Note: the mean pairwise cell-to-cell Jaccard is governed by each gene's
population frequency and is invariant to within-cell co-transcription, so operons
do not bias that statistic even when they are common; this script documents how
common they are and lists the memberships for the operon-unit robustness check.

Inputs
------
  sim_data.cPickle : a SimulationDataEcoli pickle (kb/simData.cPickle) built with
      operons ON. Provides `process.transcription.cistron_tu_mapping_matrix`,
      `cistron_data['id'|'gene_id']`, and `rna_data['id']` (the TU ids).
  pergene.tsv      : the def-5 per-gene table; genes with category == subgen
      (matched on cistron_id) define the gene set.

Outputs (written to output_dir)
-------
  subgen_operon_membership_summary.tsv   metric, count, percent.
  subgen_operon_membership_per_gene.tsv  one row per subgen gene: its TU(s),
      whether it is polycistronic, and its operon siblings (all + subgen-only).
  subgen_operon_membership_operons.tsv   one row per TU containing >= 2 subgen
      genes: the TU, its cistron count, and the subgen genes in it.
  subgen_operon_membership_run_metadata.json  provenance + the summary counts.

Usage
-----
  python subgen_operon_membership.py <sim_data.cPickle> <pergene.tsv> \
      <output_dir> [--category subgen]

Self-contained: numpy + scipy + the sim_data pickle only.
"""

import os
import sys
import csv
import json
import pickle
import argparse
import subprocess

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


def load_subgen_cistrons(path, category='subgen'):
	"""Return an ordered list of cistron_ids whose category matches.

	Delegates to sc.load_def5_categories, which accepts either canonical per-gene
	table's column spelling. See subgen_set_uniqueness.py for the details.
	"""
	return sc.load_def5_categories(path, category=category)['cistron_ids']


def git_info():
	repo = os.path.dirname(os.path.abspath(__file__))
	try:
		run = lambda a: subprocess.check_output(
			['git', '-C', repo] + a, stderr=subprocess.DEVNULL).decode().strip()
		return {'git_hash': run(['rev-parse', 'HEAD']),
			'git_branch': run(['rev-parse', '--abbrev-ref', 'HEAD'])}
	except Exception as e:
		return {'git_hash': None, 'git_branch': None, 'error': str(e)}


def main():
	ap = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	ap.add_argument('sim_data', help='kb/simData.cPickle (operons ON)')
	ap.add_argument('pergene', help='def-5 per-gene table (defines the gene set)')
	ap.add_argument('output_dir', help='directory for the output tables')
	ap.add_argument('--category', default='subgen',
		help='per-gene category defining the gene set (default: subgen)')
	args = ap.parse_args()
	os.makedirs(args.output_dir, exist_ok=True)
	prefix = os.path.join(args.output_dir, 'subgen_operon_membership')

	with open(args.sim_data, 'rb') as f:
		sim_data = pickle.load(f)
	tx = sim_data.process.transcription
	cistron_ids = np.asarray(tx.cistron_data['id'])
	cistron_gene = np.asarray(tx.cistron_data['gene_id'])
	tu_ids = np.asarray(tx.rna_data['id'])
	M = tx.cistron_tu_mapping_matrix
	M = (M.toarray() if sp.issparse(M) else np.asarray(M)) > 0   # cistron x TU

	subgen = load_subgen_cistrons(args.pergene, args.category)
	idx = {c: i for i, c in enumerate(cistron_ids)}
	missing = [c for c in subgen if c not in idx]
	sub_rows = np.array([idx[c] for c in subgen if c in idx])
	subgen_mask = np.zeros(M.shape[0], dtype=bool)
	subgen_mask[sub_rows] = True
	n_sub = sub_rows.size

	tu_size = M.sum(axis=0)                             # cistrons per TU
	tu_n_subgen = (M & subgen_mask[:, None]).sum(axis=0)  # subgen cistrons per TU
	# Precompute the cistron members of every TU once.
	tu_members = [np.where(M[:, t])[0] for t in range(M.shape[1])]

	# Per-gene membership.
	pergene_rows = []
	in_operon = np.zeros(n_sub, dtype=bool)
	shares_subgen = np.zeros(n_sub, dtype=bool)
	for k, c in enumerate(sub_rows):
		tus = np.where(M[c])[0]
		members = np.unique(np.concatenate([tu_members[t] for t in tus]))
		siblings = members[members != c]
		sib_subgen = siblings[subgen_mask[siblings]]
		in_operon[k] = bool((tu_size[tus] >= 2).any())
		shares_subgen[k] = sib_subgen.size > 0
		pergene_rows.append({
			'gene_id': cistron_gene[c],
			'cistron_id': cistron_ids[c],
			'n_TUs': int(tus.size),
			'TU_ids': ';'.join(tu_ids[tus].tolist()),
			'max_operon_size': int(tu_size[tus].max()),
			'in_polycistronic_TU': int(in_operon[k]),
			'n_operon_siblings': int(siblings.size),
			'n_subgen_siblings': int(sib_subgen.size),
			'shares_TU_with_subgen': int(shares_subgen[k]),
			'subgen_sibling_genes': ';'.join(sorted(set(cistron_gene[sib_subgen].tolist()))),
			})

	# Summary table.
	mono = int((~in_operon).sum())
	sole = int((in_operon & ~shares_subgen).sum())
	summary = [
		('subgen_genes_total', n_sub, 100.0),
		('in_polycistronic_TU', int(in_operon.sum()), 100 * in_operon.mean()),
		('shares_TU_with_another_subgen', int(shares_subgen.sum()),
			100 * shares_subgen.mean()),
		('in_operon_sole_subgen_member', sole, 100 * sole / n_sub),
		('monocistronic_own_TU_only', mono, 100 * mono / n_sub),
	]
	with open(prefix + '_summary.tsv', 'w') as f:
		w = csv.writer(f, delimiter='\t')
		w.writerow(['metric', 'count', 'percent'])
		for name, count, pct in summary:
			w.writerow([name, count, '%.1f' % pct])
	print('%-34s %6s %7s' % ('metric', 'count', 'percent'))
	for name, count, pct in summary:
		print('%-34s %6d %6.1f%%' % (name, count, pct))

	# Per-gene table.
	with open(prefix + '_per_gene.tsv', 'w') as f:
		cols = ['gene_id', 'cistron_id', 'n_TUs', 'TU_ids', 'max_operon_size',
			'in_polycistronic_TU', 'n_operon_siblings', 'n_subgen_siblings',
			'shares_TU_with_subgen', 'subgen_sibling_genes']
		w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
		w.writeheader()
		for row in pergene_rows:
			w.writerow(row)

	# Per-operon table: TUs with >= 2 subgen genes.
	multi = np.where(tu_n_subgen >= 2)[0]
	multi = multi[np.argsort(-tu_n_subgen[multi])]
	with open(prefix + '_operons.tsv', 'w') as f:
		w = csv.writer(f, delimiter='\t')
		w.writerow(['TU_id', 'n_cistrons_total', 'n_subgen_cistrons',
			'subgen_genes'])
		for t in multi:
			members = tu_members[t]
			sub_members = members[subgen_mask[members]]
			w.writerow([tu_ids[t], int(tu_size[t]), int(tu_n_subgen[t]),
				';'.join(sorted(cistron_gene[sub_members].tolist()))])

	with open(prefix + '_run_metadata.json', 'w') as f:
		json.dump({
			'inputs': {'sim_data': args.sim_data, 'pergene': args.pergene,
				'category': args.category},
			'n_subgen_genes': int(n_sub),
			'n_missing_from_sim_data': len(missing),
			'n_TUs_total': int(M.shape[1]),
			'n_TUs_with_2plus_subgen': int(multi.size),
			'summary': {name: {'count': int(c), 'percent': round(p, 2)}
				for name, c, p in summary},
			'git': git_info(),
			}, f, indent=2)
	print('\n%d TUs contain >= 2 subgen genes.' % multi.size)
	print('Wrote 4 files to %s' % args.output_dir)
	if missing:
		print('WARNING: %d subgen cistrons not found in sim_data' % len(missing))


if __name__ == '__main__':
	main()
