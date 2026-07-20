"""
Definition-5 rewrite of subgenerational_expression_table.py.

Classifies genes as subgenerational using Definition 5 (mean completed mRNA
transcripts per generation), with the canonical Form B rule: a gene is `subgen`
iff its 95% confidence interval (estimated across successful lineages) is
entirely below 1 transcript per generation. This replaces the original's
mRNA-presence-frequency definition and fixes its read-failure denominator bug by
consuming the consistent per-cell matrix from the raw extraction.

Reads the pre-computed raw extraction (run subgen_raw_extract.py first). Writes:
  <name>.tsv          all protein-coding genes with their Def-5 classification
                      and descriptive copy numbers.
  <name>_subgen.tsv   the subset classified `subgen` (Form B).
"""

import csv
import os

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		clf = sc.canonical_def5_classification(plotOutDir)
		gene_ids = clf['gene_ids']
		stats = clf['stats']
		n_genes = len(gene_ids)
		print('Classified %d genes over %d successful lineages.'
			% (n_genes, clf['n_lineages']))
		if clf['n_lineages'] == 0:
			print('No successful lineages found. Skipping.')
			return

		gene_key, cistron_ids, monomer_ids = sc.load_raw_genes(plotOutDir)
		# Descriptive copy numbers: max over successful cells only.
		is_succ = clf['is_successful']
		_, _, _, _, max_mrna = sc.load_raw_max(plotOutDir, 'mrna')
		_, _, _, _, max_prot = sc.load_raw_max(plotOutDir, 'protein')
		max_mRNA_count = max_mrna[is_succ].max(axis=0) if is_succ.any() \
			else np.zeros(n_genes)
		max_protein_count = max_prot[is_succ].max(axis=0) if is_succ.any() \
			else np.zeros(n_genes)

		cat = stats['cat']
		columns = ['gene_name', 'cistron_name', 'protein_name', 'def5_category',
			'def5_mean', 'ci_lower', 'ci_upper', 'formA_pooled_mean',
			'p_any_synth_def4', 'max_mRNA_count', 'max_protein_count']

		def row(i):
			return [
				gene_ids[i], cistron_ids[i], monomer_ids[i][:-3], cat[i],
				'%.6g' % stats['mean'][i], '%.6g' % stats['ci_low'][i],
				'%.6g' % stats['ci_high'][i], '%.6g' % clf['formA'][i],
				'%.6g' % clf['p_any_synth'][i], int(max_mRNA_count[i]),
				int(max_protein_count[i])]

		main_path = os.path.join(plotOutDir, plotOutFileName + '.tsv')
		print('Writing %s' % main_path)
		with open(main_path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(columns)
			for i in range(n_genes):
				w.writerow(row(i))

		subgen_path = os.path.join(plotOutDir, plotOutFileName + '_subgen.tsv')
		subgen_idx = np.where(cat == 'subgen')[0]
		print('Writing %s (%d subgen genes)' % (subgen_path, len(subgen_idx)))
		with open(subgen_path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(columns)
			for i in subgen_idx:
				w.writerow(row(int(i)))

		counts = sc.category_counts(cat)
		print('Definition-5 categories: '
			+ ', '.join('%s=%d' % (c, counts[c]) for c in sc.CATEGORIES))


if __name__ == '__main__':
	Plot().cli()
