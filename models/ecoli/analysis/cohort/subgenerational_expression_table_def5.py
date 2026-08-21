"""
Classifies genes as subgenerational using Definition 5 (mean completed mRNA
transcripts per generation), with the confidence intervals: a gene is `subgen`
if its 95% confidence interval (estimated across successful lineages) is
entirely below 1 transcript per generation. 

This is the per-gene table. `sc.load_def5_categories()` can read
this file or `subgen_definition5_lineage_ci_pergene_successful.tsv` interchangeably.
It also carries `p_mrna_present_def1`, which subsumes the legacy
subgenerational_expression_table.py `p_expressed` column.

Reads the pre-computed raw extraction (run subgen_raw_extract.py first). Writes:
  <name>.tsv          all protein-coding genes with their Def-5 classification
                      and descriptive copy numbers.
  <name>_subgen.tsv   the subset classified `subgen` (def5_CI).
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

		# The gene key and the synth matrix are written from the same list,
		gene_key, cistron_ids, monomer_ids = sc.load_raw_genes(plotOutDir)
		if gene_key != gene_ids:
			print('WARNING: gene column order differs between the synth matrix '
				'and the gene key; realigning the gene key to the synth-matrix '
				'order.')
			cistron_by_gene = dict(zip(gene_key, cistron_ids))
			monomer_by_gene = dict(zip(gene_key, monomer_ids))
			cistron_ids = [cistron_by_gene.get(g, '') for g in gene_ids]
			monomer_ids = [monomer_by_gene.get(g, '') for g in gene_ids]

		# Descriptive copy numbers: max over successful cells only.
		is_succ = clf['is_successful']
		_, _, _, mrna_gene_ids, max_mrna = sc.load_raw_max(plotOutDir, 'mrna')
		_, _, _, prot_gene_ids, max_prot = sc.load_raw_max(plotOutDir, 'protein')
		for which, ids in (('mrna', mrna_gene_ids), ('protein', prot_gene_ids)):
			if ids != gene_ids:
				print('WARNING: the max-%s matrix gene columns do not match the '
					'synth matrix; its max counts may be misattributed. Re-run '
					'subgen_raw_extract.py.' % which)
		max_mRNA_count = max_mrna[is_succ].max(axis=0) if is_succ.any() \
			else np.zeros(n_genes)
		max_protein_count = max_prot[is_succ].max(axis=0) if is_succ.any() \
			else np.zeros(n_genes)
		# Definition 1 (mRNA presence) recovered without touching simOut: "at least
		# one mRNA copy at some timestep this generation" 
		p_mrna_present = (max_mrna[is_succ] > 0).mean(axis=0) if is_succ.any() \
			else np.full(n_genes, np.nan)

		cat = stats['cat']
		# Column set is the 3636
		#  UNION of the two  per-gene tables, so this file
		# can be read by anything that reads either. The first four columns are
		#  `gene_name`/`cistron_name`/`protein_name` and `gene_id`/`cistron_id`/`protein_id`
		# plus `category` (subgen_definition5_lineage_ci.py's names, which the
		# standalone CLIs read). sc.load_def5_categories() accepts either.
		# monomer_name is the untrimmed protein id (e.g. FOO-MONOMER[c]);
		# protein_name is the same with the [c]/[m]... location suffix removed.
		columns = ['gene_name', 'cistron_name', 'protein_name', 'monomer_name',
			'def5_CI_category', 'def5_CI_mean', 'ci_lower', 'ci_upper',
			'def5_pooled_mean', 'p_any_synth_def4', 'p_mrna_present_def1',
			'max_mRNA_count', 'max_protein_count',
			'n_lineages', 'std', 'se',
			'gene_id', 'cistron_id', 'protein_id', 'category']

		def row(i):
			return [
				gene_ids[i], cistron_ids[i], monomer_ids[i][:-3], monomer_ids[i],
				cat[i], '%.6g' % stats['mean'][i], '%.6g' % stats['ci_low'][i],
				'%.6g' % stats['ci_high'][i], '%.6g' % clf['def5'][i],
				'%.6g' % clf['p_any_synth'][i], '%.6g' % p_mrna_present[i],
				int(max_mRNA_count[i]), int(max_protein_count[i]),
				stats['n'], '%.6g' % stats['std'][i], '%.6g' % stats['se'][i],
				gene_ids[i], cistron_ids[i], monomer_ids[i][:-3], cat[i]]

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
