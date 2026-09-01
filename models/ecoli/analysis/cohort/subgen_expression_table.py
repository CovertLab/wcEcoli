"""
Classifies genes as subgenerational using Definition 5 (mean completed mRNA
transcripts per generation), with the confidence intervals: a gene is `subgen`
if its 95% confidence interval (estimated across successful seeds) is
entirely below 1 transcript per generation. 

This is the per-gene table. `sc.load_def5_categories()` can read
this file or `subgen_seed_ci_pergene_successful.tsv` interchangeably.
It also carries two presence probabilities recovered from the extraction without
touching simOut. Both are gated on strict-successful seeds. 

Reads the pre-computed raw extraction (run subgen_extract.py first). Writes:
  <name>.tsv            all protein-coding genes with their Def-5 classification
                        and descriptive copy numbers.
  <name>_subgen.tsv     the subset classified `subgen` (def5_CI).
  <name>_frequency.pdf  the expression-frequency figure (formerly its own
                        script, transcriptFrequency_def5.py, which plotted these
                        same columns from this same classification).
"""

import csv
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_helper_functions as sc
from wholecell.analysis.analysis_tools import exportFigure


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		clf = sc.canonical_def5_classification(plotOutDir)
		gene_ids = clf['gene_ids']
		stats = clf['stats']
		n_genes = len(gene_ids)
		print('Classified %d genes over %d successful seeds.'
			% (n_genes, clf['n_seeds']))
		if clf['n_seeds'] == 0:
			print('No successful seeds found. Skipping.')
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
					'subgen_extract.py.' % which)
		max_mRNA_count = max_mrna[is_succ].max(axis=0) if is_succ.any() \
			else np.zeros(n_genes)
		max_protein_count = max_prot[is_succ].max(axis=0) if is_succ.any() \
			else np.zeros(n_genes)
		# Definition 1 (mRNA presence) recovered without touching simOut: "at least
		# one mRNA copy at some timestep this generation" 
		p_mrna_present = (max_mrna[is_succ] > 0).mean(axis=0) if is_succ.any() \
			else np.full(n_genes, np.nan)
		# Monomer-presence frequency, same idea one level down: the fraction of
		# generations in which the protein existed at all. This is the number the
		# legacy id_subgen_monomers.py existed to produce.
		p_monomer_present = (max_prot[is_succ] > 0).mean(axis=0) if is_succ.any() \
			else np.full(n_genes, np.nan)

		cat = stats['cat']
		# Column set is the
		#  UNION of the two  per-gene tables, so this file
		# can be read by anything that reads either. The first four columns are
		#  `gene_name`/`cistron_name`/`protein_name` and `gene_id`/`cistron_id`/`protein_id`
		# plus `category` (subgen_seed_ci.py's names, which the
		# standalone CLIs read). sc.load_def5_categories() accepts either.
		# monomer_name is the untrimmed protein id (e.g. FOO-MONOMER[c]);
		# protein_name is the same with the [c]/[m]... location suffix removed.
		columns = ['gene_name', 'cistron_name', 'protein_name', 'monomer_name',
			'def5_CI_category', 'def5_CI_mean', 'ci_lower', 'ci_upper',
			'def5_pooled_mean', 'p_any_synth_def4', 'p_mrna_present_def1',
			'p_monomer_present',
			'max_mRNA_count', 'max_protein_count',
			'n_seeds', 'std', 'se',
			'gene_id', 'cistron_id', 'protein_id', 'category']

		def row(i):
			return [
				gene_ids[i], cistron_ids[i], monomer_ids[i][:-3], monomer_ids[i],
				cat[i], '%.6g' % stats['mean'][i], '%.6g' % stats['ci_low'][i],
				'%.6g' % stats['ci_high'][i], '%.6g' % clf['def5'][i],
				'%.6g' % clf['p_any_synth'][i], '%.6g' % p_mrna_present[i],
				'%.6g' % p_monomer_present[i],
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

		self._plot_frequency(plotOutDir, plotOutFileName, clf, stats, cat,
			counts, metadata)

	def _plot_frequency(self, plotOutDir, plotOutFileName, clf, stats, cat,
			counts, metadata):
		"""The expression-frequency figure, adapted from Nora's
		transcriptFrequency cohort script (now in 'old_subgen_scripts').

		It plots columns this table already computed, from the same
		classification, so it lives here rather than in a script of its own.
		"""
		mean = stats['mean']
		p_any_synth = clf['p_any_synth']

		# Rank genes by mean completed transcripts per generation (ascending).
		order = np.argsort(mean)
		rank = np.arange(len(mean))

		fig, (ax_scatter, ax_hist) = plt.subplots(
			1, 2, figsize=(13, 5), gridspec_kw={'width_ratios': [2, 1]})

		# Scatter: expression frequency vs rank, colored by category
		for category in sc.CATEGORIES:
			sel = cat[order] == category
			if not np.any(sel):
				continue
			ax_scatter.scatter(
				rank[sel], p_any_synth[order][sel], s=6, alpha=0.6,
				color=sc.PALETTE[category], label=sc.CAT_LABEL[category])
		ax_scatter.axhline(1.0, color=sc.MUTED, lw=0.8, ls='--')
		ax_scatter.set_xlabel('Gene rank (by mean completed transcripts / gen)')
		ax_scatter.set_ylabel('P(>=1 completed transcript in a generation) [def 4]')
		ax_scatter.set_title('Subgenerational expression frequency (Definition 5)')
		ax_scatter.set_ylim(-0.02, 1.02)
		handles = [Patch(color=sc.PALETTE[c], label=sc.CAT_LABEL[c])
			for c in sc.CATEGORIES]
		ax_scatter.legend(handles=handles, fontsize=8, loc='lower right',
			frameon=False)

		# Histogram of the per-gene Def-5 mean (log-x for the positive tail)
		pos = mean[mean > 0]
		if pos.size:
			bins = np.logspace(
				np.log10(max(pos.min(), 1e-4)), np.log10(pos.max() + 1e-9), 40)
			ax_hist.hist(pos, bins=bins, color=sc.PALETTE['subgen'], alpha=0.8)
			ax_hist.set_xscale('log')
		ax_hist.axvline(1.0, color=sc.PALETTE['not_subgen'], lw=1.2, ls='--',
			label='1 transcript / gen')
		ax_hist.set_xlabel('Mean completed transcripts / generation')
		ax_hist.set_ylabel('Number of genes')
		ax_hist.set_title('Definition-5 rate distribution')
		ax_hist.legend(fontsize=8, frameon=False)

		fig.suptitle('%d successful seeds | subgen=%d, possibly=%d, not=%d, '
			'never=%d' % (clf['n_seeds'], counts['subgen'],
			counts['possibly_subgen'], counts['not_subgen'],
			counts['never_expressed']), fontsize=10)
		plt.tight_layout(rect=[0, 0, 1, 0.96])
		exportFigure(plt, plotOutDir, plotOutFileName + '_frequency', metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
