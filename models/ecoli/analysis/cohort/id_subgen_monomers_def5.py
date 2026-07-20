"""
Definition-5 rewrite of id_subgen_monomers.py.

The original classified genes by PROTEIN existence (near-degenerate: proteins
accumulate and rarely return to zero, so almost everything reads as "always
expressed") and, despite its name, never wrote the classification to disk. This
variant classifies each protein-coding gene by Definition 5 (completed mRNA
transcripts per generation) using the canonical Form B rule (subgen iff the 95%
CI is entirely below 1 transcript/gen across successful lineages) and WRITES the
resulting category per monomer.

Reads the pre-computed raw extraction (run subgen_raw_extract.py first). Writes:
  <name>.tsv  gene_name, cistron_name, monomer_name, def5_category, def5_mean,
              ci_lower, ci_upper, prob_any_synth (all genes/monomers).
"""

import csv
import os

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		clf = sc.canonical_def5_classification(plotOutDir)
		gene_ids = clf['gene_ids']
		stats = clf['stats']
		n_genes = len(gene_ids)
		print('Classified %d monomers over %d successful lineages.'
			% (n_genes, clf['n_lineages']))
		if clf['n_lineages'] == 0:
			print('No successful lineages found. Skipping.')
			return

		_, cistron_ids, monomer_ids = sc.load_raw_genes(plotOutDir)
		cat = stats['cat']

		out_path = os.path.join(plotOutDir, plotOutFileName + '.tsv')
		print('Writing %s' % out_path)
		with open(out_path, 'w') as f:
			w = csv.writer(f, delimiter='\t')
			w.writerow(['gene_name', 'cistron_name', 'monomer_name',
				'def5_category', 'def5_mean', 'ci_lower', 'ci_upper',
				'prob_any_synth'])
			for i in range(n_genes):
				w.writerow([
					gene_ids[i], cistron_ids[i], monomer_ids[i], cat[i],
					'%.6g' % stats['mean'][i], '%.6g' % stats['ci_low'][i],
					'%.6g' % stats['ci_high'][i], '%.6g' % clf['p_any_synth'][i]])

		counts = sc.category_counts(cat)
		print('Definition-5 categories: '
			+ ', '.join('%s=%d' % (c, counts[c]) for c in sc.CATEGORIES))


if __name__ == '__main__':
	Plot().cli()
