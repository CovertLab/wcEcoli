"""
Definition-5 rewrite of transcriptFrequency.py.

The original was disabled (early `return`) and operated on legacy TU-level
BulkMolecules counts. This variant produces the intended frequency figure using
Definition 5: for every protein-coding gene it plots the Definition-4 expression
frequency (fraction of successful cell-generations with >= 1 completed
transcript) against the gene's rank by mean completed transcripts per
generation, colored by the canonical def5_CI subgen category. A companion
histogram shows the distribution of the per-gene Def-5 mean.

Reads the pre-computed raw extraction (run subgen_raw_extract.py first).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.cohort import subgen_common as sc
from wholecell.analysis.analysis_tools import exportFigure


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		clf = sc.canonical_def5_classification(plotOutDir)
		stats = clf['stats']
		if clf['n_lineages'] == 0:
			print('No successful lineages found. Skipping.')
			return

		mean = stats['mean']
		cat = stats['cat']
		p_any_synth = clf['p_any_synth']

		# Rank genes by mean completed transcripts per generation (ascending).
		order = np.argsort(mean)
		rank = np.arange(len(mean))

		fig, (ax_scatter, ax_hist) = plt.subplots(
			1, 2, figsize=(13, 5), gridspec_kw={'width_ratios': [2, 1]})

		# --- Scatter: expression frequency vs rank, colored by category ---
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

		# --- Histogram of the per-gene Def-5 mean (log-x for the positive tail) ---
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

		counts = sc.category_counts(cat)
		fig.suptitle('%d successful lineages | subgen=%d, possibly=%d, not=%d, '
			'never=%d' % (clf['n_lineages'], counts['subgen'],
			counts['possibly_subgen'], counts['not_subgen'],
			counts['never_expressed']), fontsize=10)
		plt.tight_layout(rect=[0, 0, 1, 0.96])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
