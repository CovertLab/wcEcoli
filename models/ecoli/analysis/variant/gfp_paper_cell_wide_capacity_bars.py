"""
Cell-wide RNAP and growth capacity, high GFP over no GFP.

# Figure 5H

Four bars, each the mean over surviving high-GFP lineages of the lineage's
window mean over the same seed's no-GFP lineage, with a t-based 95% CI:
	RNAP formed           RNAP complexes assembled per second
	                      (ComplexationListener complexationEvents, summed
	                      over the reactions that produce full RNAP)
	total RNAP            active RNAP (unique molecules) + free RNAP (bulk)
	transcription budget  all RNAP initiations per second (RnapData
	                      didInitialize)
	growth rate           Mass.instantaneous_growth_rate, a specific rate
	                      (per unit mass), which is what the storyboard's
	                      "mass accumulation rate" bar showed

While decision D4 is open a second version (suffix DOUBLING_TIME_PART)
replaces the last bar with doubling time, high GFP over no GFP. That ratio
reads above 1 where the others read below it, because a longer doubling time
is the same slowdown as a lower growth rate.

The supporting claim: transcription and translation capacity fall further
than gene copies, so machinery proteins drop more than dosage alone predicts.
"""

from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

CAPACITY_BARS = (
	('RNAP\nformed', 'rnap_formed_per_s', style.COLORS['rnap']),
	('total\nRNAP', 'total_rnap', style.COLORS['rnap']),
	('transcription\nbudget', 'rnap_initiations_per_s', style.COLORS['rnap']),
)
# (file suffix, last bar): one version per growth metric
VERSIONS = (
	(None, ('growth\nrate', 'growth_rate_per_h', style.INK_MUTED)),
	(C.DOUBLING_TIME_PART, ('doubling\ntime', 'doubling_time_min',
		style.INK_MUTED)),
)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		high = gfp_paper_data.nearest_design(data.designs,
			C.FIG5_PROTEOME_SHARES[-1])
		for part, last in VERSIONS:
			bars = CAPACITY_BARS + (last,)
			ratios = [gfp_paper_data.seed_matched_ratio(data.lineages,
				int(high['variant']), column) for _, column, _ in bars]

			style.apply()
			fig, ax = plt.subplots(figsize=style.panel_size('single', 0.7))
			style.ratio_bars(ax, [label for label, _, _ in bars],
				[r['mean'] for r in ratios], [r['ci95'] for r in ratios],
				[color for _, _, color in bars])
			ax.set_ylabel('high GFP / no GFP')

			style.save_panel(fig, plotOutDir,
				style.panel_name(__file__, part), script=__file__,
				batch=inputDir, sim_data_path=simDataFile,
				n=min(r['n'] for r in ratios),
				notes='High GFP = variant %d (promoter %g, RBS %g, GFP %.3f of'
					' proteome). Seed-matched ratios over %d lineages: %s.'
					% (high['variant'], high['promoter'], high['rbs'],
					high['gfp_proteome_share'], min(r['n'] for r in ratios),
					', '.join('%s %.3f +/- %.3f' % (label.replace('\n', ' '),
					r['mean'], r['ci95']) for (label, _, _), r in zip(bars,
					ratios))),
				**data.provenance([gfp_paper_data.CONTROL_VARIANT,
				int(high['variant'])]))


if __name__ == '__main__':
	Plot().cli()
