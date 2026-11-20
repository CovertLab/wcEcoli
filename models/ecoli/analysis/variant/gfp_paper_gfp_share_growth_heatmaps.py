"""
GFP share of the proteome and growth across the promoter x RBS grid.

# Figure 2B

Two heatmaps over the induced designs, each cell a window mean over the
design's surviving lineages; hatched cells are designs where some lineages
did not complete, and a hatched cell with a dash is a design that lost every
lineage. Left: GFP share of proteome (%), on the GFP ramp. Right:
growth, in two versions while decision D4 is open:
	growth rate     relative to the seed-matched no-GFP control (plain name)
	doubling time   minutes (suffix DOUBLING_TIME_PART)
Growth uses a neutral ramp, darker for slower growth, because blue and green
already stand for RNAP and GFP.

The supporting claim: promoter and RBS strength set GFP output along
different paths, and the highest-producing designs grow slowest.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0


def share_label(value):
	"""GFP share in percent: one decimal below 1%, whole numbers above."""
	return '%.1f' % value if value < 1 else '%.0f' % value


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs.copy()
		designs['gfp_share_pct'] = PERCENT * designs['gfp_proteome_share']
		induced = designs[designs['promoter'] > 0]
		control_dt = float(designs.loc[designs['variant']
			== gfp_paper_data.CONTROL_VARIANT, 'doubling_time_min'].iloc[0])

		versions = (
			(None, 'relative_growth_rate', 'growth rate relative to no GFP',
				(min(0.2, np.nanmin(induced['relative_growth_rate'])), 1.1),
				True, lambda v: '%.2f' % v),
			(C.DOUBLING_TIME_PART, 'doubling_time_min', 'doubling time (min)',
				(min(control_dt, np.nanmin(induced['doubling_time_min'])),
				max(150.0, np.nanmax(induced['doubling_time_min']))),
				False, lambda v: '%.0f' % v),
		)
		for part, column, label, (lo, hi), reverse, fmt in versions:
			style.apply()
			fig, (left, right) = plt.subplots(1, 2,
				figsize=style.panel_size('double', 0.42))
			_, _, hatched = style.design_heatmap(left, designs,
				'gfp_share_pct',
				style.ramp_cmap(style.GFP_RAMP), 0,
				max(75.0, float(np.nanmax(induced['gfp_share_pct']))),
				fmt=share_label)
			left.set_title('GFP share of proteome (%)')
			style.design_heatmap(right, designs, column,
				style.ramp_cmap(style.GROWTH_RAMP, reverse=reverse), lo, hi,
				fmt=fmt)
			right.set_title(label)
			if hatched:
				style.legend_below(fig, [style.hatched_legend_handle()],
					ncol=1)

			variants = sorted(int(v) for v in designs['variant'])
			style.save_panel(fig, plotOutDir,
				style.panel_name(__file__, part), script=__file__,
				batch=inputDir, sim_data_path=simDataFile,
				n=int(len(induced)),
				notes='Window means over surviving lineages. Right: %s.'
					' Designs with no surviving lineage are blank: %s.'
					% (label, gfp_paper_data.lost_designs(induced)),
				**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
