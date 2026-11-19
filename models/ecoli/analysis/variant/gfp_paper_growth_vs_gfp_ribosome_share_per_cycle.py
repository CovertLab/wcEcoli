"""
Growth rate against GFP's share of active ribosomes, one point per cell
cycle.

# Figure 3H

Every divided cell cycle in the analysis window of every lineage, so the
scatter shows cycle-to-cycle variation within a design as well as the
difference between designs. Growth is the cycle's mean growth rate over the
seed-matched no-GFP control lineage's window mean (decision D4, growth rate
rather than doubling time). Cycles from lineages that later failed are drawn
hollow. The dashed line is the control (1).

Designs are coloured along the GFP ramp in order of their mean share of
active ribosomes, so colour darkens with burden, and the control is grey. The Spearman correlation over induced cycles is
recorded in the sidecar, not drawn.

The supporting claim: the larger the share of the ribosome pool a cell gives
to gfp, the slower it grows.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		from scipy import stats
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		cycles = data.window_cycles()
		cycles = cycles[np.isfinite(cycles['relative_growth_rate'])
			& np.isfinite(cycles['gfp_ribosome_share'])]
		survived = data.lineages.set_index(['variant', 'seed'])['survived']
		cycles = cycles.assign(lineage_failed=[not survived[(v, s)]
			for v, s in zip(cycles['variant'], cycles['seed'])])

		# Order designs by how much of the ribosome pool they take, so the
		# ramp darkens with burden and the legend reads in the same order.
		share = data.designs.set_index('variant')['gfp_ribosome_share']
		induced = sorted((v for v in cycles['variant'].unique()
			if v != gfp_paper_data.CONTROL_VARIANT),
			key=lambda v: (np.nan_to_num(share.get(v, np.nan), nan=np.inf), v))
		ramp = style.gfp_ramp(len(induced) + 1)[1:]
		colors = {gfp_paper_data.CONTROL_VARIANT: style.COLORS['control']}
		colors.update(dict(zip(induced, ramp)))

		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.75))
		for variant in [gfp_paper_data.CONTROL_VARIANT] + induced:
			group = cycles[cycles['variant'] == variant]
			if group.empty:
				continue
			promoter, rbs = group['promoter'].iloc[0], group['rbs'].iloc[0]
			label = 'no GFP' if variant == gfp_paper_data.CONTROL_VARIANT \
				else '%s promoter, %s' % (C.promoter_label(promoter),
				C.rbs_label(rbs))
			style.design_scatter(ax, group['gfp_ribosome_share'],
				group['relative_growth_rate'], colors[variant],
				failed=group['lineage_failed'], size=9, alpha=0.8,
				label=label, zorder=2)
		ax.axhline(1.0, color=style.INK_MUTED, linestyle=(0, (4, 3)),
			linewidth=0.8, zorder=1)
		ax.set_xlim(left=0)
		ax.set_ylim(bottom=0)
		ax.set_xlabel('GFP share of active ribosomes')
		ax.set_ylabel('growth rate relative to no GFP')
		style.legend_with_failed(ax, bool(cycles['lineage_failed'].any()),
			loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		on = cycles[cycles['variant'] != gfp_paper_data.CONTROL_VARIANT]
		if len(on) > 2:
			rho, p = stats.spearmanr(on['gfp_ribosome_share'],
				on['relative_growth_rate'])
			correlation = 'Spearman rho %.3f (p %.2g) over %d induced cycles' \
				% (rho, p, len(on))
		else:
			correlation = 'too few induced cycles for a correlation'

		style.save_panel(fig, plotOutDir, style.panel_name(__file__),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(len(cycles)),
			notes='One point per divided cycle in generations %d-%d. Cycles'
				' whose seed has no control lineage are omitted. %s.'
				% (C.ANALYSIS_FIRST_GEN, C.ANALYSIS_LAST_GEN, correlation),
			**data.provenance())


if __name__ == '__main__':
	Plot().cli()
