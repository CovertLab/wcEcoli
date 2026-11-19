"""
Growth against GFP's share of active ribosomes, one point per cell cycle.

# Figure 3H

Every divided cell cycle in the analysis window of every lineage, so the
scatter shows cycle-to-cycle variation within a design as well as the
difference between designs. Cycles from lineages that later failed are drawn
hollow. Two versions are written, since decision D4 is open:
	growth rate     the cycle's mean growth rate over the seed-matched no-GFP
	                lineage's window mean; dashed line at 1. Cycles whose seed
	                has no control lineage are left out.
	doubling time   the cycle's length in minutes, on a log axis as in the
	                storyboard; dashed line at the control's mean. Suffix
	                DOUBLING_TIME_PART.

Designs are coloured along the GFP ramp in order of their mean share of
active ribosomes, so colour darkens with burden, and the control is grey. The
Spearman correlation over induced cycles is recorded in each sidecar, not
drawn.

The supporting claim: the larger the share of the ribosome pool a cell gives
to gfp, the slower it grows.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style


DOUBLING_TIME_TICKS = (30, 50, 75, 100, 150, 200, 300)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		cycles = data.window_cycles()
		cycles = cycles[np.isfinite(cycles['gfp_ribosome_share'])]
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
		order = [gfp_paper_data.CONTROL_VARIANT] + induced

		control = cycles[cycles['variant'] == gfp_paper_data.CONTROL_VARIANT]
		versions = (
			(None, 'relative_growth_rate', 'growth rate relative to no GFP',
				1.0, False),
			(C.DOUBLING_TIME_PART, 'doubling_time_min', 'doubling time (min)',
				float(control['doubling_time_min'].mean()), True),
		)
		for part, column, ylabel, reference, log in versions:
			self._version(cycles, order, colors, part, column, ylabel,
				reference, log, data, inputDir, plotOutDir, simDataFile)

	@staticmethod
	def _version(cycles, order, colors, part, column, ylabel, reference, log,
			data, inputDir, plotOutDir, simDataFile):
		"""Draw and save one growth metric's version of the panel."""
		from scipy import stats
		cycles = cycles[np.isfinite(cycles[column])]
		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.75))
		for variant in order:
			group = cycles[cycles['variant'] == variant]
			if group.empty:
				continue
			promoter, rbs = group['promoter'].iloc[0], group['rbs'].iloc[0]
			label = 'no GFP' if variant == gfp_paper_data.CONTROL_VARIANT \
				else '%s promoter, %s' % (C.promoter_label(promoter),
				C.rbs_label(rbs))
			style.design_scatter(ax, group['gfp_ribosome_share'],
				group[column], colors[variant],
				failed=group['lineage_failed'], size=9, alpha=0.8,
				label=label, zorder=2)
		if np.isfinite(reference):
			ax.axhline(reference, color=style.INK_MUTED,
				linestyle=(0, (4, 3)), linewidth=0.8, zorder=1)
		ax.set_xlim(left=0)
		if log:
			style.log_axis(ax, 'y', 30, 150, data=cycles[column].values,
				ticks=DOUBLING_TIME_TICKS)
		else:
			ax.set_ylim(bottom=0)
		ax.set_xlabel('GFP share of active ribosomes')
		ax.set_ylabel(ylabel)
		style.legend_with_failed(ax, bool(cycles['lineage_failed'].any()),
			loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)

		on = cycles[cycles['variant'] != gfp_paper_data.CONTROL_VARIANT]
		if len(on) > 2:
			rho, p = stats.spearmanr(on['gfp_ribosome_share'], on[column])
			correlation = 'Spearman rho %.3f (p %.2g) over %d induced cycles' \
				% (rho, p, len(on))
		else:
			correlation = 'too few induced cycles for a correlation'

		style.save_panel(fig, plotOutDir, style.panel_name(__file__, part),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(len(cycles)),
			notes='%s. One point per divided cycle in generations %d-%d.'
				' Reference line %.3f. %s.' % (ylabel, C.ANALYSIS_FIRST_GEN,
				C.ANALYSIS_LAST_GEN, reference, correlation),
			**data.provenance())


if __name__ == '__main__':
	Plot().cli()
