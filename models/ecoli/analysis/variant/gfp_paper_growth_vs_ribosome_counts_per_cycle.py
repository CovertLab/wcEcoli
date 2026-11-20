"""
Growth against ribosomes per cell, one point per cell cycle.

# Figure 2F

Every divided cell cycle in the analysis window, for the no-GFP control and
the designs nearest each FIG5_PROTEOME_SHARES, so the same "GFP 14%" and
"GFP 35%" as Figure 5. Ribosomes are the cycle's mean total (active plus
free 70S-forming pairs). Cycles from lineages that later failed are hollow.
Two versions while decision D4 is open:
	growth rate     cycle growth rate over the seed-matched control lineage's
	                window mean; dashed line at 1 (plain name)
	doubling time   cycle length in minutes, log axis; dashed line at the
	                control mean (suffix DOUBLING_TIME_PART)

The supporting claim: fewer ribosomes means slower growth, along one curve
whatever the GFP level, so ribosome count rather than GFP itself sets growth.
The Spearman correlation over all plotted cycles is in the sidecar.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
DOUBLING_TIME_TICKS = (30, 50, 75, 100, 150, 200, 300)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		control = gfp_paper_data.CONTROL_VARIANT
		chosen = [gfp_paper_data.nearest_design(data.designs, share)
			for share in C.FIG5_PROTEOME_SHARES]
		variants = list(dict.fromkeys([control] + [int(d['variant'])
			for d in chosen]))
		labels = {control: 'no GFP'}
		labels.update({int(d['variant']): 'GFP %.0f%% of proteome'
			% (PERCENT * d['gfp_proteome_share']) for d in chosen})
		steps = style.GFP_RAMP[1:4]
		colors = {control: style.COLORS['control']}
		colors.update(dict(zip(variants[1:], [steps[int(round(i))] for i in
			np.linspace(0, len(steps) - 1, max(len(variants) - 1, 1))])))

		cycles = data.window_cycles()
		cycles = cycles[cycles['variant'].isin(variants)
			& np.isfinite(cycles['total_ribosomes'])]
		survived = data.lineages.set_index(['variant', 'seed'])['survived']
		cycles = cycles.assign(lineage_failed=[not survived[(v, s)]
			for v, s in zip(cycles['variant'], cycles['seed'])])
		control_dt = float(cycles[cycles['variant'] == control][
			'doubling_time_min'].mean())

		for part, column, ylabel, reference, log in (
				(None, 'relative_growth_rate',
					'growth rate relative to no GFP', 1.0, False),
				(C.DOUBLING_TIME_PART, 'doubling_time_min',
					'doubling time (min)', control_dt, True)):
			self._version(cycles, variants, labels, colors, part, column,
				ylabel, reference, log, data, inputDir, plotOutDir,
				simDataFile)

	@staticmethod
	def _version(cycles, variants, labels, colors, part, column, ylabel,
			reference, log, data, inputDir, plotOutDir, simDataFile):
		"""Draw and save one growth metric's version of the panel."""
		from scipy import stats
		cycles = cycles[np.isfinite(cycles[column])]
		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.8))
		for variant in variants:
			group = cycles[cycles['variant'] == variant]
			if group.empty:
				continue
			style.design_scatter(ax, group['total_ribosomes'], group[column],
				colors[variant], failed=group['lineage_failed'], size=9,
				alpha=0.7, label=labels[variant], zorder=2)
		if np.isfinite(reference):
			ax.axhline(reference, color=style.INK_MUTED,
				linestyle=(0, (4, 3)), linewidth=0.8, zorder=1)
		ax.set_xlim(0, max(30000.0, float(cycles['total_ribosomes'].max())
			* 1.05) if len(cycles) else 30000.0)
		ax.xaxis.set_major_formatter(plt.FuncFormatter(
			lambda v, _: '{:,.0f}'.format(v)))
		if log:
			style.log_axis(ax, 'y', 30, 150, data=cycles[column].values,
				ticks=DOUBLING_TIME_TICKS)
		else:
			ax.set_ylim(bottom=0)
		ax.set_xlabel('ribosomes per cell')
		ax.set_ylabel(ylabel)
		style.legend_with_failed(ax, bool(cycles['lineage_failed'].any()),
			below=True)
		if len(cycles) > 2:
			rho, p = stats.spearmanr(cycles['total_ribosomes'], cycles[column])
			correlation = 'Spearman rho %.3f (p %.2g) over %d cycles' % (
				rho, p, len(cycles))
		else:
			correlation = 'too few cycles for a correlation'
		style.save_panel(fig, plotOutDir, style.panel_name(__file__, part),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(len(cycles)),
			notes='%s against total ribosomes, one point per divided cycle in'
				' generations %d-%d. Designs: %s. %s.' % (ylabel,
				C.ANALYSIS_FIRST_GEN, C.ANALYSIS_LAST_GEN, ', '.join(
				'%d = %s' % (v, labels[v]) for v in variants), correlation),
			**data.provenance(variants))


if __name__ == '__main__':
	Plot().cli()
