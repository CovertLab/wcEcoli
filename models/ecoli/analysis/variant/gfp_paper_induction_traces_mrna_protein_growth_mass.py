"""
gfp mRNA, GFP protein, growth and cell mass through induction, every lineage.

# Figure 1B

Four tracks against generation (growth rate smoothed as every paper trace
is, see gfp_paper_data.GROWTH_TRACE_SMOOTH_STEPS) for the showcase design (nearest
SHOWCASE_PROTEOME_SHARE): one representative lineage in dark green, the other
lineages in light green, and the no-GFP control as its median and 10th-90th
percentile band across lineages (growth and mass only; the control has no
gfp). Each cycle is split into N_POINTS equal fractions and averaged within
each, so lineages line up generation by generation. The dashed line is
induction.

The representative lineage is the surviving lineage whose window-mean growth,
relative to its seed-matched control, is closest to the median of all
surviving lineages: typical rather than chosen by eye.

Divergence bar (fixed before the data, see gfp_paper_style.paired_divergence):
per generation, each seed's cycle growth with GFP minus the same seed's
control, tested against zero across seeds (paired t-test), Holm-corrected
across the generations. The bar starts at the first generation that is
significant and stays significant to the end; the adjusted p-values are in
the sidecar.

Two versions while decision D4 is open:
	growth rate     instantaneous growth rate (1/h) track (plain name)
	doubling time   per-cycle doubling time (min) as a step per generation,
	                and the divergence test run on doubling time (suffix
	                DOUBLING_TIME_PART)
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

N_POINTS = 40
SERIES = ('gfp_mrna', 'gfp_protein', 'growth_rate', 'cell_mass')
REPRESENTATIVE = style.GFP_RAMP[3]
OTHERS = style.GFP_RAMP[1]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		showcase = gfp_paper_data.nearest_design(data.designs,
			C.SHOWCASE_PROTEOME_SHARE)
		design = int(showcase['variant'])
		control = gfp_paper_data.CONTROL_VARIANT
		traces = gfp_paper_data.cycle_traces(data, [control, design],
			plotOutDir, self.cpus, n_points=N_POINTS, series=SERIES)
		lineages = data.lineages[data.lineages['variant'] == design]
		alive = lineages[lineages['survived']
			& np.isfinite(lineages['relative_growth_rate'])]
		if alive.empty:
			alive = lineages[lineages['survived']]
		median = alive['relative_growth_rate'].median()
		representative = int(alive.loc[(alive['relative_growth_rate']
			- median).abs().idxmin(), 'seed']) if len(alive) else None

		for part, column, ylabel in (
				(None, 'growth_rate_per_h', 'growth rate\n(1/h)'),
				(C.DOUBLING_TIME_PART, 'doubling_time_min',
					'doubling time\n(min)')):
			test = self._divergence(data, design, column)
			self._version(data, traces, design, representative, part,
				column, ylabel, test, showcase, inputDir, plotOutDir,
				simDataFile)

	@staticmethod
	def _divergence(data, design, column):
		"""Paired per-generation test of design against control, by seed."""
		cycles = data.cycles[data.cycles['divided']]
		control = cycles[cycles['variant'] == gfp_paper_data.CONTROL_VARIANT
			].pivot_table(index='seed', columns='generation', values=column)
		treated = cycles[cycles['variant'] == design].pivot_table(
			index='seed', columns='generation', values=column)
		seeds = sorted(set(control.index) & set(treated.index))
		generations = list(range(C.N_GENERATIONS))
		control = control.reindex(index=seeds, columns=generations)
		treated = treated.reindex(index=seeds, columns=generations)
		if len(seeds) < 2:
			return None
		return style.paired_divergence(treated.values, control.values)

	@staticmethod
	def _version(data, traces, design, representative, part, column, ylabel,
			test, showcase, inputDir, plotOutDir, simDataFile):
		control_id = gfp_paper_data.CONTROL_VARIANT
		style.apply()
		fig, axes = plt.subplots(4, 1, sharex=True, figsize=style.panel_size(
			'double', 0.5))
		tracks = (('gfp_mrna', 'GFP mRNA\n(copies)', False),
			('gfp_protein', 'GFP protein\n(copies)', False),
			('growth', ylabel, True),
			('cell_mass', 'cell mass\n(fg)', True))
		design_traces = traces[traces['variant'] == design]
		control_traces = traces[traces['variant'] == control_id]
		for ax, (name, label, band) in zip(axes, tracks):
			if name == 'growth' and part is not None:
				Plot._doubling_steps(ax, data, design, representative)
			else:
				series = 'growth_rate' if name == 'growth' else name
				if band:
					stats = control_traces.groupby('x')[series].quantile(
						[0.1, 0.5, 0.9]).unstack()
					ax.fill_between(stats.index, stats[0.1], stats[0.9],
						color=style.COLORS['control'], alpha=0.25,
						linewidth=0, zorder=1)
					ax.plot(stats.index, stats[0.5],
						color=style.COLORS['control'], linewidth=0.9,
						zorder=2)
				for seed, group in design_traces.groupby('seed'):
					is_rep = seed == representative
					ax.plot(group['x'], group[series], color=REPRESENTATIVE
						if is_rep else OTHERS, linewidth=1.0 if is_rep else 0.5,
						alpha=1.0 if is_rep else 0.5, zorder=4 if is_rep
						else 3)
			ax.axvline(C.NEW_GENE_INDUCTION_GEN, color=style.INK_MUTED,
				linestyle=(0, (4, 3)), linewidth=0.8, zorder=0)
			ax.set_ylim(bottom=0)
			ax.set_ylabel(label)
			if name == 'gfp_protein':
				ax.yaxis.set_major_formatter(plt.FuncFormatter(
					lambda v, _: '%gM' % (v / 1e6) if v >= 1e6 else
					'{:,.0f}'.format(v)))
			elif name in ('gfp_mrna', 'cell_mass'):
				ax.yaxis.set_major_formatter(plt.FuncFormatter(
					lambda v, _: '{:,.0f}'.format(v)))
			if name == 'growth' and test is not None and (
					test['first_sustained'] is not None):
				start = test['first_sustained']
				top = ax.get_ylim()[1]
				ax.plot([start, C.N_GENERATIONS], [top, top], color=style.INK,
					linewidth=2.0, solid_capstyle='butt', clip_on=False,
					zorder=5)
				worst = float(np.max(test['p_adjusted'][start:]))
				ax.annotate('p < %s from generation %d (paired, Holm)'
					% (_p_text(worst), start), xy=(C.N_GENERATIONS, top),
					xytext=(0, 3), textcoords='offset points', ha='right',
					va='bottom', fontsize=style.FONT_SIZES['annotation'],
					color=style.INK)
		axes[0].annotate('induction', xy=(C.NEW_GENE_INDUCTION_GEN, 1),
			xycoords=('data', 'axes fraction'), xytext=(3, -2),
			textcoords='offset points', va='top',
			fontsize=style.FONT_SIZES['annotation'], color=style.INK_MUTED)
		axes[-1].set_xlim(0, C.N_GENERATIONS)
		axes[-1].set_xticks(range(0, C.N_GENERATIONS + 1, 4))
		axes[-1].set_xlabel('generation')
		fig.legend(handles=[
			Line2D([], [], color=REPRESENTATIVE, label='GFP, one lineage'),
			Line2D([], [], color=OTHERS, label='other GFP lineages'),
			Patch(facecolor=style.COLORS['control'], alpha=0.4,
				label='no GFP, median and 10th-90th percentile')],
			loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.06))

		if test is None:
			test_note = 'Divergence test not run: fewer than 2 seed pairs.'
		else:
			test_note = ('Divergence on %s: first sustained generation %s;'
				' Holm-adjusted p per generation: %s.' % (column,
				test['first_sustained'], ', '.join('%d:%.2g' % (g, p)
				for g, p in enumerate(test['p_adjusted']))))
		style.save_panel(fig, plotOutDir, style.panel_name(__file__, part),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(design_traces['seed'].nunique()),
			notes='Showcase = variant %d (GFP %.3f of proteome);'
				' representative lineage seed %s. %d within-cycle bins. %s'
				% (design, showcase['gfp_proteome_share'], representative,
				N_POINTS, test_note),
			**data.provenance([control_id, design]))

	@staticmethod
	def _doubling_steps(ax, data, design, representative):
		"""Per-cycle doubling time, a flat step across each generation."""
		cycles = data.cycles[data.cycles['divided']]
		control = cycles[cycles['variant'] == gfp_paper_data.CONTROL_VARIANT]
		stats = control.groupby('generation')['doubling_time_min'].quantile(
			[0.1, 0.5, 0.9]).unstack()
		edges = np.append(stats.index.values, stats.index.values[-1] + 1)
		ax.fill_between(edges, np.append(stats[0.1], stats[0.1].iloc[-1]),
			np.append(stats[0.9], stats[0.9].iloc[-1]), step='post',
			color=style.COLORS['control'], alpha=0.25, linewidth=0, zorder=1)
		ax.step(edges, np.append(stats[0.5], stats[0.5].iloc[-1]),
			where='post', color=style.COLORS['control'], linewidth=0.9,
			zorder=2)
		for seed, group in cycles[cycles['variant'] == design].groupby('seed'):
			group = group.sort_values('generation')
			is_rep = seed == representative
			x = np.append(group['generation'].values,
				group['generation'].values[-1] + 1)
			y = np.append(group['doubling_time_min'].values,
				group['doubling_time_min'].values[-1])
			ax.step(x, y, where='post', color=REPRESENTATIVE if is_rep
				else OTHERS, linewidth=1.0 if is_rep else 0.5,
				alpha=1.0 if is_rep else 0.5, zorder=4 if is_rep else 3)


def _p_text(p):
	"""A p-value threshold for the bar label: the next round number above."""
	for threshold in (0.001, 0.01, 0.05):
		if p < threshold:
			return '%g' % threshold
	return '%.2g' % p


if __name__ == '__main__':
	Plot().cli()
