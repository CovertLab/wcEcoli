"""
What changes in the cell with GFP, as small bar charts and as one chart of
relative changes.

# Figure 1E

Small multiples: each quantity for the no-GFP control and the showcase design
(nearest SHOWCASE_PROTEOME_SHARE), means over surviving lineages of each
lineage's window mean, with t-based 95% CIs. The order follows the paper
(GFP itself, then growth and size, then the machinery of Figure 2, the
enzymes of Figure 4 and the replication of Figure 5) and is set by BARS, so
it can be changed in one place. Metabolic enzymes are a concentration: every
FBA catalyst molecule over cell volume, in micromolar.

The relative version ('relative') puts every reduction in one chart: each
quantity in the showcase lineage over the same seed's control, mean with a
t-based 95% CI, dashed line at 1.

Both are written in growth-rate and doubling-time versions while decision D4
is open (suffix DOUBLING_TIME_PART); only the growth bar differs.
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
GROWTH = ('growth rate', '(1/h)', 'growth_rate_per_h', 1.0)
DOUBLING = ('doubling time', '(min)', 'doubling_time_min', 1.0)
# (title, unit line, lineage column, scale); None marks the growth slot
BARS = (
	('GFP mRNA', '(copies)', 'gfp_mrna', 1.0),
	('GFP protein', '(copies)', 'gfp_monomers', 1.0),
	('GFP share of', 'proteome (%)', 'gfp_proteome_share', PERCENT),
	None,
	('cell mass', '(fg)', 'cellMass', 1.0),
	('ribosomes', '(counts)', 'total_ribosomes', 1.0),
	('RNAP', '(counts)', 'total_rnap', 1.0),
	('metabolic enzymes', '(µM)', 'enzyme_uM', 1.0),
	('replication', 'origins (per cell)', 'oric', 1.0),
)
# The relative chart shows the quantities that exist without GFP, with short
# labels so six bars fit one column.
RELATIVE_LABELS = {
	'growth_rate_per_h': 'growth\nrate', 'doubling_time_min': 'doubling\ntime',
	'cellMass': 'cell\nmass', 'total_ribosomes': 'ribo-\nsomes',
	'total_rnap': 'RNAP', 'enzyme_uM': 'enzymes\n(\u00b5M)',
	'oric': 'origins',
}
RELATIVE_COLORS = {
	'growth_rate_per_h': style.INK_MUTED, 'doubling_time_min': style.INK_MUTED,
	'cellMass': style.INK_MUTED, 'total_ribosomes': style.COLORS['ribosome'],
	'total_rnap': style.COLORS['rnap'], 'enzyme_uM': style.COLORS['enzyme'],
	'oric': style.COLORS['replication'],
}


def _mean_ci(values):
	values = np.asarray(values, dtype=float)
	values = values[np.isfinite(values)]
	if values.size == 0:
		return np.nan, np.nan
	return float(values.mean()), gfp_paper_data._ci95(values)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		showcase = gfp_paper_data.nearest_design(data.designs,
			C.SHOWCASE_PROTEOME_SHARE)
		design = int(showcase['variant'])
		control = gfp_paper_data.CONTROL_VARIANT
		alive = data.lineages[data.lineages['survived']]
		label = 'GFP %.0f%% of proteome' % (
			PERCENT * showcase['gfp_proteome_share'])
		colors = {control: style.COLORS['control'], design: style.GFP_RAMP[3]}
		provenance = data.provenance([control, design])

		for growth, part in ((GROWTH, None), (DOUBLING, C.DOUBLING_TIME_PART)):
			bars = [growth if b is None else b for b in BARS]
			self._multiples(bars, alive, control, design, colors, label, part,
				data, inputDir, plotOutDir, simDataFile, provenance)
			self._relative(bars, data, design, label, part, inputDir,
				plotOutDir, simDataFile, provenance)

	@staticmethod
	def _multiples(bars, alive, control, design, colors, label, part, data,
			inputDir, plotOutDir, simDataFile, provenance):
		style.apply()
		ncols = 3
		nrows = int(np.ceil(len(bars) / ncols))
		fig, axes = plt.subplots(nrows, ncols, figsize=style.panel_size(
			'single', 1.15))
		notes = []
		for ax, (title, unit, column, scale) in zip(axes.flat, bars):
			for k, variant in enumerate((control, design)):
				mean, ci = _mean_ci(scale * alive[alive['variant']
					== variant][column].values)
				ax.bar(k, mean, width=0.7, color=colors[variant], yerr=ci
					if np.isfinite(ci) else None, error_kw=dict(
					ecolor=style.INK, elinewidth=0.6, capsize=1.5), zorder=2)
				notes.append('%s %s %.4g +/- %.3g' % (title, 'control' if
					variant == control else 'GFP', mean, ci))
			ax.set_title('%s\n%s' % (title, unit), fontsize=style.FONT_SIZES[
				'tick'])
			ax.set_xticks([])
			ax.set_xlim(-0.6, 1.6)
			ax.set_ylim(bottom=0)
			ax.yaxis.set_major_formatter(plt.FuncFormatter(
				lambda v, _: '%gM' % (v / 1e6) if v >= 1e6 else '%g' % v))
			ax.locator_params(axis='y', nbins=3)
		for ax in list(axes.flat)[len(bars):]:
			ax.set_visible(False)
		style.legend_below(fig, [Patch(color=colors[control], label='no GFP'),
			Patch(color=colors[design], label=label)], ncol=2)
		style.save_panel(fig, plotOutDir, style.panel_name(__file__, part),
			script=__file__, batch=inputDir, sim_data_path=simDataFile,
			n=int(alive[alive['variant'] == design].shape[0]),
			notes='Means over surviving lineages, t-based 95%% CI. Showcase ='
				' variant %d. %s.' % (design, '; '.join(notes)), **provenance)

	@staticmethod
	def _relative(bars, data, design, label, part, inputDir, plotOutDir,
			simDataFile, provenance):
		relative = [b for b in bars if b[2] in RELATIVE_COLORS]
		ratios = [gfp_paper_data.seed_matched_ratio(data.lineages, design,
			column) for _, _, column, _ in relative]
		style.apply()
		fig, ax = plt.subplots(figsize=style.panel_size('single', 0.7))
		style.ratio_bars(ax, [RELATIVE_LABELS[column]
			for _, _, column, _ in relative],
			[r['mean'] for r in ratios], [r['ci95'] for r in ratios],
			[RELATIVE_COLORS[column] for _, _, column, _ in relative])
		ax.set_ylabel('%s / no GFP' % label.replace(' of proteome', ''))
		name_part = 'relative' if part is None else 'relative_%s' % part
		style.save_panel(fig, plotOutDir, style.panel_name(__file__,
			name_part), script=__file__, batch=inputDir,
			sim_data_path=simDataFile, n=min(r['n'] for r in ratios),
			notes='Seed-matched ratios, showcase over control, mean with'
				' t-based 95%% CI: %s.' % ', '.join('%s %.3f +/- %.3f'
				% (title, r['mean'], r['ci95']) for (title, _, _, _), r in
				zip(relative, ratios)), **provenance)


if __name__ == '__main__':
	Plot().cli()
