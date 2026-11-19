"""
Growth against GFP's share of the proteome, one point per design.

# Figure 2C

Each point is a design's window mean over its surviving lineages; hollow
points are designs where some lineages failed. A straight line is fitted
(least squares) to growth rate relative to the seed-matched control against
GFP share, using the no-GFP control (1 at 0%, by definition) and every design
where all lineages survived; it is solid over the range it was fitted to and
dashed beyond. The slope and intercept are in the
sidecar. Two versions while decision D4 is open:
	growth rate     relative growth rate, the fitted line, a dashed line at 1
	doubling time   minutes (suffix DOUBLING_TIME_PART); the same fit drawn
	                as control doubling time / fitted relative growth rate,
	                and a dashed line at the control mean

Experimental points are added once the data exist (EXPERIMENT_FILE); until
then the panel shows simulation only and says so in the sidecar.

The supporting claim: growth falls steadily as GFP takes more of the
proteome, over a range of reductions consistent with experiments.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis import gfp_paper_constants as C
from wholecell.analysis import gfp_paper_data
from wholecell.analysis import gfp_paper_style as style

PERCENT = 100.0
# Columns gfp_share_pct, growth_rel_wt; drawn in the experiment colour when
# present. Path is relative to the repository root.
EXPERIMENT_FILE = os.path.join('runscripts', 'gfp_paper', 'data',
	'growth_vs_gfp_share_experiments.tsv')


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		data = gfp_paper_data.load(inputDir, simDataFile, metadata,
			plotOutDir, self.cpus)
		designs = data.designs[data.designs['promoter'] > 0]
		designs = designs[np.isfinite(designs['gfp_proteome_share'])
			& np.isfinite(designs['relative_growth_rate'])]
		x = PERCENT * designs['gfp_proteome_share'].values
		complete = ~designs['failed'].values
		fit_x = np.append(0.0, x[complete])
		fit_y = np.append(1.0, designs['relative_growth_rate'].values[complete])
		if fit_x.size >= 2:
			slope, intercept = np.polyfit(fit_x, fit_y, 1)
		else:
			slope = intercept = np.nan
		control_dt = float(data.designs.loc[data.designs['variant']
			== gfp_paper_data.CONTROL_VARIANT, 'doubling_time_min'].iloc[0])
		x_fit_end = float(x[complete].max()) if complete.any() else np.nan
		x_end = max(75.0, float(np.nanmax(x)) if len(x) else 75.0)

		versions = (
			(None, 'relative_growth_rate', 'growth rate relative to no GFP',
				1.0, lambda rel: rel),
			(C.DOUBLING_TIME_PART, 'doubling_time_min', 'doubling time (min)',
				control_dt, lambda rel: control_dt / rel),
		)
		for part, column, ylabel, reference, transform in versions:
			style.apply()
			fig, ax = plt.subplots(figsize=style.panel_size('single', 0.8))
			color = style.COLORS['gfp']
			if np.isfinite(slope):
				solid = np.linspace(0, x_fit_end, 50)
				dashed = np.linspace(x_fit_end, x_end, 50)
				with np.errstate(divide='ignore', invalid='ignore'):
					ys = transform(intercept + slope * solid)
					yd = transform(intercept + slope * dashed)
				ax.plot(solid, ys, color=color, linewidth=1.2, zorder=2,
					label='simulation fit')
				keep = np.isfinite(yd) & (yd > 0)
				ax.plot(dashed[keep], yd[keep], color=color, linewidth=1.0,
					linestyle=(0, (4, 3)), zorder=2)
			style.design_scatter(ax, x, designs[column].values, color,
				failed=designs['failed'].values, label='simulation', zorder=3)
			ax.axhline(reference, color=style.INK_MUTED,
				linestyle=(0, (4, 3)), linewidth=0.8, zorder=1)
			experiments = self._experiments()
			if experiments is not None:
				ax.scatter(experiments['gfp_share_pct'], transform(
					experiments['growth_rel_wt']), s=14,
					color=style.COLORS['enzyme'], zorder=3,
					label='experiments')
			ax.set_xlim(0, x_end)
			if part is None:
				ax.set_ylim(0, max(1.2, float(np.nanmax(designs[column]))
					* 1.05))
			else:
				ax.set_ylim(0, float(np.nanmax(designs[column])) * 1.15)
			ax.set_xlabel('GFP share of proteome (%)')
			ax.set_ylabel(ylabel)
			style.legend_with_failed(ax, bool(designs['failed'].any()),
				loc='lower left' if part is None else 'upper left')

			style.save_panel(fig, plotOutDir,
				style.panel_name(__file__, part), script=__file__,
				batch=inputDir, sim_data_path=simDataFile,
				n=int(len(designs)),
				notes='%s. Fit on the control and %d complete designs:'
					' relative growth ='
					' %.4f + %.5f x GFP share (%%). %s'
					% (ylabel, int(complete.sum()), intercept, slope,
					'Experimental points included.' if experiments
					is not None else 'No experimental data yet (%s).'
					% EXPERIMENT_FILE),
				**data.provenance(sorted(int(v) for v in designs['variant'])))

	@staticmethod
	def _experiments():
		"""The experimental table, or None until it exists."""
		import pandas as pd
		from wholecell.utils import filepath
		path = os.path.join(filepath.ROOT_PATH, EXPERIMENT_FILE)
		if not os.path.exists(path):
			return None
		return pd.read_csv(path, sep='\t')


if __name__ == '__main__':
	Plot().cli()
