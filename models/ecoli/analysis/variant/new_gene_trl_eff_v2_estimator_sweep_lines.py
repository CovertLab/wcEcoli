"""
Line-plot summary of the new_gene_trl_eff_v2_estimator_sweep: three metrics
(completion rate, average doubling time, average dry mass) plotted as a
function of GFP translation efficiency, with one line per block (gfp_only,
gap, smoothed_max).  Wildtype is drawn as a horizontal reference line in each
panel.

Sister plot to new_gene_trl_eff_v2_estimator_sweep_comparison.py (which shows
the same data as bars).  This view compares the kcat estimators head-to-head
at matching translation efficiencies, where the metabolic stress comes from
GFP ribosome competition (multiplier fixed at 1.0).

Adapted from kcat_v2_estimator_sweep_lines.py.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_v2_estimator_sweep import (
	BLOCK_ESTIMATORS,
	BLOCK_SHORT_NAMES,
	TRL_EFF_VALUES,
	N_PER_BLOCK,
	is_wildtype,
	variant_to_block,
)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


MAX_DOUBLING_TIME_MIN = 300

# One color per block (index into BLOCK_ESTIMATORS / BLOCK_SHORT_NAMES).
BLOCK_COLORS = ('#7f7f7f', '#1f77b4', '#2ca02c')  # gfp_only, gap, smoothed_max
WILDTYPE_COLOR = '#222222'


def _collect_metrics(plot, variant_indexes, n_total_gens, ignore_first_n_gens):
	completion = {}
	avg_dt = {}
	avg_mass = {}
	for vi in variant_indexes:
		n_gen0 = len(plot.ap.get_cells(
			variant=[vi], generation=[0], only_successful=True))
		n_last = len(plot.ap.get_cells(
			variant=[vi], generation=[n_total_gens - 1],
			only_successful=True))
		completion[vi] = (n_last / n_gen0) if n_gen0 > 0 else 0.0

		# Per-gen mean doubling time
		mean_dts = np.full(n_total_gens, np.nan)
		for gen in range(n_total_gens):
			cells = plot.ap.get_cells(
				variant=[vi], generation=[gen], only_successful=True)
			times = []
			for cell_path in cells:
				try:
					t = TableReader(
						os.path.join(cell_path, 'simOut', 'Main')
						).readColumn('time', squeeze=True)
					dt_min = (t[-1] - t[0]) / 60.
					if dt_min <= MAX_DOUBLING_TIME_MIN:
						times.append(dt_min)
				except Exception:
					continue
			if times:
				mean_dts[gen] = float(np.mean(times))
		later = mean_dts[ignore_first_n_gens:]
		avg_dt[vi] = (float(np.nanmean(later))
			if np.any(np.isfinite(later)) else np.nan)

		# Average dry mass
		all_cells = plot.ap.get_cells(
			variant=[vi],
			generation=np.arange(ignore_first_n_gens, n_total_gens),
			only_successful=True)
		cell_means = []
		for cell_path in all_cells:
			try:
				dm = TableReader(
					os.path.join(cell_path, 'simOut', 'Mass')
					).readColumn('dryMass', squeeze=True)
				cell_means.append(float(np.mean(dm)))
			except Exception:
				continue
		avg_mass[vi] = float(np.mean(cell_means)) if cell_means else np.nan

	return completion, avg_dt, avg_mass


def _series_by_block(metric, variant_indexes):
	"""Return {block_idx: list-of-value-by-trl_eff-index}."""
	by_block = {b: [np.nan] * N_PER_BLOCK
		for b in range(len(BLOCK_ESTIMATORS))}
	for vi in variant_indexes:
		if is_wildtype(vi):
			continue
		block_idx, te_idx = variant_to_block(vi)
		by_block[block_idx][te_idx] = metric.get(vi, np.nan)
	return by_block


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)

		completion, avg_dt, avg_mass = _collect_metrics(
			self, variant_indexes, n_total_gens, ignore_first_n_gens)

		# Wildtype values for reference lines
		wt_indexes = [vi for vi in variant_indexes if is_wildtype(vi)]
		wt_completion = (completion[wt_indexes[0]]
			if wt_indexes else np.nan)
		wt_dt = avg_dt[wt_indexes[0]] if wt_indexes else np.nan
		wt_mass = avg_mass[wt_indexes[0]] if wt_indexes else np.nan

		comp_series = _series_by_block(completion, variant_indexes)
		dt_series = _series_by_block(avg_dt, variant_indexes)
		mass_series = _series_by_block(avg_mass, variant_indexes)

		x = np.array(TRL_EFF_VALUES)

		fig, axes = plt.subplots(1, 3, figsize=(16, 5))

		def _plot_metric(ax, series, ylabel, title, wt_value, wt_fmt):
			for block_idx in range(len(BLOCK_ESTIMATORS)):
				ax.plot(x, series[block_idx],
					marker='o', lw=1.5,
					color=BLOCK_COLORS[block_idx],
					label=BLOCK_SHORT_NAMES[block_idx])
			if np.isfinite(wt_value):
				ax.axhline(wt_value, color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=f'wildtype ({wt_fmt.format(wt_value)})')
			ax.set_xlabel('GFP translation efficiency')
			ax.set_ylabel(ylabel)
			ax.set_title(title)

		# Panel 1: completion rate
		_plot_metric(
			axes[0], comp_series, 'Completion rate',
			'Completion rate vs translation efficiency',
			wt_completion, '{:.0%}')
		axes[0].set_ylim(0, 1.05)
		axes[0].legend(fontsize=8, loc='lower left')

		# Panel 2: average doubling time
		_plot_metric(
			axes[1], dt_series, 'Avg doubling time (min)',
			f'Avg doubling time vs translation efficiency\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1})',
			wt_dt, '{:.1f} min')
		axes[1].legend(fontsize=8, loc='upper left')

		# Panel 3: average dry mass
		_plot_metric(
			axes[2], mass_series, 'Avg dry mass (fg)',
			f'Avg dry mass vs translation efficiency\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1})',
			wt_mass, '{:.2f} fg')
		axes[2].legend(fontsize=8, loc='lower left')

		fig.suptitle(
			'new_gene_trl_eff_v2_estimator_sweep: estimator comparison',
			fontsize=12)
		fig.tight_layout(rect=[0, 0, 1, 0.94])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
