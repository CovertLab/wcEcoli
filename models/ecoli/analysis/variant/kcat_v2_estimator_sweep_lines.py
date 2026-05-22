"""
Line-plot summary of the kcat_v2_estimator_sweep: three metrics
(completion rate, average doubling time, average dry mass) plotted as a
function of multiplier, with one line per estimator family.  Wildtype is
drawn as a horizontal reference line in each panel.

Sister plot to kcat_v2_estimator_sweep_comparison.py (which shows the same
data as bars).  This view is the cleanest way to compare estimators
head-to-head at matching multipliers.
"""

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.kcat_v2_estimator_sweep import (
	ESTIMATOR_LABELS,
	ESTIMATOR_SHORT_NAMES,
	MULTIPLIERS,
	N_PER_ESTIMATOR,
	is_wildtype,
	variant_to_estimator,
)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


MAX_DOUBLING_TIME_MIN = 300

ESTIMATOR_BASE_COLORS = {
	'v2_gap':          '#1f77b4',
	'v2_smoothed_max': '#2ca02c',
	'v2_drop_top_20':  '#ff7f0e',
}
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


def _series_by_estimator(metric, variant_indexes):
	"""Return {estimator_label: list-of-(multiplier, value)-by-multiplier}."""
	by_est = {label: [None] * N_PER_ESTIMATOR for label in ESTIMATOR_LABELS}
	for vi in variant_indexes:
		if is_wildtype(vi):
			continue
		est_idx, mult_idx = variant_to_estimator(vi)
		by_est[ESTIMATOR_LABELS[est_idx]][mult_idx] = metric.get(vi, np.nan)
	# Replace None with NaN for plotting safety
	for label, series in by_est.items():
		by_est[label] = [v if v is not None else np.nan for v in series]
	return by_est


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

		comp_series = _series_by_estimator(completion, variant_indexes)
		dt_series = _series_by_estimator(avg_dt, variant_indexes)
		mass_series = _series_by_estimator(avg_mass, variant_indexes)

		x = np.array(MULTIPLIERS)

		fig, axes = plt.subplots(1, 3, figsize=(16, 5))

		# Panel 1: completion rate
		ax = axes[0]
		for label in ESTIMATOR_LABELS:
			ax.plot(x, comp_series[label],
				marker='o', lw=1.5,
				color=ESTIMATOR_BASE_COLORS[label],
				label=ESTIMATOR_SHORT_NAMES[ESTIMATOR_LABELS.index(label)])
		if np.isfinite(wt_completion):
			ax.axhline(wt_completion, color=WILDTYPE_COLOR, ls='--', lw=1.0,
				label=f'wildtype ({wt_completion:.0%})')
		ax.set_xlabel('Multiplier')
		ax.set_ylabel('Completion rate')
		ax.set_ylim(0, 1.05)
		ax.set_title('Completion rate vs multiplier')
		ax.invert_xaxis()  # tighter bound (smaller multiplier) on the right
		ax.legend(fontsize=8, loc='lower right')

		# Panel 2: average doubling time
		ax = axes[1]
		for label in ESTIMATOR_LABELS:
			ax.plot(x, dt_series[label],
				marker='o', lw=1.5,
				color=ESTIMATOR_BASE_COLORS[label],
				label=ESTIMATOR_SHORT_NAMES[ESTIMATOR_LABELS.index(label)])
		if np.isfinite(wt_dt):
			ax.axhline(wt_dt, color=WILDTYPE_COLOR, ls='--', lw=1.0,
				label=f'wildtype ({wt_dt:.1f} min)')
		ax.set_xlabel('Multiplier')
		ax.set_ylabel('Avg doubling time (min)')
		ax.set_title(
			f'Avg doubling time vs multiplier\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1})')
		ax.invert_xaxis()
		ax.legend(fontsize=8, loc='upper right')

		# Panel 3: average dry mass
		ax = axes[2]
		for label in ESTIMATOR_LABELS:
			ax.plot(x, mass_series[label],
				marker='o', lw=1.5,
				color=ESTIMATOR_BASE_COLORS[label],
				label=ESTIMATOR_SHORT_NAMES[ESTIMATOR_LABELS.index(label)])
		if np.isfinite(wt_mass):
			ax.axhline(wt_mass, color=WILDTYPE_COLOR, ls='--', lw=1.0,
				label=f'wildtype ({wt_mass:.2f} fg)')
		ax.set_xlabel('Multiplier')
		ax.set_ylabel('Avg dry mass (fg)')
		ax.set_title(
			f'Avg dry mass vs multiplier\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1})')
		ax.invert_xaxis()
		ax.legend(fontsize=8, loc='lower right')

		fig.suptitle(
			'kcat_v2_estimator_sweep: estimator comparison',
			fontsize=12)
		fig.tight_layout(rect=[0, 0, 1, 0.94])
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
