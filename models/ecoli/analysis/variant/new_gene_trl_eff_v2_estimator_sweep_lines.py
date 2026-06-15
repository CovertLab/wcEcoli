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

Produces two PDFs:
  * <name>.pdf            -- 3 panels (completion, doubling time, dry mass),
                            block means only.
  * <name>_errorbars.pdf  -- 2 panels (doubling time, dry mass) with mean +/-
                            SD across seeds as error bars.

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


def _seed_of(cell_path):
	"""Seed directory component of a cell path (.../SEED/generation_*/*)."""
	parts = cell_path.rstrip('/').split(os.sep)
	for i, part in enumerate(parts):
		if part.startswith('generation_') and i > 0:
			return parts[i - 1]
	return parts[-3] if len(parts) >= 3 else cell_path


def _mean_std(by_seed, variant_indexes):
	"""Per-variant (mean, std) across the per-seed values (sample SD, ddof=1).

	std is 0.0 when only one seed has data and NaN when none do (so error
	bars vanish exactly where the mean is undefined).
	"""
	mean = {}
	std = {}
	for vi in variant_indexes:
		vals = [v for v in by_seed.get(vi, []) if np.isfinite(v)]
		if vals:
			mean[vi] = float(np.mean(vals))
			std[vi] = float(np.std(vals, ddof=1)) if len(vals) >= 2 else 0.0
		else:
			mean[vi] = np.nan
			std[vi] = np.nan
	return mean, std


def _collect_metrics(plot, variant_indexes, n_total_gens, ignore_first_n_gens):
	completion = {}
	avg_dt = {}
	avg_mass = {}
	# Per-seed values within the averaging window, for error bars.
	dt_by_seed = {}
	mass_by_seed = {}
	for vi in variant_indexes:
		n_gen0 = len(plot.ap.get_cells(
			variant=[vi], generation=[0], only_successful=True))
		n_last = len(plot.ap.get_cells(
			variant=[vi], generation=[n_total_gens - 1],
			only_successful=True))
		completion[vi] = (n_last / n_gen0) if n_gen0 > 0 else 0.0

		# Per-gen mean doubling time (+ per-seed values in the window)
		mean_dts = np.full(n_total_gens, np.nan)
		seed_dt = {}
		for gen in range(n_total_gens):
			cells = plot.ap.get_cells(
				variant=[vi], generation=[gen], only_successful=True)
			times = []
			for cell_path in cells:
				try:
					t = TableReader(
						os.path.join(cell_path, 'simOut', 'Main')
						).readColumn('time', squeeze=True)
				except Exception:
					continue
				dt_min = (t[-1] - t[0]) / 60.
				if dt_min <= MAX_DOUBLING_TIME_MIN:
					times.append(dt_min)
					if gen >= ignore_first_n_gens:
						seed_dt.setdefault(_seed_of(cell_path), []).append(
							dt_min)
			if times:
				mean_dts[gen] = float(np.mean(times))
		later = mean_dts[ignore_first_n_gens:]
		avg_dt[vi] = (float(np.nanmean(later))
			if np.any(np.isfinite(later)) else np.nan)
		dt_by_seed[vi] = [float(np.mean(v)) for v in seed_dt.values() if v]

		# Average dry mass (+ per-seed values)
		all_cells = plot.ap.get_cells(
			variant=[vi],
			generation=np.arange(ignore_first_n_gens, n_total_gens),
			only_successful=True)
		cell_means = []
		seed_mass = {}
		for cell_path in all_cells:
			try:
				dm = TableReader(
					os.path.join(cell_path, 'simOut', 'Mass')
					).readColumn('dryMass', squeeze=True)
			except Exception:
				continue
			m = float(np.mean(dm))
			cell_means.append(m)
			seed_mass.setdefault(_seed_of(cell_path), []).append(m)
		avg_mass[vi] = float(np.mean(cell_means)) if cell_means else np.nan
		mass_by_seed[vi] = [float(np.mean(v)) for v in seed_mass.values() if v]

	return completion, avg_dt, avg_mass, dt_by_seed, mass_by_seed


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

		completion, avg_dt, avg_mass, dt_by_seed, mass_by_seed = (
			_collect_metrics(
				self, variant_indexes, n_total_gens, ignore_first_n_gens))

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

		# --- Second PDF: doubling time and dry mass with error bars ---
		dt_mean, dt_std = _mean_std(dt_by_seed, variant_indexes)
		mass_mean, mass_std = _mean_std(mass_by_seed, variant_indexes)
		dt_mean_series = _series_by_block(dt_mean, variant_indexes)
		dt_std_series = _series_by_block(dt_std, variant_indexes)
		mass_mean_series = _series_by_block(mass_mean, variant_indexes)
		mass_std_series = _series_by_block(mass_std, variant_indexes)
		wt_dt_mean = dt_mean[wt_indexes[0]] if wt_indexes else np.nan
		wt_dt_std = dt_std[wt_indexes[0]] if wt_indexes else np.nan
		wt_mass_mean = mass_mean[wt_indexes[0]] if wt_indexes else np.nan
		wt_mass_std = mass_std[wt_indexes[0]] if wt_indexes else np.nan

		def _plot_errorbar(ax, mean_series, std_series, ylabel, title,
				wt_mean, wt_std, wt_fmt):
			for block_idx in range(len(BLOCK_ESTIMATORS)):
				ax.errorbar(
					x, np.array(mean_series[block_idx], dtype=float),
					yerr=np.array(std_series[block_idx], dtype=float),
					marker='o', lw=1.5, capsize=3,
					color=BLOCK_COLORS[block_idx],
					label=BLOCK_SHORT_NAMES[block_idx])
			if np.isfinite(wt_mean):
				ax.axhline(wt_mean, color=WILDTYPE_COLOR, ls='--', lw=1.0,
					label=f'wildtype ({wt_fmt.format(wt_mean)})')
				if np.isfinite(wt_std) and wt_std > 0:
					ax.fill_between(
						[float(np.min(x)), float(np.max(x))],
						wt_mean - wt_std, wt_mean + wt_std,
						color=WILDTYPE_COLOR, alpha=0.08)
			ax.set_xlabel('GFP translation efficiency')
			ax.set_ylabel(ylabel)
			ax.set_title(title)

		fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
		_plot_errorbar(
			axes2[0], dt_mean_series, dt_std_series,
			'Avg doubling time (min)',
			f'Avg doubling time vs translation efficiency\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1}; '
			f'mean ± SD across seeds)',
			wt_dt_mean, wt_dt_std, '{:.1f} min')
		axes2[0].legend(fontsize=8, loc='upper left')
		_plot_errorbar(
			axes2[1], mass_mean_series, mass_std_series,
			'Avg dry mass (fg)',
			f'Avg dry mass vs translation efficiency\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1}; '
			f'mean ± SD across seeds)',
			wt_mass_mean, wt_mass_std, '{:.2f} fg')
		axes2[1].legend(fontsize=8, loc='lower left')

		fig2.suptitle(
			'new_gene_trl_eff_v2_estimator_sweep: doubling time & dry mass '
			'(mean ± SD across seeds)',
			fontsize=12)
		fig2.tight_layout(rect=[0, 0, 1, 0.93])
		exportFigure(plt, plotOutDir, plotOutFileName + '_errorbars', metadata)

		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
