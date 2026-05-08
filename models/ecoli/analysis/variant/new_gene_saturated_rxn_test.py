"""
Per-condition summary for the saturated-reaction causality test variant
(`new_gene_saturated_rxn_test`).

Aggregates doubling time and average ppGpp concentration across seeds for
each of the 17 conditions, emits a structured CSV, and plots:
  1. Bar chart: doubling time per condition with seed std error bars,
     color-coded by group (no-gfp / gfp+unconstrained / gfp+full kcat /
     gfp+4-rxn / gfp+16-rxn / no-gfp+constrained), and a dashed reference
     line at the no-gfp-unconstrained dt.
  2. dt vs trl_eff scatter for the four gfp constraint scopes
     (unconstrained / full / 4-rxn / 16-rxn) -- the primary causality plot.

Inputs come straight from simOut tables; no dependency on any prior analysis
script.
"""

import csv
import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_saturated_rxn_test import (
	VARIANTS,
	N_VARIANTS,
	TRL_EFF_LEVELS,
	EXPRESSION_FACTOR,
	_condition_label,
)
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)

# Use last 8 generations after new gene induction at gen 8 (24 total gens).
N_GENS_TO_USE = 8

# Group encoding used for both the bar-plot color and the scatter series.
# Maps (gfp, rxn_set) -> (group_name, base_color).
GROUP_STYLE = {
	(False, None):     ('no-gfp unconstrained',     '#999999'),
	(True,  None):     ('gfp unconstrained',        '#1f77b4'),  # blue
	(True,  'full'):   ('gfp + full kcat 1.0x',     '#9467bd'),  # purple
	(True,  'rxns4'):  ('gfp + 4-rxn 1.0x',         '#ff7f0e'),  # orange
	(True,  'rxns16'): ('gfp + 16-rxn 1.0x',        '#d62728'),  # red
	(False, 'rxns4'):  ('no-gfp + 4-rxn',           '#a1d99b'),  # light green
	(False, 'rxns16'): ('no-gfp + 16-rxn',          '#31a354'),  # dark green
}


def _shade(base_hex, depth):
	"""Return a slightly darker hex shade for trl_eff depth in [0, n-1]."""
	r, g, b = (int(base_hex[i:i + 2], 16) for i in (1, 3, 5))
	# Darken progressively for higher depth.
	scale = 1.0 - 0.18 * depth
	r, g, b = (max(0, int(c * scale)) for c in (r, g, b))
	return f'#{r:02x}{g:02x}{b:02x}'


def _bar_color(index):
	"""Pick a shade for one variant index, deepening with trl_eff."""
	gfp, rxn_set, multiplier, trl_eff = VARIANTS[index]
	_label, base = GROUP_STYLE[(gfp, rxn_set)]
	# Shading by trl_eff for gfp variants, by multiplier for no-gfp constrained.
	if gfp and trl_eff is not None:
		depth = TRL_EFF_LEVELS.index(trl_eff)
	elif (not gfp) and rxn_set is not None and multiplier == 0.4:
		depth = 1  # darken the 0.4x bar slightly vs the 1.0x bar
	else:
		depth = 0
	return _shade(base, depth)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variant_indexes = self.ap.get_variants()
		available = set(variant_indexes)
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - N_GENS_TO_USE, 0)
		gen_range = np.arange(ignore_first_n_gens, n_total_gens)

		dt_mean = np.full(N_VARIANTS, np.nan)
		dt_std = np.full(N_VARIANTS, np.nan)
		dt_n = np.zeros(N_VARIANTS, dtype=int)
		ppgpp_mean = np.full(N_VARIANTS, np.nan)

		for vi in range(N_VARIANTS):
			if vi not in available:
				continue
			cells = self.ap.get_cells(
				variant=[vi],
				generation=gen_range,
				only_successful=True)
			if len(cells) == 0:
				continue

			dt = read_stacked_columns(
				cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			dt = np.atleast_1d(dt)
			dt_mean[vi] = float(np.mean(dt))
			dt_std[vi] = float(np.std(dt))
			dt_n[vi] = int(dt.size)

			ppgpp = read_stacked_columns(
				cells, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True, ignore_exception=True)
			if ppgpp is not None and ppgpp.size > 0:
				ppgpp_mean[vi] = float(np.mean(ppgpp))

		# -------------------------- CSV ---------------------------- #
		csv_path = os.path.join(plotOutDir, 'saturated_rxn_test_summary.csv')
		with open(csv_path, 'w', newline='') as f:
			writer = csv.writer(f)
			writer.writerow([
				'variant', 'condition', 'gfp', 'rxn_set', 'multiplier',
				'trl_eff', 'dt_mean_min', 'dt_std_min', 'n_cells',
				'avg_ppgpp_uM'])
			for vi in range(N_VARIANTS):
				gfp, rxn_set, multiplier, trl_eff = VARIANTS[vi]
				writer.writerow([
					vi,
					_condition_label(vi),
					int(gfp),
					rxn_set if rxn_set is not None else '',
					'' if multiplier is None else f'{multiplier}',
					'' if trl_eff is None else f'{trl_eff}',
					'' if np.isnan(dt_mean[vi]) else f'{dt_mean[vi]:.4f}',
					'' if np.isnan(dt_std[vi]) else f'{dt_std[vi]:.4f}',
					int(dt_n[vi]),
					'' if np.isnan(ppgpp_mean[vi]) else f'{ppgpp_mean[vi]:.6f}',
				])

		# Reference dt for the "no-gfp unconstrained" baseline (idx 0).
		ref_dt = dt_mean[0] if not np.isnan(dt_mean[0]) else np.nan

		# ------------------------ Figure --------------------------- #
		fig, (ax_bar, ax_scatter) = plt.subplots(
			2, 1, figsize=(max(14, 0.7 * N_VARIANTS + 4), 12))

		# --- Bar plot --- #
		xs = np.arange(N_VARIANTS)
		colors = [_bar_color(vi) for vi in range(N_VARIANTS)]
		bars = ax_bar.bar(
			xs, np.nan_to_num(dt_mean, nan=0.0),
			yerr=np.nan_to_num(dt_std, nan=0.0),
			color=colors, edgecolor='black', linewidth=0.4, capsize=3)
		# Hide bars for variants with no data (dt_mean still NaN).
		for vi, bar in enumerate(bars):
			if np.isnan(dt_mean[vi]):
				bar.set_alpha(0.2)

		ax_bar.set_xticks(xs)
		ax_bar.set_xticklabels(
			[_condition_label(vi) for vi in range(N_VARIANTS)],
			rotation=90, ha='right', fontsize=8)
		ax_bar.set_ylabel('Doubling time (min)')
		ax_bar.set_title(
			f'Saturated-rxn causality test: dt per condition '
			f'(last {N_GENS_TO_USE} gens of {n_total_gens})')
		if not np.isnan(ref_dt):
			ax_bar.axhline(ref_dt, linestyle='--', color='gray', linewidth=1,
				label=f'no-gfp unconstrained ({ref_dt:.1f} min)')

		# Group legend: one swatch per (gfp, rxn_set) tuple actually present.
		group_handles = []
		seen_groups = set()
		for vi in range(N_VARIANTS):
			gfp, rxn_set, _m, _t = VARIANTS[vi]
			key = (gfp, rxn_set)
			if key in seen_groups:
				continue
			seen_groups.add(key)
			group_label, base = GROUP_STYLE[key]
			group_handles.append(plt.Rectangle(
				(0, 0), 1, 1, color=base, label=group_label))
		ax_bar.legend(handles=group_handles, fontsize=8, loc='upper left',
			ncol=2)

		# --- dt vs trl_eff scatter (gfp variants only) --- #
		series = {
			'unconstrained': [],
			'full kcat 1.0x': [],
			'4-rxn 1.0x': [],
			'16-rxn 1.0x': [],
		}
		series_color = {
			'unconstrained': GROUP_STYLE[(True, None)][1],
			'full kcat 1.0x': GROUP_STYLE[(True, 'full')][1],
			'4-rxn 1.0x': GROUP_STYLE[(True, 'rxns4')][1],
			'16-rxn 1.0x': GROUP_STYLE[(True, 'rxns16')][1],
		}
		for vi in range(N_VARIANTS):
			gfp, rxn_set, _multiplier, trl_eff = VARIANTS[vi]
			if not gfp or trl_eff is None or np.isnan(dt_mean[vi]):
				continue
			if rxn_set is None:
				series['unconstrained'].append((trl_eff, dt_mean[vi], dt_std[vi]))
			elif rxn_set == 'full':
				series['full kcat 1.0x'].append((trl_eff, dt_mean[vi], dt_std[vi]))
			elif rxn_set == 'rxns4':
				series['4-rxn 1.0x'].append((trl_eff, dt_mean[vi], dt_std[vi]))
			elif rxn_set == 'rxns16':
				series['16-rxn 1.0x'].append((trl_eff, dt_mean[vi], dt_std[vi]))

		for name, points in series.items():
			if not points:
				continue
			points.sort()
			xs_, ys_, ye_ = zip(*points)
			ax_scatter.errorbar(
				xs_, ys_, yerr=ye_,
				marker='o', linewidth=2, markersize=8,
				capsize=3, color=series_color[name], label=name)

		if not np.isnan(ref_dt):
			ax_scatter.axhline(
				ref_dt, linestyle='--', color='gray', linewidth=1,
				label=f'no-gfp unconstrained ({ref_dt:.1f} min)')
		ax_scatter.set_xlabel('Translation efficiency')
		ax_scatter.set_ylabel('Doubling time (min)')
		ax_scatter.set_title(
			'dt vs trl_eff for gfp variants -- '
			'subset constraints vs full constraint vs unconstrained')
		ax_scatter.set_xticks(TRL_EFF_LEVELS)
		ax_scatter.legend(fontsize=9, loc='upper left')
		ax_scatter.grid(True, alpha=0.3)

		fig.tight_layout()
		exportFigure(fig, plotOutDir, plotOutFileName, metadata)
		plt.close(fig)


if __name__ == '__main__':
	Plot().cli()
