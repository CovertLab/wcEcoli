"""
Summary plot for the kcat_v2_estimator_sweep variant comparing doubling-time
impacts across three v2 kcat estimators (gap, smoothed_max, drop_top_20) at a
uniform multiplier grid, plus a wildtype control.

Adapted from new_gene_trl_eff_sweep_comparison.py.  Produces four panels:

  1. Completion rate per variant (bar; colored by estimator family, lighter
     shades for smaller multipliers).
  2. Mean doubling time per generation (heatmap; variants x generations).
  3. Average doubling time over the last 8 generations (bar).
  4. Average dry mass over the last 8 generations (bar).

Also writes doubling_times_by_lineage.csv with one row per (variant, seed)
containing each generation's doubling time.
"""

import csv
import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import matplotlib.ticker as mticker

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


# Cells longer than this are considered timed-out and excluded from doubling
# time statistics.
MAX_DOUBLING_TIME_MIN = 300

# One base color per estimator family.  Multipliers shade lighter -> darker
# toward smaller multipliers (more aggressive trimming), so within a family
# the leftmost (multiplier = 1.0) bar is the lightest.
ESTIMATOR_BASE_COLORS = {
	'v2_gap':         (0.122, 0.467, 0.706),  # blue   (#1f77b4)
	'v2_smoothed_max': (0.173, 0.627, 0.173),  # green  (#2ca02c)
	'v2_drop_top_20': (1.000, 0.498, 0.055),  # orange (#ff7f0e)
}
WILDTYPE_COLOR = (0.1, 0.1, 0.1)  # near-black


def _shade(base_rgb, multiplier):
	"""Mix base color with white based on multiplier.

	multiplier = 1.0 -> 80% white + 20% base (lightest)
	multiplier = 0.1 -> 0% white + 100% base (darkest)
	"""
	# Map multiplier in [0.1, 1.0] to mix factor in [0.0, 0.8].
	frac = max(0.0, min(0.8, (multiplier - 0.1) * (0.8 / 0.9)))
	r, g, b = base_rgb
	return (
		r + frac * (1 - r),
		g + frac * (1 - g),
		b + frac * (1 - b),
	)


def _variant_color(index):
	if is_wildtype(index):
		return WILDTYPE_COLOR
	est_idx, mult_idx = variant_to_estimator(index)
	return _shade(
		ESTIMATOR_BASE_COLORS[ESTIMATOR_LABELS[est_idx]],
		MULTIPLIERS[mult_idx],
	)


def _variant_label(index):
	if is_wildtype(index):
		return 'WT'
	est_idx, mult_idx = variant_to_estimator(index)
	pct = int(round(MULTIPLIERS[mult_idx] * 100))
	return f'{ESTIMATOR_SHORT_NAMES[est_idx]}\nx{pct}%'


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant_indexes = self.ap.get_variants()
		n_variants = len(variant_indexes)
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)

		completion_rates = {}
		dt_by_gen = {}
		avg_doubling_time = {}
		avg_dry_mass = {}
		lineage_dts = []

		for vi in variant_indexes:
			# Completion rate
			n_gen0 = len(self.ap.get_cells(
				variant=[vi], generation=[0], only_successful=True))
			n_last = len(self.ap.get_cells(
				variant=[vi], generation=[n_total_gens - 1],
				only_successful=True))
			completion_rates[vi] = (n_last / n_gen0) if n_gen0 > 0 else 0.0

			# Per-generation doubling time
			seed_dts = {}
			mean_dts = np.full(n_total_gens, np.nan)
			for gen in range(n_total_gens):
				cells = self.ap.get_cells(
					variant=[vi], generation=[gen], only_successful=True)
				if len(cells) == 0:
					continue
				times = []
				for cell_path in cells:
					sim_out = os.path.join(cell_path, 'simOut')
					try:
						t = TableReader(
							os.path.join(sim_out, 'Main')
							).readColumn('time', squeeze=True)
						dt_min = (t[-1] - t[0]) / 60.
					except Exception:
						continue
					parts = cell_path.rstrip('/').split('/')
					seed = parts[-3]
					seed_dts.setdefault(seed, {})[gen] = dt_min
					if dt_min <= MAX_DOUBLING_TIME_MIN:
						times.append(dt_min)
				if times:
					mean_dts[gen] = np.mean(times)
			dt_by_gen[vi] = mean_dts

			later_dts = mean_dts[ignore_first_n_gens:]
			avg_doubling_time[vi] = (
				np.nanmean(later_dts) if np.any(np.isfinite(later_dts))
				else np.nan)

			for seed in sorted(seed_dts.keys()):
				row = {'variant': vi, 'seed': seed}
				for gen in range(n_total_gens):
					row[f'gen_{gen}_dt_min'] = seed_dts[seed].get(gen, np.nan)
				lineage_dts.append(row)

			# Average dry mass over the last 8 gens
			all_cells = self.ap.get_cells(
				variant=[vi],
				generation=np.arange(ignore_first_n_gens, n_total_gens),
				only_successful=True)
			cell_means = []
			for cell_path in all_cells:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					dm = TableReader(
						os.path.join(sim_out, 'Mass')
						).readColumn('dryMass', squeeze=True)
					cell_means.append(np.mean(dm))
				except Exception:
					pass
			avg_dry_mass[vi] = np.mean(cell_means) if cell_means else np.nan

		# Per-lineage CSV
		if lineage_dts:
			gen_cols = [f'gen_{g}_dt_min' for g in range(n_total_gens)]
			csv_path = os.path.join(plotOutDir, 'doubling_times_by_lineage.csv')
			with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
				writer = csv.DictWriter(
					fh, fieldnames=['variant', 'seed'] + gen_cols)
				writer.writeheader()
				writer.writerows(lineage_dts)
			print(f'Wrote {csv_path}')

		# Figure
		fig = plt.figure(figsize=(max(14, n_variants * 0.55), 20))
		gs = fig.add_gridspec(4, 1, hspace=0.5,
							  height_ratios=[1, 2.5, 1, 1])

		x = np.arange(n_variants)
		labels = [_variant_label(vi) for vi in variant_indexes]
		colors = [_variant_color(vi) for vi in variant_indexes]

		# Panel 1: completion rate
		ax1 = fig.add_subplot(gs[0])
		bars = ax1.bar(x, [completion_rates[vi] for vi in variant_indexes],
					   color=colors, edgecolor='white', linewidth=0.5)
		ax1.axhline(1.0, color='k', lw=0.8, ls='--')
		ax1.set_ylim(0, 1.12)
		ax1.set_ylabel('Completion rate', fontsize=10)
		ax1.set_title(
			'Fraction of seeds completing all generations '
			'(kcat_v2_estimator_sweep)',
			fontsize=11)
		ax1.set_xticks(x)
		ax1.set_xticklabels(labels, fontsize=7, rotation=90)
		ax1.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
		for bar, vi in zip(bars, variant_indexes):
			rate = completion_rates[vi]
			ax1.text(bar.get_x() + bar.get_width() / 2,
					 rate + 0.02, f'{rate:.0%}',
					 ha='center', va='bottom', fontsize=6)

		# Panel 2: heatmap of mean dt by generation
		ax2 = fig.add_subplot(gs[1])
		heatmap_data = np.array([dt_by_gen[vi] for vi in variant_indexes])
		finite = heatmap_data[np.isfinite(heatmap_data)]
		if finite.size:
			vmin = float(np.nanmin(finite))
			vmax = min(float(np.nanmax(finite)), MAX_DOUBLING_TIME_MIN)
		else:
			vmin, vmax = 0, MAX_DOUBLING_TIME_MIN
		im = ax2.imshow(
			heatmap_data, aspect='auto', origin='upper', cmap='RdYlGn_r',
			vmin=vmin, vmax=vmax, interpolation='nearest')
		cbar = fig.colorbar(im, ax=ax2, fraction=0.03, pad=0.02)
		cbar.set_label('Mean doubling time (min)', fontsize=9)
		ax2.set_xlabel('Generation', fontsize=10)
		ax2.set_ylabel('Variant', fontsize=10)
		ax2.set_title(
			'Mean doubling time per generation\n'
			f'(cells > {MAX_DOUBLING_TIME_MIN} min excluded; gray = no data)',
			fontsize=11)
		ax2.set_xticks(np.arange(n_total_gens))
		ax2.set_xticklabels(np.arange(n_total_gens), fontsize=7)
		ax2.set_yticks(np.arange(n_variants))
		ax2.set_yticklabels(labels, fontsize=7)
		nan_mask = np.isnan(heatmap_data)
		for (vi_idx, gen_idx) in zip(*np.where(nan_mask)):
			ax2.add_patch(plt.Rectangle(
				(gen_idx - 0.5, vi_idx - 0.5), 1, 1,
				color='#cccccc', zorder=2))

		# Panel 3: average doubling time
		ax3 = fig.add_subplot(gs[2])
		dts = [avg_doubling_time[vi] for vi in variant_indexes]
		ax3.bar(x, dts, color=colors, edgecolor='white', linewidth=0.5)
		finite_dts = [d for d in dts if np.isfinite(d)]
		if finite_dts:
			ax3.set_ylim(min(finite_dts) * 0.97, max(finite_dts) * 1.03)
		ax3.set_ylabel('Average doubling time (min)', fontsize=10)
		ax3.set_title(
			f'Mean doubling time averaged over seeds\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1}; '
			f'cells > {MAX_DOUBLING_TIME_MIN} min excluded)',
			fontsize=11)
		ax3.set_xticks(x)
		ax3.set_xticklabels(labels, fontsize=7, rotation=90)

		# Panel 4: average dry mass
		ax4 = fig.add_subplot(gs[3])
		masses = [avg_dry_mass[vi] for vi in variant_indexes]
		ax4.bar(x, masses, color=colors, edgecolor='white', linewidth=0.5)
		finite_masses = [m for m in masses if np.isfinite(m)]
		if finite_masses:
			ax4.set_ylim(min(finite_masses) * 0.97, max(finite_masses) * 1.03)
		ax4.set_ylabel('Average dry mass (fg)', fontsize=10)
		ax4.set_title(
			f'Mean dry mass averaged over timesteps and seeds\n'
			f'(gens {ignore_first_n_gens}–{n_total_gens - 1})',
			fontsize=11)
		ax4.set_xticks(x)
		ax4.set_xticklabels(labels, fontsize=7, rotation=90)

		# Legend: one entry per estimator + wildtype
		legend_elements = [Patch(facecolor=WILDTYPE_COLOR, label='wildtype')]
		for est_idx, label in enumerate(ESTIMATOR_LABELS):
			legend_elements.append(Patch(
				facecolor=ESTIMATOR_BASE_COLORS[label],
				label=f'{ESTIMATOR_SHORT_NAMES[est_idx]} '
				f'(darker = smaller multiplier)'))
		fig.legend(handles=legend_elements, loc='lower center',
				   ncol=min(len(legend_elements), 4), fontsize=8,
				   bbox_to_anchor=(0.5, -0.01))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
