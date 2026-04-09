"""
Summary plot for the new_gene_trl_eff_sweep variant comparing doubling time
impacts with and without kcat constraints at matching translation efficiencies.

Produces four panels:

  1. Completion rate  -- fraction of seeds that reached the final generation,
     one bar per variant.  Blue hues for no-kcat, red/orange hues for kcat.

  2. Mean doubling time by generation  -- heatmap (variants x generations).
     Cells that timed out (> MAX_DOUBLING_TIME_MIN minutes) are excluded.

  3. Average doubling time  -- paired bars (no-kcat vs kcat) at each
     translation efficiency value.

  4. Average dry mass  -- mean dry mass (fg) averaged over all timesteps and
     seeds for each variant, using only the last 8 generations.
"""

import csv
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_sweep import (
	TRL_EFF_VALUES,
	EXPRESSION_FACTOR,
	N_TRL_EFF,
	KCAT_HALF_START,
)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
import os

# Cells longer than this are considered timed-out and excluded from doubling
# time statistics.
MAX_DOUBLING_TIME_MIN = 300


def _variant_label(index):
	"""Return a short human-readable label for a variant index."""
	is_kcat = (index >= KCAT_HALF_START)
	local_index = index - KCAT_HALF_START if is_kcat else index
	kcat_tag = '+kcat' if is_kcat else ''

	if local_index == 0:
		return f'Control\n(KO{", " + kcat_tag if kcat_tag else ""})'

	trl_eff = TRL_EFF_VALUES[local_index - 1]
	return f'Trl {trl_eff}\n{kcat_tag}' if kcat_tag else f'Trl {trl_eff}'


def _variant_color(index):
	"""Return a color for a variant index.
	Blue hues for no-kcat half, red/orange hues for kcat half."""
	is_kcat = (index >= KCAT_HALF_START)
	local_index = index - KCAT_HALF_START if is_kcat else index

	if local_index == 0:
		return '#888888' if not is_kcat else '#555555'

	shade = 0.3 + 0.6 * (local_index / N_TRL_EFF)
	if is_kcat:
		return cm.Reds(shade)
	else:
		return cm.Blues(shade)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variant_indexes = self.ap.get_variants()
		n_variants = len(variant_indexes)
		n_total_gens = self.ap.n_generation

		# Use only the last 8 generations for average doubling time and
		# average dry mass panels.
		ignore_first_n_gens = max(n_total_gens - 8, 0)

		# ------------------------------------------------------------------ #
		# Collect per-variant statistics                                       #
		# ------------------------------------------------------------------ #
		completion_rates = {}            # variant -> float [0, 1]
		dt_by_gen = {}                   # variant -> array (n_gens,) of means
		avg_doubling_time = {}           # variant -> float (mean dt, min)
		avg_dry_mass = {}                # variant -> float (mean dry mass, fg)
		lineage_dts = []                 # list of dicts for per-lineage CSV

		for vi in variant_indexes:
			# -- completion rate ------------------------------------------
			n_gen0 = len(self.ap.get_cells(
				variant=[vi], generation=[0], only_successful=True))
			n_last = len(self.ap.get_cells(
				variant=[vi], generation=[n_total_gens - 1],
				only_successful=True))
			completion_rates[vi] = (n_last / n_gen0) if n_gen0 > 0 else 0.0

			# -- doubling time by generation --------------------------------
			seed_dts = {}  # seed -> {gen: dt_min}
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
						t = TableReader(os.path.join(sim_out, 'Main')
							).readColumn('time', squeeze=True)
						dt_min = (t[-1] - t[0]) / 60.
					except Exception:
						continue
					parts = cell_path.rstrip('/').split('/')
					seed = parts[-3]
					if seed not in seed_dts:
						seed_dts[seed] = {}
					seed_dts[seed][gen] = dt_min
					if dt_min <= MAX_DOUBLING_TIME_MIN:
						times.append(dt_min)
				if times:
					mean_dts[gen] = np.mean(times)
			dt_by_gen[vi] = mean_dts

			# -- average doubling time (skip early gens) --------------------
			later_dts = mean_dts[ignore_first_n_gens:]
			avg_doubling_time[vi] = (
				np.nanmean(later_dts) if np.any(np.isfinite(later_dts))
				else np.nan)

			# Build per-lineage rows for CSV
			for seed in sorted(seed_dts.keys()):
				row = {'variant': vi, 'seed': seed}
				for gen in range(n_total_gens):
					row[f'gen_{gen}_dt_min'] = seed_dts[seed].get(gen, np.nan)
				lineage_dts.append(row)

			# -- average dry mass ------------------------------------------
			all_cells = self.ap.get_cells(
				variant=[vi],
				generation=np.arange(ignore_first_n_gens, n_total_gens),
				only_successful=True)
			cell_means = []
			for cell_path in all_cells:
				sim_out = os.path.join(cell_path, 'simOut')
				try:
					dm = TableReader(os.path.join(sim_out, 'Mass')
						).readColumn('dryMass', squeeze=True)
					cell_means.append(np.mean(dm))
				except Exception:
					pass
			avg_dry_mass[vi] = np.mean(cell_means) if cell_means else np.nan

		# ------------------------------------------------------------------ #
		# Write per-lineage doubling times to CSV                              #
		# ------------------------------------------------------------------ #
		if lineage_dts:
			gen_cols = [f'gen_{g}_dt_min' for g in range(n_total_gens)]
			csv_path = os.path.join(plotOutDir, 'doubling_times_by_lineage.csv')
			with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
				writer = csv.DictWriter(
					fh, fieldnames=['variant', 'seed'] + gen_cols)
				writer.writeheader()
				writer.writerows(lineage_dts)
			print(f'Wrote {csv_path}')

		# ------------------------------------------------------------------ #
		# Plot                                                                 #
		# ------------------------------------------------------------------ #
		fig = plt.figure(figsize=(max(14, n_variants * 0.9), 20))
		gs = fig.add_gridspec(4, 1, hspace=0.45,
							  height_ratios=[1, 2, 1, 1])

		x = np.arange(n_variants)
		labels = [_variant_label(vi) for vi in variant_indexes]
		colors = [_variant_color(vi) for vi in variant_indexes]

		# ---- Panel 1: completion rate ------------------------------------ #
		ax1 = fig.add_subplot(gs[0])
		bars = ax1.bar(x, [completion_rates[vi] for vi in variant_indexes],
					   color=colors, edgecolor='white', linewidth=0.5)
		ax1.axhline(1.0, color='k', lw=0.8, ls='--')
		ax1.set_ylim(0, 1.12)
		ax1.set_ylabel('Completion rate', fontsize=10)
		ax1.set_title('Fraction of seeds completing all generations', fontsize=11)
		ax1.set_xticks(x)
		ax1.set_xticklabels(labels, fontsize=7, rotation=0)
		ax1.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
		for bar, vi in zip(bars, variant_indexes):
			rate = completion_rates[vi]
			ax1.text(bar.get_x() + bar.get_width() / 2,
					 rate + 0.02, f'{rate:.0%}',
					 ha='center', va='bottom', fontsize=7)

		# ---- Panel 2:  mean doubling time heatmap ----------------------- #
		ax2 = fig.add_subplot(gs[1])
		heatmap_data = np.array([dt_by_gen[vi] for vi in variant_indexes])
		vmin = np.nanmin(heatmap_data)
		vmax = min(np.nanmax(heatmap_data), MAX_DOUBLING_TIME_MIN)
		im = ax2.imshow(
			heatmap_data,
			aspect='auto', origin='upper',
			cmap='RdYlGn_r',
			vmin=vmin, vmax=vmax,
			interpolation='nearest')
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

		# ---- Panel 3: average doubling time ------------------------------- #
		ax3 = fig.add_subplot(gs[2])
		dts = [avg_doubling_time[vi] for vi in variant_indexes]
		ax3.bar(x, dts, color=colors, edgecolor='white', linewidth=0.5)
		finite_dts = [d for d in dts if np.isfinite(d)]
		if finite_dts:
			ymin = min(finite_dts) * 0.97
			ymax = max(finite_dts) * 1.03
			ax3.set_ylim(ymin, ymax)
		ax3.set_ylabel('Average doubling time (min)', fontsize=10)
		ax3.set_title(
			f'Mean doubling time averaged over all seeds\n'
			f'(gens {ignore_first_n_gens}\u2013{n_total_gens - 1}; '
			f'cells > {MAX_DOUBLING_TIME_MIN} min excluded)',
			fontsize=11)
		ax3.set_xticks(x)
		ax3.set_xticklabels(labels, fontsize=7, rotation=0)

		# ---- Panel 4: average dry mass ----------------------------------- #
		ax4 = fig.add_subplot(gs[3])
		masses = [avg_dry_mass[vi] for vi in variant_indexes]
		ax4.bar(x, masses, color=colors, edgecolor='white', linewidth=0.5)
		finite_masses = [m for m in masses if np.isfinite(m)]
		if finite_masses:
			ymin = min(finite_masses) * 0.97
			ymax = max(finite_masses) * 1.03
			ax4.set_ylim(ymin, ymax)
		ax4.set_ylabel('Average dry mass (fg)', fontsize=10)
		ax4.set_title(
			f'Mean dry mass averaged over all timesteps and seeds\n'
			f'(gens {ignore_first_n_gens}\u2013{n_total_gens - 1})',
			fontsize=11)
		ax4.set_xticks(x)
		ax4.set_xticklabels(labels, fontsize=7, rotation=0)

		# Legend
		from matplotlib.patches import Patch
		legend_elements = [
			Patch(facecolor='#888888', label='Control (KO, no kcat)'),
			Patch(facecolor='#555555', label='Control (KO, +kcat)'),
			Patch(facecolor=cm.Blues(0.6), label=f'No kcat (exp 10^{EXPRESSION_FACTOR - 1})'),
			Patch(facecolor=cm.Reds(0.6), label=f'+kcat (exp 10^{EXPRESSION_FACTOR - 1})'),
		]
		fig.legend(handles=legend_elements, loc='lower center',
				   ncol=len(legend_elements), fontsize=8,
				   bbox_to_anchor=(0.5, -0.01))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
