"""
Compare actual flux against multiple kcat constraint levels, grouped by
behavioral category.

Covers ALL reactions with smoothed_max kcat >= 1.0 (not just buffered).
Classifies each into behavioral profiles (always saturated, spiking,
sparse, normal) and shows top N per category with multi-quantile
constraint lines.

Each reaction page has two panels:
  Top:    Actual flux (blue) vs five constraint lines (max, smoothed_max,
          smoothed_max*1.1, smoothed_max_buffered, p999).
  Bottom: Raw enzyme (catalyst) counts over time.

Designed to run on wildtype sims where no kcat constraints are active but
all flux/enzyme data is available.

Output: multi-page PDF (with category header pages) + summary to stdout.
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Only include reactions with kcat >= this threshold.
MIN_KCAT_ESTIMATE = 1.0

# Number of representative reactions per behavioral category.
N_PER_CATEGORY = 30

# Classification thresholds (same as kcat_counterfactual_bounds).
SPIKE_RATIO = 5.0
SATURATED_THRESHOLD = 0.8
SPARSE_ZERO_FRACTION = 0.5

# Ignore first N generations (initialization transient).
IGNORE_FIRST_N_GENS = 4

# Quantiles to plot, with display properties.
# Each entry: (quantile_key, label_template, color, dashes_or_linestyle, linewidth)
QUANTILE_STYLES = [
	('max',                  'max (kcat={:.2f})',         'green',  '-',           1.8),
	('smoothed_max',         'smoothed_max (kcat={:.2f})','red',    '-',           1.4),
	('smoothed_max_x1.1',   'SM x 1.1 (kcat={:.2f})',   'orange', (5, 3),        1.4),
	('smoothed_max_buffered','SM buffered (kcat={:.2f})', 'purple', (8, 3, 2, 3), 1.2),
	('p999',                 'p999 (kcat={:.2f})',        'brown',  (2, 2),        1.0),
	]

CATEGORY_ORDER = ['always_saturated', 'spiking', 'sparse', 'normal']
CATEGORY_LABELS = {
	'always_saturated': 'Always Saturated (mean flux/bound > {:.1f})'.format(SATURATED_THRESHOLD),
	'spiking': 'Spiking (max/median > {:.0f}x)'.format(SPIKE_RATIO),
	'sparse': 'Sparse (>{:.0%} zero-flux timesteps)'.format(SPARSE_ZERO_FRACTION),
	'normal': 'Normal',
}


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism

		# --- Build target reaction set ------------------------------------
		# All (rxn_id, cat_id) pairs from smoothed_max with kcat >= threshold.
		sm_estimates = metabolism.kcat_estimates.get('smoothed_max', {})
		targets = {
			pair for pair, kcat in sm_estimates.items()
			if kcat >= MIN_KCAT_ESTIMATE
		}

		if not targets:
			print('No target reactions found.')
			return

		# --- Collect kcat values per quantile per target ------------------
		target_kcats = {}
		for pair in sorted(targets):
			kcats = {}
			for qkey, _, _, _, _ in QUANTILE_STYLES:
				if qkey == 'smoothed_max_x1.1':
					val = sm_estimates.get(pair)
					if val is not None:
						kcats[qkey] = val * 1.1
				else:
					val = metabolism.kcat_estimates.get(qkey, {}).get(pair)
					if val is not None:
						kcats[qkey] = val
			target_kcats[pair] = kcats

		# --- Get cell paths -----------------------------------------------
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Build index lookups from first cell --------------------------
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		valid_targets = []
		for rxn_id, cat_id in sorted(targets):
			if rxn_id in rxn_id_to_idx and cat_id in cat_id_to_idx:
				valid_targets.append((rxn_id, cat_id))

		if not valid_targets:
			print('No valid target reactions found in listener.')
			return

		n_targets = len(valid_targets)
		rxn_indexes = np.array([rxn_id_to_idx[r] for r, _ in valid_targets])
		cat_indexes = np.array([cat_id_to_idx[c] for _, c in valid_targets])

		cell_density = sim_data.constants.cell_density

		# --- Accumulate time traces per cell ------------------------------
		time_all = []
		flux_all = []
		cat_counts_all = []
		bounds_all = {qkey: [] for qkey, _, _, _, _ in QUANTILE_STYLES}
		gen_boundary_times = []

		time_offset = 0.0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')

			try:
				t_raw = TableReader(
					os.path.join(sim_out, 'Main')
				).readColumn('time', squeeze=True)[1:]

				counts_to_molar = TableReader(
					os.path.join(sim_out, 'EnzymeKinetics')
				).readColumn('countsToMolar', squeeze=True)[1:]

				cell_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('cellMass', squeeze=True)[1:]

				dry_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('dryMass', squeeze=True)[1:]

				cat_counts = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('catalyst_counts',
							 indices=cat_indexes,
							 squeeze=False)[1:]

				rxn_fluxes_raw = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('reactionFluxes',
							 indices=rxn_indexes,
							 squeeze=False)[1:]

			except Exception as e:
				print(f'Ignored exception reading {sim_out}: {e!r}')
				continue

			t_sec = t_raw - t_raw[0]
			t_h = (t_sec + time_offset) / 3600.
			gen_boundary_times.append(time_offset / 3600.)
			time_offset += (t_raw[-1] - t_raw[0])

			conversion_coeffs = (dry_mass / cell_mass
								 * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			for ti in range(n_targets):
				pair = valid_targets[ti]
				conc = cat_counts[:, ti] * counts_to_molar
				for qkey, _, _, _, _ in QUANTILE_STYLES:
					kcat_val = target_kcats[pair].get(qkey)
					if kcat_val is not None:
						if len(bounds_all[qkey]) <= len(time_all):
							bounds_all[qkey].append(
								np.full((len(t_h), n_targets), np.nan))
						bounds_all[qkey][-1][:, ti] = kcat_val * conc

			time_all.append(t_h)
			flux_all.append(fluxes)
			cat_counts_all.append(cat_counts)

		if not time_all:
			print('No cells processed successfully.')
			return

		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)
		cat_counts_cat = np.vstack(cat_counts_all)
		bounds_cat = {}
		for qkey in bounds_all:
			if bounds_all[qkey]:
				bounds_cat[qkey] = np.vstack(bounds_all[qkey])

		# --- Classify reactions into behavioral profiles ------------------
		# Classification uses smoothed_max bound (the primary constraint).
		categories = {}  # category -> list of (ti, sort_key)
		for ti in range(n_targets):
			f = flux_cat[:, ti]
			if 'smoothed_max' not in bounds_cat:
				continue
			b = bounds_cat['smoothed_max'][:, ti]
			valid = (b > 0) & np.isfinite(f) & np.isfinite(b)

			if not valid.any():
				continue

			f_valid = f[valid]
			b_valid = b[valid]

			median_f = np.median(np.abs(f_valid))
			max_f = np.max(np.abs(f_valid))
			frac_zero = np.mean(np.abs(f_valid) < 1e-12)

			with np.errstate(divide='ignore', invalid='ignore'):
				ratio = f_valid / b_valid
			ratio = ratio[np.isfinite(ratio)]
			mean_ratio = np.mean(ratio) if len(ratio) > 0 else 0.0

			if median_f > 0 and max_f / median_f > SPIKE_RATIO:
				cat = 'spiking'
				sort_key = max_f / median_f
			elif mean_ratio > SATURATED_THRESHOLD:
				cat = 'always_saturated'
				sort_key = mean_ratio
			elif frac_zero > SPARSE_ZERO_FRACTION:
				cat = 'sparse'
				sort_key = frac_zero
			else:
				cat = 'normal'
				sort_key = mean_ratio

			categories.setdefault(cat, []).append((ti, sort_key))

		# Sort each category by sort_key descending, pick top N
		selected = {}
		for cat, entries in categories.items():
			entries.sort(key=lambda x: -x[1])
			selected[cat] = [ti for ti, _ in entries[:N_PER_CATEGORY]]

		# --- Print summary ------------------------------------------------
		print('\n' + '=' * 100)
		print('Constraint comparison by category')
		print('-' * 100)
		for cat in CATEGORY_ORDER:
			if cat not in categories:
				continue
			total = len(categories[cat])
			shown = len(selected.get(cat, []))
			print(f'  {cat:<20s}: {total:4d} reactions total, showing top {shown}')
			for ti in selected.get(cat, []):
				rxn_id, cat_id = valid_targets[ti]
				is_buf = (rxn_id, cat_id) in metabolism.kcat_buffered_reactions
				print(f'    - {rxn_id:<50s}  {"[buffered]" if is_buf else ""}')
		print('=' * 100 + '\n')

		# --- Multi-page PDF -----------------------------------------------
		pdf_path = os.path.join(plotOutDir, plotOutFileName + '.pdf')
		n_pages = 0
		with PdfPages(pdf_path) as pdf:
			for cat in CATEGORY_ORDER:
				if cat not in selected or not selected[cat]:
					continue

				# Category header page
				fig = plt.figure(figsize=(14, 2))
				fig.text(0.5, 0.5, CATEGORY_LABELS.get(cat, cat),
						 ha='center', va='center', fontsize=18,
						 fontweight='bold')
				fig.text(0.5, 0.25,
						 f'{len(categories[cat])} reactions in this category, '
						 f'showing top {len(selected[cat])}',
						 ha='center', va='center', fontsize=12, color='gray')
				pdf.savefig(fig)
				plt.close(fig)
				n_pages += 1

				for ti in selected[cat]:
					rxn_id, cat_id = valid_targets[ti]
					is_buffered = (rxn_id, cat_id) in metabolism.kcat_buffered_reactions
					kcats = target_kcats[(rxn_id, cat_id)]

					fig, (ax_top, ax_bot) = plt.subplots(
						2, 1, figsize=(14, 8), sharex=True,
						gridspec_kw={'height_ratios': [2, 1]})

					# Title
					kcat_strs = []
					for qkey, _, _, _, _ in QUANTILE_STYLES:
						if qkey in kcats:
							kcat_strs.append(f'{qkey}={kcats[qkey]:.2f}')
					fig.suptitle(
						f'{rxn_id}  /  {cat_id}  [{cat}]\n'
						f'Buffered: {"yes" if is_buffered else "no"}  |  '
						+ '  '.join(kcat_strs),
						fontsize=10)

					# Top panel: constraint lines first, then flux on top
					for qkey, label_tmpl, color, dashes, lw in QUANTILE_STYLES:
						if qkey in bounds_cat and qkey in kcats:
							kwargs = dict(
								lw=lw, color=color, alpha=0.85,
								label=label_tmpl.format(kcats[qkey]))
							if isinstance(dashes, tuple):
								kwargs['dashes'] = dashes
							else:
								kwargs['ls'] = dashes
							ax_top.plot(
								time_cat, bounds_cat[qkey][:, ti], **kwargs)

					ax_top.plot(time_cat, flux_cat[:, ti],
								lw=1.6, color='blue', alpha=0.9, label='flux',
								zorder=10)

					ax_top.set_ylabel('mmol/g DCW/h', fontsize=9)
					ax_top.legend(fontsize=7, loc='upper right', ncol=2)

					# Bottom panel: enzyme counts
					ax_bot.plot(time_cat, cat_counts_cat[:, ti],
								lw=0.8, color='teal', alpha=0.85)
					ax_bot.set_ylabel('Enzyme counts (molecules)', fontsize=9)
					ax_bot.set_xlabel('Time (h)', fontsize=9)

					# Generation boundaries
					for gt in gen_boundary_times[1:]:
						ax_top.axvline(gt, color='gray', lw=0.5, ls='--',
									   alpha=0.5)
						ax_bot.axvline(gt, color='gray', lw=0.5, ls='--',
									   alpha=0.5)

					for ax in (ax_top, ax_bot):
						ax.set_xlim(time_cat[0], time_cat[-1])
						ax.tick_params(labelsize=7)

					plt.tight_layout()
					pdf.savefig(fig)
					plt.close(fig)
					n_pages += 1

		print(f'PDF written to {pdf_path}')
		print(f'Total pages: {n_pages}')


if __name__ == '__main__':
	Plot().cli()
