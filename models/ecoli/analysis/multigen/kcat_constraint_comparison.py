"""
Compare actual flux against multiple kcat constraint levels for buffered
reactions and high-saturation reactions from variant 2 seeds 3/5.

For each reaction of interest, produces a two-panel page:
  Top:    Actual flux (blue) vs five constraint lines (max, smoothed_max,
          smoothed_max*1.1, smoothed_max_buffered, p999).
  Bottom: Raw enzyme (catalyst) counts over time.

Designed to run on wildtype sims where no kcat constraints are active but
all flux/enzyme data is available.

Output: multi-page PDF + summary to stdout.
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

# Only include buffered reactions with kcat >= this threshold.
MIN_KCAT_ESTIMATE = 1.0

# Extra reactions to include (high saturation in variant 2 seeds 3/5).
EXTRA_REACTIONS = [
	('PNPOXI-RXN', 'PDXH-CPLX[c]'),
	('CTPSYN-RXN', 'CTPSYN-CPLX[c]'),
	('L-ASPARTATE-OXID-RXN', 'L-ASPARTATE-OXID-MONOMER[c]'),
	('CHORISMATEMUT-RXN__CHORISMUTPREPHENDEHYDROG-CPLX', 'CHORISMUTPREPHENDEHYDROG-CPLX[c]'),
	('ASPCARBTRANS-RXN', 'ASPCARBTRANS-CPLX[c]'),
	('GLYOHMETRANS-RXN-SER/THF//GLY/METHYLENE-THF/WATER.33.', 'GLYOHMETRANS-CPLX[c]'),
	('DXPREDISOM-RXN', 'DXPREDISOM-CPLX[i]'),
	]

# Quantiles to plot, with display properties.
# Each entry: (quantile_key, label_template, color, linestyle)
QUANTILE_STYLES = [
	('max',                  'max (kcat={:.2f})',         'green',  '-'),
	('smoothed_max',         'smoothed_max (kcat={:.2f})','red',    '-'),
	('smoothed_max_x1.1',   'SM x 1.1 (kcat={:.2f})',   'orange', '--'),
	('smoothed_max_buffered','SM buffered (kcat={:.2f})', 'purple', '-'),
	('p999',                 'p999 (kcat={:.2f})',        'brown',  '-'),
	]


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism

		# --- Build target reaction set ------------------------------------
		# 1) Buffered reactions with kcat >= MIN_KCAT_ESTIMATE
		targets = set()
		sm_buffered = metabolism.kcat_estimates.get('smoothed_max_buffered', {})
		for pair in metabolism.kcat_buffered_reactions:
			kcat_val = sm_buffered.get(pair)
			if kcat_val is not None and kcat_val >= MIN_KCAT_ESTIMATE:
				targets.add(pair)

		# 2) Extra unbuffered reactions
		for pair in EXTRA_REACTIONS:
			targets.add(pair)

		if not targets:
			print('No target reactions found.')
			return

		# --- Collect kcat values per quantile per target ------------------
		# target_kcats: {(rxn, cat): {quantile_key: kcat_value}}
		target_kcats = {}
		for pair in sorted(targets):
			kcats = {}
			for qkey, _, _, _ in QUANTILE_STYLES:
				if qkey == 'smoothed_max_x1.1':
					val = metabolism.kcat_estimates.get('smoothed_max', {}).get(pair)
					if val is not None:
						kcats[qkey] = val * 1.1
				else:
					val = metabolism.kcat_estimates.get(qkey, {}).get(pair)
					if val is not None:
						kcats[qkey] = val
			target_kcats[pair] = kcats

		# --- Get cell paths -----------------------------------------------
		cell_paths = self.ap.get_cells(only_successful=True)
		if len(cell_paths) == 0:
			print('Skipping analysis -- no successful cells found.')
			return

		# --- Build index lookups from first cell --------------------------
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		# Filter targets to those present in listener
		valid_targets = []
		for rxn_id, cat_id in sorted(targets):
			if rxn_id in rxn_id_to_idx and cat_id in cat_id_to_idx:
				valid_targets.append((rxn_id, cat_id))
			else:
				print(f'Warning: {rxn_id} / {cat_id} not found in listener, skipping.')

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
		# bounds_all: {quantile_key: [arrays]}
		bounds_all = {qkey: [] for qkey, _, _, _ in QUANTILE_STYLES}
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

			# Time in hours, offset-corrected
			t_sec = t_raw - t_raw[0]
			t_h = (t_sec + time_offset) / 3600.
			gen_boundary_times.append(time_offset / 3600.)
			time_offset += (t_raw[-1] - t_raw[0])

			# Reaction fluxes -> mmol/g DCW/h
			conversion_coeffs = (dry_mass / cell_mass
								 * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			# Bounds for each quantile: kcat * counts * countsToMolar
			# -> same units as flux (mmol/g DCW/h)
			for ti in range(n_targets):
				pair = valid_targets[ti]
				conc = cat_counts[:, ti] * counts_to_molar  # molar
				for qkey, _, _, _ in QUANTILE_STYLES:
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

		# Concatenate across cells
		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)
		cat_counts_cat = np.vstack(cat_counts_all)
		bounds_cat = {}
		for qkey in bounds_all:
			if bounds_all[qkey]:
				bounds_cat[qkey] = np.vstack(bounds_all[qkey])

		# --- Compute sorting metric: mean saturation vs smoothed_max ------
		mean_sat = np.full(n_targets, np.nan)
		for ti in range(n_targets):
			if 'smoothed_max' in bounds_cat:
				bound = bounds_cat['smoothed_max'][:, ti]
				valid = (bound > 0) & np.isfinite(flux_cat[:, ti])
				if valid.any():
					mean_sat[ti] = np.nanmean(
						flux_cat[valid, ti] / bound[valid])

		# Sort by saturation descending (NaN last)
		sort_order = np.argsort(np.where(np.isnan(mean_sat), -1, mean_sat))[::-1]

		# --- Print summary ------------------------------------------------
		print('\n' + '=' * 110)
		print(f'{"Rank":>4}  {"Reaction":<40}  {"Catalyst":<25}  '
			  f'{"Buffered":>8}  {"Mean Sat":>8}')
		print('-' * 110)
		for rank, ti in enumerate(sort_order):
			rxn_id, cat_id = valid_targets[ti]
			is_buffered = (rxn_id, cat_id) in metabolism.kcat_buffered_reactions
			sat_str = f'{mean_sat[ti]:.4f}' if np.isfinite(mean_sat[ti]) else 'N/A'
			print(f'{rank + 1:4d}  {rxn_id:<40}  {cat_id:<25}  '
				  f'{"yes" if is_buffered else "no":>8}  {sat_str:>8}')
		print('=' * 110 + '\n')

		# --- Multi-page PDF -----------------------------------------------
		pdf_path = os.path.join(plotOutDir, plotOutFileName + '.pdf')
		with PdfPages(pdf_path) as pdf:
			for ti in sort_order:
				rxn_id, cat_id = valid_targets[ti]
				is_buffered = (rxn_id, cat_id) in metabolism.kcat_buffered_reactions
				kcats = target_kcats[(rxn_id, cat_id)]

				fig, (ax_top, ax_bot) = plt.subplots(
					2, 1, figsize=(14, 8), sharex=True,
					gridspec_kw={'height_ratios': [2, 1]})

				# Title
				kcat_strs = []
				for qkey, _, _, _ in QUANTILE_STYLES:
					if qkey in kcats:
						kcat_strs.append(f'{qkey}={kcats[qkey]:.2f}')
				fig.suptitle(
					f'{rxn_id}  /  {cat_id}\n'
					f'Buffered: {"yes" if is_buffered else "no"}  |  '
					+ '  '.join(kcat_strs),
					fontsize=10)

				# Top panel: flux vs constraint lines
				ax_top.plot(time_cat, flux_cat[:, ti],
							lw=0.8, color='blue', alpha=0.85, label='flux')

				for qkey, label_tmpl, color, ls in QUANTILE_STYLES:
					if qkey in bounds_cat and qkey in kcats:
						ax_top.plot(
							time_cat, bounds_cat[qkey][:, ti],
							lw=0.8, color=color, ls=ls, alpha=0.8,
							label=label_tmpl.format(kcats[qkey]))

				ax_top.set_ylabel('mmol/g DCW/h', fontsize=9)
				ax_top.legend(fontsize=7, loc='upper right', ncol=2)

				# Bottom panel: enzyme counts
				ax_bot.plot(time_cat, cat_counts_cat[:, ti],
							lw=0.8, color='teal', alpha=0.85)
				ax_bot.set_ylabel('Enzyme counts (molecules)', fontsize=9)
				ax_bot.set_xlabel('Time (h)', fontsize=9)

				# Generation boundaries
				for gt in gen_boundary_times[1:]:
					ax_top.axvline(gt, color='gray', lw=0.5, ls='--', alpha=0.5)
					ax_bot.axvline(gt, color='gray', lw=0.5, ls='--', alpha=0.5)

				for ax in (ax_top, ax_bot):
					ax.set_xlim(time_cat[0], time_cat[-1])
					ax.tick_params(labelsize=7)

				plt.tight_layout()
				pdf.savefig(fig)
				plt.close(fig)

		print(f'PDF written to {pdf_path}')
		print(f'Total pages: {len(valid_targets)}')


if __name__ == '__main__':
	Plot().cli()
