"""
Time-series analysis of kcat upper bound saturation across generations.

For each kcat-constrained reaction, tracks how the saturation ratio
(flux / bound) evolves over time across generations.  Designed to diagnose
which constraints are causing lineage failure by revealing reactions that
are persistently at their bound and getting worse.

Outputs:
  1. A ranked summary table printed to stdout (mean saturation, fraction of
     timesteps at >90% and >99% saturation, and generation trend).
  2. A multi-page PDF with two-panel plots for the top N most-saturated
     reactions showing flux vs bound and saturation ratio over time.

Uses only_successful=False to include the last cell before failure.
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

# Number of top-saturated reactions to plot in the PDF.
N_REACTIONS_TO_SHOW = 50

# Ignore first N generations (initialization transient).
IGNORE_FIRST_N_GENS = 4


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		if not hasattr(metabolism, 'selected_kcat_estimates'):
			print('Skipping: sim_data does not have selected_kcat_estimates.'
				  ' Run with the kcat_estimate_scale variant.')
			return

		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=False)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Build per-reaction lookup ------------------------------------
		fba_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		# rxn_to_pairs: rxn_id -> [(cat_listener_idx, kcat), ...]
		rxn_to_pairs = {}
		for (rxn_id, cat_id), kcat in metabolism.selected_kcat_estimates.items():
			if rxn_id not in rxn_id_to_idx or cat_id not in cat_id_to_idx:
				continue
			rxn_to_pairs.setdefault(rxn_id, []).append(
				(cat_id_to_idx[cat_id], kcat))

		if not rxn_to_pairs:
			print('Skipping: no kcat-constrained reactions found in listener IDs.')
			return

		constrained_rxn_ids = list(rxn_to_pairs.keys())
		rxn_indexes = np.array([rxn_id_to_idx[r] for r in constrained_rxn_ids])

		all_cat_indexes = sorted({ci for pairs in rxn_to_pairs.values()
								  for ci, _ in pairs})
		cat_idx_local = {ci: li for li, ci in enumerate(all_cat_indexes)}

		cell_density = sim_data.constants.cell_density
		n_rxns = len(constrained_rxn_ids)

		# --- Accumulate per-cell time traces ------------------------------
		# For the summary we collect ratios per cell; for the PDF we need
		# time, flux, bound arrays concatenated across cells.
		time_all = []           # hours
		flux_all = []           # (T, n_rxns)
		bound_all = []          # (T, n_rxns)
		gen_boundary_times = [] # hours
		cell_generation = []    # generation index per cell

		# Per-reaction ratio accumulator for summary stats
		ratio_accumulator = [[] for _ in range(n_rxns)]
		# Per-reaction per-generation mean ratio for trend computation
		gen_ratio_means = [[] for _ in range(n_rxns)]  # list of (gen, mean)

		time_offset = 0.0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			gen_idx = self.ap.get_cell_generation(cell_path)

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
							 indices=np.array(all_cat_indexes),
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
			cell_generation.append(gen_idx)

			# Reaction fluxes -> mmol/g DCW/h
			conversion_coeffs = (dry_mass / cell_mass
								 * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			# Bounds -> mmol/g DCW/h: max over catalysts of kcat * [E]
			bounds = np.zeros_like(fluxes)
			for ri, rxn_id in enumerate(constrained_rxn_ids):
				for cat_idx, kcat in rxn_to_pairs[rxn_id]:
					li = cat_idx_local[cat_idx]
					conc = cat_counts[:, li] * counts_to_molar
					bounds[:, ri] = np.maximum(bounds[:, ri], kcat * conc)

			time_all.append(t_h)
			flux_all.append(fluxes)
			bound_all.append(bounds)

			# Collect ratios for summary stats
			for ri in range(n_rxns):
				valid = (bounds[:, ri] > 0) & np.isfinite(fluxes[:, ri])
				if valid.any():
					ratios = fluxes[valid, ri] / bounds[valid, ri]
					ratio_accumulator[ri].append(ratios)
					gen_ratio_means[ri].append((gen_idx, np.mean(ratios)))

		if not time_all:
			print('No cells processed successfully.')
			return

		# --- Compute per-reaction summary statistics ----------------------
		mean_sat = np.full(n_rxns, np.nan)
		frac_90 = np.full(n_rxns, np.nan)
		frac_99 = np.full(n_rxns, np.nan)
		trend_slope = np.full(n_rxns, np.nan)

		for ri in range(n_rxns):
			if not ratio_accumulator[ri]:
				continue
			ratios = np.concatenate(ratio_accumulator[ri])
			mean_sat[ri] = np.mean(ratios)
			frac_90[ri] = np.mean(ratios > 0.9)
			frac_99[ri] = np.mean(ratios > 0.99)

			# Trend: simple linear regression of per-generation mean ratio
			# vs generation index
			if len(gen_ratio_means[ri]) >= 2:
				gens_arr = np.array([g for g, _ in gen_ratio_means[ri]])
				means_arr = np.array([m for _, m in gen_ratio_means[ri]])
				if len(np.unique(gens_arr)) >= 2:
					coeffs = np.polyfit(gens_arr, means_arr, 1)
					trend_slope[ri] = coeffs[0]

		# Sort by mean saturation descending
		has_data = ~np.isnan(mean_sat)
		sort_order = np.argsort(mean_sat[has_data])[::-1]
		ranked_idxs = np.where(has_data)[0][sort_order]

		# --- Print summary table ------------------------------------------
		print('\n' + '=' * 100)
		print(f'{"Rank":>4}  {"Reaction ID":<55}  {"Mean Sat":>8}  '
			  f'{">90%":>6}  {">99%":>6}  {"Trend":>8}')
		print('-' * 100)
		for rank, ri in enumerate(ranked_idxs):
			trend_str = (f'{trend_slope[ri]:+.4f}'
						 if np.isfinite(trend_slope[ri]) else '   N/A')
			print(f'{rank + 1:4d}  {constrained_rxn_ids[ri]:<55}  '
				  f'{mean_sat[ri]:8.4f}  {frac_90[ri]:6.2%}  '
				  f'{frac_99[ri]:6.2%}  {trend_str:>8}')
		print('=' * 100 + '\n')

		# --- Write CSV summary --------------------------------------------
		csv_path = os.path.join(plotOutDir, plotOutFileName + '_summary.csv')
		with open(csv_path, 'w') as f:
			f.write('rank,reaction_id,mean_saturation,frac_above_90pct,'
					'frac_above_99pct,trend_slope\n')
			for rank, ri in enumerate(ranked_idxs):
				f.write(f'{rank + 1},{constrained_rxn_ids[ri]},'
						f'{mean_sat[ri]:.6f},{frac_90[ri]:.6f},'
						f'{frac_99[ri]:.6f},'
						f'{trend_slope[ri]:.6f}\n')
		print(f'Summary CSV written to {csv_path}')

		# --- Multi-page PDF of top saturated reactions --------------------
		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)
		bound_cat = np.vstack(bound_all)

		top_idxs = ranked_idxs[:N_REACTIONS_TO_SHOW]

		quantile = getattr(metabolism, 'kcat_estimate_quantile', 'unknown')
		multiplier = getattr(metabolism, 'kcat_estimate_multiplier', '?')

		pdf_path = os.path.join(plotOutDir, plotOutFileName + '.pdf')
		with PdfPages(pdf_path) as pdf:
			for ri in top_idxs:
				rxn_id = constrained_rxn_ids[ri]
				fig, (ax_top, ax_bot) = plt.subplots(
					2, 1, figsize=(12, 6), sharex=True)

				fig.suptitle(
					f'{rxn_id}\n'
					f'quantile: {quantile}, multiplier: {multiplier}  |  '
					f'mean sat: {mean_sat[ri]:.3f},  '
					f'>90%: {frac_90[ri]:.1%},  >99%: {frac_99[ri]:.1%}',
					fontsize=10)

				# Top panel: flux and bound
				ax_top.plot(time_cat, flux_cat[:, ri],
							lw=0.6, color='#1f77b4', alpha=0.85, label='flux')
				ax_top.plot(time_cat, bound_cat[:, ri],
							lw=0.6, color='#d62728', alpha=0.85, label='bound')
				ax_top.set_ylabel('mmol/g DCW/h', fontsize=9)
				ax_top.legend(fontsize=8, loc='upper right')

				# Bottom panel: saturation ratio
				with np.errstate(divide='ignore', invalid='ignore'):
					ratio = flux_cat[:, ri] / bound_cat[:, ri]
				ratio[~np.isfinite(ratio)] = np.nan
				ax_bot.plot(time_cat, ratio,
							lw=0.6, color='#2ca02c', alpha=0.85)
				ax_bot.axhline(1.0, color='k', lw=0.8, ls='--', alpha=0.5)
				ax_bot.set_ylabel('flux / bound', fontsize=9)
				ax_bot.set_xlabel('Time (h)', fontsize=9)
				ax_bot.set_ylim(-0.05, max(1.2, np.nanmax(ratio) * 1.1)
								if np.any(np.isfinite(ratio)) else 1.2)

				# Generation boundaries on both panels
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


if __name__ == '__main__':
	Plot().cli()
