"""
Per-reaction breakdown of kcat upper bound saturation across multiple
generations.

For each kcat-constrained reaction, computes two statistics over all
timesteps where the bound is positive (enzyme present):

  - Saturation rate: fraction of timesteps where flux > SATURATION_THRESHOLD
    of the bound (i.e. the bound is actively limiting).
  - Violation rate: fraction of timesteps where flux exceeds the bound
    (flux / bound > 1.0), which should be zero if constraints are working.

Produces two horizontal bar charts, ranked by saturation rate:

  Left:  Saturation rate per reaction (fraction of timesteps with
         flux > SATURATION_THRESHOLD * bound).
  Right: Mean flux/bound ratio per reaction (with std error bars).

Reactions are labeled as "reaction_id / catalyst_id".  Only the top
N_REACTIONS_TO_SHOW reactions by saturation rate are displayed.

Skips gracefully if the simulation was not run with the kcat_estimate_scale
variant (i.e., sim_data.process.metabolism.selected_kcat_estimates is absent).
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

IGNORE_FIRST_N_GENS = 4

# Bound is considered "active" when flux exceeds this fraction of the bound.
SATURATION_THRESHOLD = 0.9

# Show only the top N reactions by saturation rate to keep the plot readable.
N_REACTIONS_TO_SHOW = 40


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
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Build per-reaction lookup ------------------------------------
		fba_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
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

		# Accumulators: per reaction, collect all flux/bound ratios
		ratio_accumulator = [[] for _ in range(n_rxns)]  # list of arrays

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
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
							 squeeze=False)[1:]   # (T, n_cats_needed)

				rxn_fluxes_raw = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('reactionFluxes',
							 indices=rxn_indexes,
							 squeeze=False)[1:]   # (T, n_rxns)

			except Exception as e:
				print(f'Ignored exception reading {sim_out}: {e!r}')
				continue

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

			# Collect flux/bound ratios for timesteps where bound > 0
			for ri in range(n_rxns):
				valid = (bounds[:, ri] > 0) & np.isfinite(fluxes[:, ri])
				if valid.any():
					ratio_accumulator[ri].append(
						fluxes[valid, ri] / bounds[valid, ri])

		# --- Compute per-reaction statistics ------------------------------
		saturation_rates = np.full(n_rxns, np.nan)
		violation_rates = np.full(n_rxns, np.nan)
		mean_ratios = np.full(n_rxns, np.nan)
		std_ratios = np.full(n_rxns, np.nan)
		n_timesteps = np.zeros(n_rxns, dtype=int)

		for ri in range(n_rxns):
			if not ratio_accumulator[ri]:
				continue
			ratios = np.concatenate(ratio_accumulator[ri])
			n_timesteps[ri] = len(ratios)
			saturation_rates[ri] = np.mean(ratios > SATURATION_THRESHOLD)
			violation_rates[ri] = np.mean(ratios > 1.0)
			mean_ratios[ri] = np.mean(ratios)
			std_ratios[ri] = np.std(ratios) / np.sqrt(len(ratios))

		# Sort by saturation rate, keep top N with any data
		has_data = ~np.isnan(saturation_rates)
		sort_order = np.argsort(saturation_rates[has_data])[::-1]
		ranked_idxs = np.where(has_data)[0][sort_order]
		ranked_idxs = ranked_idxs[:N_REACTIONS_TO_SHOW]

		# --- Plot ---------------------------------------------------------
		n_show = len(ranked_idxs)
		fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, max(6, n_show * 0.35)))

		y = np.arange(n_show)
		labels = [
			f'{constrained_rxn_ids[ri]}'
			for ri in ranked_idxs
		]

		quantile = getattr(metabolism, 'kcat_estimate_quantile', 'unknown')
		multiplier = getattr(metabolism, 'kcat_estimate_multiplier', '?')
		fig.suptitle(
			f'kcat bound saturation by reaction  —  quantile: {quantile}, '
			f'multiplier: {multiplier}  —  gens {IGNORE_FIRST_N_GENS}+  '
			f'(saturation threshold: flux > {SATURATION_THRESHOLD} × bound)',
			fontsize=11)

		# Left panel: saturation rate bars, colored by violation rate
		sat = saturation_rates[ranked_idxs]
		viol = violation_rates[ranked_idxs]
		bar_colors = ['#d62728' if v > 0 else '#1f77b4' for v in viol]
		ax1.barh(y, sat, color=bar_colors, edgecolor='none', height=0.7)
		ax1.axvline(0, color='k', lw=0.5)
		ax1.set_xlim(0, 1.05)
		ax1.set_xlabel(f'Fraction of timesteps with flux > {SATURATION_THRESHOLD} × bound',
					   fontsize=9)
		ax1.set_yticks(y)
		ax1.set_yticklabels(labels, fontsize=6)
		ax1.invert_yaxis()
		ax1.set_title('Bound saturation rate\n(red = has violations > 1.0)', fontsize=10)
		ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0%}'))

		# Right panel: mean flux/bound ratio with std error bars
		mean_r = mean_ratios[ranked_idxs]
		std_r = std_ratios[ranked_idxs]
		n_ts = n_timesteps[ranked_idxs]
		ax2.barh(y, mean_r, xerr=std_r, color='#ff7f0e',
				 edgecolor='none', height=0.7, error_kw=dict(lw=0.8, capsize=2))
		ax2.axvline(1.0, color='k', lw=1, ls='--', label='flux = bound')
		ax2.set_xlabel('Mean flux / bound  (± SE)', fontsize=9)
		ax2.set_yticks(y)
		ax2.set_yticklabels([], fontsize=6)
		ax2.invert_yaxis()
		ax2.set_title('Mean flux/bound ratio\n(dashed = bound exactly met)', fontsize=10)
		ax2.legend(fontsize=8)

		# Annotate n timesteps on right panel
		for yi, (ri, n) in enumerate(zip(ranked_idxs, n_ts)):
			ax2.text(ax2.get_xlim()[1] * 0.99, yi,
					 f'n={n:,}', va='center', ha='right', fontsize=5, color='gray')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
