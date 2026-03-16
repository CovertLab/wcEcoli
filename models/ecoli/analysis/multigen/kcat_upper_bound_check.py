"""
Checks whether kcat-based reaction flux upper bounds are being respected in
simulations run with the kcat_estimate_scale variant.

For each kcat-constrained reaction at every timestep, the per-reaction upper
bound is reconstructed from listener data:

    bound (mmol/g DCW/h) = max over catalysts of (kcat * [E])
                         = max(kcat_i * catalyst_count_i * counts_to_molar)

This is then compared against the actual reaction flux (also in mmol/g DCW/h).

Produces two panels:
  1. Log-log scatter of actual flux vs computed bound with the 1:1 line.
     All points should fall on or below the diagonal.
  2. Histogram of flux / bound ratios (where bound > 0).
     All values should be <= 1.

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

IGNORE_FIRST_N_GENS = 0


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

		# --- Build per-reaction lookup from selected_kcat_estimates ----------
		# Group by reaction: reaction_id -> list of (catalyst_idx, kcat)
		# catalyst_idx indexes into FBAResults/catalyst_counts.
		fba_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids  = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids  = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		# Deduplicate into per-reaction structure; skip unknown IDs
		rxn_to_pairs = {}   # rxn_id -> [(cat_idx, kcat), ...]
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

		# All unique catalyst indexes needed (may be fewer than n_reactions)
		all_cat_indexes = sorted({ci for pairs in rxn_to_pairs.values()
								  for ci, _ in pairs})
		cat_idx_local = {ci: li for li, ci in enumerate(all_cat_indexes)}

		cell_density = sim_data.constants.cell_density

		# --- Accumulate (flux, bound) pairs one cell at a time ---------------
		all_fluxes = []   # mmol/g DCW/h
		all_bounds = []   # mmol/g DCW/h

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				counts_to_molar = TableReader(
					os.path.join(sim_out, 'EnzymeKinetics')
				).readColumn('countsToMolar', squeeze=True)[1:]   # (T,)

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

			# Convert reaction fluxes to mmol/g DCW/h  (same as kcat_estimations)
			conversion_coeffs = (dry_mass / cell_mass
								 * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))  # (T,)
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)   # (T, n_rxns)
			fluxes[np.isinf(fluxes)] = np.nan

			# Compute per-reaction upper bound = max(kcat_i * [E_i]) over catalysts
			# bound in mmol/g DCW/h: kcat (L/g/h) * [E] (mmol/L) = mmol/g/h
			bounds = np.zeros_like(fluxes)   # (T, n_rxns)
			for ri, rxn_id in enumerate(constrained_rxn_ids):
				for cat_idx, kcat in rxn_to_pairs[rxn_id]:
					li = cat_idx_local[cat_idx]
					conc = cat_counts[:, li] * counts_to_molar   # mmol/L  (T,)
					candidate = kcat * conc                       # mmol/g/h (T,)
					bounds[:, ri] = np.maximum(bounds[:, ri], candidate)

			all_fluxes.append(fluxes)
			all_bounds.append(bounds)

		if not all_fluxes:
			print('No cells processed successfully.')
			return

		fluxes = np.vstack(all_fluxes)   # (total_T, n_rxns)
		bounds = np.vstack(all_bounds)

		# Flatten and filter to timesteps where bound > 0 (enzyme present)
		f_flat = fluxes.ravel()
		b_flat = bounds.ravel()
		valid = (b_flat > 0) & np.isfinite(f_flat) & np.isfinite(b_flat)
		f_valid = f_flat[valid]
		b_valid = b_flat[valid]
		ratios = f_valid / b_valid

		n_violations = int(np.sum(ratios > 1.0))
		n_total = len(ratios)
		pct_active = 100 * np.mean(ratios > 0.9)  # bound is "active" (flux > 90% of bound)

		# --- Plot ------------------------------------------------------------
		fig, axes = plt.subplots(1, 2, figsize=(14, 6))

		# Panel 1: log-log scatter flux vs bound
		ax = axes[0]
		colors = np.where(ratios > 1.0, '#d62728', '#1f77b4')
		ax.scatter(b_valid, f_valid, c=colors, alpha=0.15, s=4, rasterized=True)
		lims = [
			min(b_valid[b_valid > 0].min(), f_valid[f_valid > 0].min()) * 0.5,
			max(b_valid.max(), f_valid.max()) * 2,
		]
		ax.plot(lims, lims, 'k--', lw=1, label='flux = bound (1:1)')
		ax.set_xlim(lims)
		ax.set_ylim(lims)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlabel('Computed upper bound (mmol/g DCW/h)', fontsize=11)
		ax.set_ylabel('Actual flux (mmol/g DCW/h)', fontsize=11)
		ax.set_title(
			f'Flux vs kcat upper bound\n'
			f'({n_violations}/{n_total} violations >1.0  |  '
			f'{pct_active:.1f}% timesteps with flux >90% of bound)',
			fontsize=10)
		ax.legend(fontsize=9)

		quantile = getattr(metabolism, 'kcat_estimate_quantile', 'unknown')
		multiplier = getattr(metabolism, 'kcat_estimate_multiplier', '?')
		fig.suptitle(
			f'kcat upper bound check  —  quantile: {quantile},  '
			f'multiplier: {multiplier}  —  '
			f'{len(cell_paths)} cells, gens {IGNORE_FIRST_N_GENS}+',
			fontsize=12)

		# Panel 2: histogram of flux/bound ratios
		ax2 = axes[1]
		bins = np.logspace(np.log10(max(ratios.min(), 1e-4)), np.log10(max(ratios.max(), 1.1)), 80)
		ax2.hist(ratios, bins=bins, color='#1f77b4', edgecolor='none', alpha=0.8)
		ax2.axvline(1.0, color='k', ls='--', lw=1.5, label='flux = bound')
		ax2.set_xscale('log')
		ax2.set_xlabel('flux / bound', fontsize=11)
		ax2.set_ylabel('Timestep-reaction count', fontsize=11)
		ax2.set_title(
			f'Distribution of flux/bound ratios\n'
			f'(n = {n_total:,} timestep-reaction pairs where bound > 0)',
			fontsize=10)
		ax2.legend(fontsize=9)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
