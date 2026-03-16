"""
Time-trace plot for kcat-constrained reaction fluxes, computed upper bounds,
and enzyme concentrations across multiple generations.

For each (reaction_id, catalyst_id) pair in REACTION_CATALYST_PAIRS (or, if
that list is empty, the first N_PAIRS_TO_PLOT pairs sorted by reaction_id from
sim_data.process.metabolism.selected_kcat_estimates), two panels are shown
per row:

  Left:  Reaction flux and kcat upper bound overlaid (mmol/g DCW/h)
  Right: Enzyme concentration [E] (mmol/L)

The x-axis is simulated time in hours, concatenated across all generations of a
single seed lineage (seed 0 by default).  Vertical dashed gray lines mark
generation boundaries.

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

# Hardcode specific pairs here to override automatic selection.
# Each entry is (reaction_id, catalyst_id).  Leave empty to use
# the first N_PAIRS_TO_PLOT pairs from selected_kcat_estimates.
REACTION_CATALYST_PAIRS = [
	# ('FADSYN-RXN', 'FADSYN-CPLX[c]'),
]

# Maximum number of pairs to plot when REACTION_CATALYST_PAIRS is empty.
N_PAIRS_TO_PLOT = 100

# Ignore first N generations (initialization transient).
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

		# Build pair list
		if REACTION_CATALYST_PAIRS:
			pairs_to_plot = list(REACTION_CATALYST_PAIRS)
		else:
			all_pairs = sorted(metabolism.selected_kcat_estimates.keys())
			pairs_to_plot = all_pairs[:N_PAIRS_TO_PLOT]

		if not pairs_to_plot:
			print('Skipping: no pairs to plot.')
			return

		# Single-lineage cell paths (seed 0), post-initialization gens only
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# --- Build index lookups from the first cell ----------------------
		fba_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

		# Resolve indexes; drop pairs not present in listeners
		valid_pairs = []      # (rxn_id, cat_id)
		rxn_indexes = []
		cat_indexes = []
		kcats = []
		for rxn_id, cat_id in pairs_to_plot:
			if rxn_id not in rxn_id_to_idx:
				print(f'Warning: {rxn_id} not in FBA reactionIDs — skipping.')
				continue
			if cat_id not in cat_id_to_idx:
				print(f'Warning: {cat_id} not in catalyst_ids — skipping.')
				continue
			# Check if valid key in selected_kcat_estimates
			if (rxn_id, cat_id) not in metabolism.selected_kcat_estimates:
				print(f'Warning: ({rxn_id}, {cat_id}) not in selected_kcat_estimates — skipping.')
				continue
			valid_pairs.append((rxn_id, cat_id))
			rxn_indexes.append(rxn_id_to_idx[rxn_id])
			cat_indexes.append(cat_id_to_idx[cat_id])
			kcats.append(metabolism.selected_kcat_estimates[(rxn_id, cat_id)])

		if not valid_pairs:
			print('No valid pairs found in listener IDs.')
			return

		n_pairs = len(valid_pairs)
		rxn_indexes = np.array(rxn_indexes)
		cat_indexes = np.array(cat_indexes)
		cell_density = sim_data.constants.cell_density

		# --- Accumulate time traces one cell at a time --------------------
		# Each list element is a 1-D array over timesteps for that cell.
		time_all = []          # hours, offset-corrected
		flux_all = []          # (T, n_pairs) mmol/g DCW/h
		bound_all = []         # (T, n_pairs) mmol/g DCW/h
		conc_all = []          # (T, n_pairs) mmol/L
		gen_boundary_times = []  # hours at which each new generation starts

		time_offset = 0.0

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				t_raw = TableReader(
					os.path.join(sim_out, 'Main')
				).readColumn('time', squeeze=True)[1:]   # drop t=0 duplicate

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
							 squeeze=False)[1:]   # (T, n_pairs)

				rxn_fluxes_raw = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('reactionFluxes',
							 indices=rxn_indexes,
							 squeeze=False)[1:]   # (T, n_pairs)

			except Exception as e:
				print(f'Ignored exception reading {sim_out}: {e!r}')
				time_offset += 0.0
				continue

			# Normalize time to hours, offset by previous cells
			t_sec = t_raw - t_raw[0]   # start each cell at 0 within cell
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

			# Enzyme concentrations -> mmol/L
			concs = cat_counts * counts_to_molar[:, np.newaxis]   # (T, n_pairs)

			# Upper bounds -> mmol/g DCW/h  (kcat [L/g/h] * [E] [mmol/L])
			bounds = np.array(kcats)[np.newaxis, :] * concs   # (T, n_pairs)

			time_all.append(t_h)
			flux_all.append(fluxes)
			bound_all.append(bounds)
			conc_all.append(concs)

		if not time_all:
			print('No cells processed successfully.')
			return

		time_cat = np.concatenate(time_all)
		flux_cat = np.vstack(flux_all)    # (total_T, n_pairs)
		bound_cat = np.vstack(bound_all)
		conc_cat = np.vstack(conc_all)

		# --- Plot ---------------------------------------------------------
		# One row per pair, two columns: flux+bound overlaid | enzyme conc
		fig, axes = plt.subplots(
			n_pairs, 2,
			figsize=(12, 2.5 * n_pairs),
			sharex=True)

		# Ensure axes is always 2-D
		if n_pairs == 1:
			axes = axes[np.newaxis, :]

		quantile = getattr(metabolism, 'kcat_estimate_quantile', 'unknown')
		multiplier = getattr(metabolism, 'kcat_estimate_multiplier', '?')
		fig.suptitle(
			f'kcat flux traces  —  quantile: {quantile},  multiplier: {multiplier}  '
			f'—  gens {IGNORE_FIRST_N_GENS}+',
			fontsize=12, y=1.002)

		axes[0, 0].set_title('Reaction flux and kcat upper bound (mmol/g DCW/h)', fontsize=10)
		axes[0, 1].set_title('Enzyme conc [E] (mmol/L)', fontsize=10)

		for ri, (rxn_id, cat_id) in enumerate(valid_pairs):
			ax_flux, ax_conc = axes[ri]

			# Left panel: flux (blue) and bound (red) overlaid
			ax_flux.plot(time_cat, flux_cat[:, ri],
						 lw=0.7, color='#1f77b4', alpha=0.85, label='flux')
			ax_flux.plot(time_cat, bound_cat[:, ri],
						 lw=0.7, color='#d62728', alpha=0.85, label='bound')

			# Right panel: enzyme concentration (green)
			ax_conc.plot(time_cat, conc_cat[:, ri],
						 lw=0.7, color='#2ca02c', alpha=0.85)

			# Generation boundary lines on both panels
			for gt in gen_boundary_times[1:]:  # skip t=0
				ax_flux.axvline(gt, color='gray', lw=0.6, ls='--', alpha=0.5)
				ax_conc.axvline(gt, color='gray', lw=0.6, ls='--', alpha=0.5)

			for ax in (ax_flux, ax_conc):
				ax.set_xlim(time_cat[0], time_cat[-1])
				ax.tick_params(labelsize=7)

			# Row label and legend on left panel
			ax_flux.set_ylabel(f'{rxn_id}\n{cat_id}', fontsize=7, labelpad=4)
			if ri == 0:
				ax_flux.legend(fontsize=7, loc='upper right')

		# x-axis label on bottom row only
		for ax in axes[-1]:
			ax.set_xlabel('Time (h)', fontsize=9)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
