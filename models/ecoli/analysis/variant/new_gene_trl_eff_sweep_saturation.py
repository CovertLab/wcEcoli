"""
Heatmap of kcat upper-bound saturation rates across
new_gene_trl_eff_sweep variants.

Only shows kcat-constrained variants (indices KCAT_HALF_START to end).
No-kcat variants are excluded since they have no kcat constraints to saturate.

For each kcat-constrained reaction (rows) and kcat variant (columns), shows the
fraction of timesteps where flux > 0.9 * bound, pooled over all seeds and the
last 8 generations. Only timesteps where the bound is positive (enzyme present)
are counted.

Also exports a CSV of the saturation-rate matrix.
"""

import csv
import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_trl_eff_sweep import (
	TRL_EFF_VALUES,
	EXPRESSION_FACTOR,
	N_TRL_EFF,
	KCAT_HALF_START,
)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Bound is considered "active" when flux exceeds this fraction of the bound.
SATURATION_THRESHOLD = 0.9


def _variant_label(index):
	"""Return a short human-readable label for a kcat variant index."""
	local_index = index - KCAT_HALF_START

	if local_index == 0:
		return 'Control\n(KO + kcat)'

	trl_eff = TRL_EFF_VALUES[local_index - 1]
	return f'Trl {trl_eff}\n+kcat'


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)

		# Only keep kcat variants (indices >= KCAT_HALF_START)
		kcat_variants = [vi for vi in variant_indexes
						 if vi >= KCAT_HALF_START]
		if not kcat_variants:
			print('Skipping: no kcat-constrained variants found.')
			return

		n_variants = len(kcat_variants)

		# ------------------------------------------------------------------ #
		# Determine constrained reactions from the first kcat variant          #
		# ------------------------------------------------------------------ #
		kb_path = self.ap.get_variant_kb(kcat_variants[0])
		with open(kb_path, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		if not hasattr(metabolism, 'selected_kcat_estimates'):
			print('Skipping: sim_data does not have selected_kcat_estimates.')
			return

		# Get listener IDs from a representative cell to build index maps.
		ref_cells = self.ap.get_cells(
			variant=[kcat_variants[0]], generation=[0], only_successful=True)
		if len(ref_cells) == 0:
			print('Skipping: no cells found for reference variant.')
			return

		fba_reader = TableReader(
			os.path.join(ref_cells[0], 'simOut', 'FBAResults'))
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

		all_cat_indexes = sorted(
			{ci for pairs in rxn_to_pairs.values() for ci, _ in pairs})
		cat_idx_local = {ci: li for li, ci in enumerate(all_cat_indexes)}

		cell_density = sim_data.constants.cell_density
		n_rxns = len(constrained_rxn_ids)

		# ------------------------------------------------------------------ #
		# Per-variant: compute saturation rates                                #
		# ------------------------------------------------------------------ #
		saturation_matrix = np.full((n_rxns, n_variants), np.nan)

		for vi_col, vi in enumerate(kcat_variants):
			cells = self.ap.get_cells(
				variant=[vi],
				generation=np.arange(ignore_first_n_gens, n_total_gens),
				only_successful=True)
			if len(cells) == 0:
				continue

			saturated_count = np.zeros(n_rxns, dtype=int)
			valid_count = np.zeros(n_rxns, dtype=int)

			for cell_path in cells:
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
								 squeeze=False)[1:]

					rxn_fluxes_raw = TableReader(
						os.path.join(sim_out, 'FBAResults')
					).readColumn('reactionFluxes',
								 indices=rxn_indexes,
								 squeeze=False)[1:]
				except Exception as e:
					print(f'Ignored exception reading {sim_out}: {e!r}')
					continue

				# Reaction fluxes -> mmol/g DCW/h
				conversion_coeffs = (
					dry_mass / cell_mass
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
						bounds[:, ri] = np.maximum(
							bounds[:, ri], kcat * conc)

				# Count saturated and valid timesteps per reaction.
				for ri in range(n_rxns):
					valid = (bounds[:, ri] > 0) & np.isfinite(fluxes[:, ri])
					n_valid = valid.sum()
					if n_valid == 0:
						continue
					valid_count[ri] += n_valid
					saturated_count[ri] += (
						fluxes[valid, ri]
						> SATURATION_THRESHOLD * bounds[valid, ri]).sum()

			has_data = valid_count > 0
			saturation_matrix[has_data, vi_col] = (
				saturated_count[has_data] / valid_count[has_data])

		# ------------------------------------------------------------------ #
		# Sort rows by mean saturation rate (descending)                       #
		# ------------------------------------------------------------------ #
		mean_sat = np.nanmean(saturation_matrix, axis=1)
		mean_sat[np.isnan(mean_sat)] = -1
		sort_order = np.argsort(mean_sat)[::-1]

		saturation_matrix = saturation_matrix[sort_order]
		sorted_rxn_ids = [constrained_rxn_ids[i] for i in sort_order]
		mean_sat = mean_sat[sort_order]

		# Drop reactions with no data at all.
		has_any = mean_sat >= 0
		saturation_matrix = saturation_matrix[has_any]
		sorted_rxn_ids = [r for r, h in zip(sorted_rxn_ids, has_any) if h]

		if len(sorted_rxn_ids) == 0:
			print('Skipping: no reactions had valid data across variants.')
			return

		n_show = len(sorted_rxn_ids)

		# ------------------------------------------------------------------ #
		# Export CSV                                                           #
		# ------------------------------------------------------------------ #
		csv_path = os.path.join(plotOutDir, 'kcat_saturation_rates.csv')
		with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
			writer = csv.writer(fh)
			header = ['reaction_id'] + [
				f'variant_{vi}' for vi in kcat_variants]
			writer.writerow(header)
			for ri, rxn_id in enumerate(sorted_rxn_ids):
				row = [rxn_id] + [
					f'{saturation_matrix[ri, ci]:.4f}'
					if np.isfinite(saturation_matrix[ri, ci]) else 'NaN'
					for ci in range(n_variants)]
				writer.writerow(row)
		print(f'Wrote {csv_path}')

		# ------------------------------------------------------------------ #
		# Plot heatmap                                                         #
		# ------------------------------------------------------------------ #
		fig_height = max(6, n_show * 0.3 + 2)
		fig_width = max(12, n_variants * 0.7 + 3)
		fig, ax = plt.subplots(figsize=(fig_width, fig_height))

		cmap = plt.cm.YlOrRd.copy()
		cmap.set_bad(color='#cccccc')

		im = ax.imshow(
			saturation_matrix * 100,
			aspect='auto',
			origin='upper',
			cmap=cmap,
			vmin=0, vmax=100,
			interpolation='nearest')

		cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
		cbar.set_label(
			f'Saturation rate (% timesteps with flux > '
			f'{SATURATION_THRESHOLD} \u00d7 bound)',
			fontsize=9)

		variant_labels = [_variant_label(vi) for vi in kcat_variants]
		ax.set_xticks(np.arange(n_variants))
		ax.set_xticklabels(variant_labels, fontsize=7, rotation=45, ha='right')
		ax.set_yticks(np.arange(n_show))
		ax.set_yticklabels(sorted_rxn_ids, fontsize=6)

		ax.set_xlabel('Variant (kcat half only)', fontsize=10)
		ax.set_ylabel('Reaction (sorted by mean saturation rate)', fontsize=10)
		ax.set_title(
			'Kcat bound saturation rate by reaction and variant\n'
			f'(last {min(8, n_total_gens)} gens; '
			f'expression factor {EXPRESSION_FACTOR}; gray = no data)',
			fontsize=11)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
