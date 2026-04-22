"""
Heatmap of kcat upper-bound saturation rates across
new_gene_trl_eff_sweep variants.

For each kcat-constrained category (1..3 in KCAT_MULTIPLIERS), emits one
heatmap subplot. Rows are constrained reactions and columns are all variants
that belong to that category (including controls). Cells show the fraction
of timesteps where flux > SATURATION_THRESHOLD * bound, pooled over all seeds
and the last 8 generations. Cells where saturation is undefined (e.g. cat 0
variants, or enzymes never present) are drawn as light gray. Reactions where
every variant has zero flux across all categories are dropped.

Also exports one CSV per category of the saturation-rate matrix.
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
	KCAT_MULTIPLIERS,
	category_label,
	variant_to_category,
)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# Bound is considered "active" when flux exceeds this fraction of the bound.
SATURATION_THRESHOLD = 0.9


def _variant_label(index):
	"""Return a short human-readable label for a variant index."""
	cat_idx, local_idx = variant_to_category(index)
	cat = category_label(cat_idx)
	if local_idx == 0:
		return f'Control\n({cat})'
	trl_eff = TRL_EFF_VALUES[local_idx - 1]
	return f'Trl {trl_eff}\n{cat}'


def _compute_saturation(ap, variants, rxn_to_pairs, constrained_rxn_ids,
		rxn_indexes, all_cat_indexes, cat_idx_local, cell_density,
		n_total_gens, ignore_first_n_gens, flux_sum):
	"""Return (n_rxns, n_variants) saturation matrix for variants.

	`flux_sum` (n_rxns,) is updated with the total absolute flux observed
	across these variants so callers can drop always-zero reactions.
	"""
	n_rxns = len(constrained_rxn_ids)
	n_variants = len(variants)
	saturation_matrix = np.full((n_rxns, n_variants), np.nan)

	for vi_col, vi in enumerate(variants):
		cells = ap.get_cells(
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

			conversion_coeffs = (
				dry_mass / cell_mass
				* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS))
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (rxn_fluxes_raw / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			flux_sum += np.nansum(np.abs(fluxes), axis=0)

			bounds = np.zeros_like(fluxes)
			for ri, rxn_id in enumerate(constrained_rxn_ids):
				for cat_idx, kcat in rxn_to_pairs[rxn_id]:
					li = cat_idx_local[cat_idx]
					conc = cat_counts[:, li] * counts_to_molar
					bounds[:, ri] = np.maximum(
						bounds[:, ri], kcat * conc)

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

	return saturation_matrix


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variant_indexes = self.ap.get_variants()
		n_total_gens = self.ap.n_generation
		ignore_first_n_gens = max(n_total_gens - 8, 0)

		# Group variants by kcat category. We iterate constrained categories
		# (1..len(KCAT_MULTIPLIERS)-1). Cat 0 variants included in the
		# heatmap as NaN/gray if they're present (they're not here because
		# we plot per category, but the spec calls for NaN treatment).
		constrained_cats = [
			c for c, m in enumerate(KCAT_MULTIPLIERS) if m is not None]

		variants_by_cat = {c: [] for c in constrained_cats}
		for vi in variant_indexes:
			cat, _ = variant_to_category(vi)
			if cat in variants_by_cat:
				variants_by_cat[cat].append(vi)

		constrained_cats = [c for c in constrained_cats if variants_by_cat[c]]
		if not constrained_cats:
			print('Skipping: no kcat-constrained variants found.')
			return

		# ------------------------------------------------------------------ #
		# Determine constrained reactions from the first available kcat     #
		# variant                                                            #
		# ------------------------------------------------------------------ #
		ref_variant = variants_by_cat[constrained_cats[0]][0]
		kb_path = self.ap.get_variant_kb(ref_variant)
		with open(kb_path, 'rb') as f:
			sim_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		if not hasattr(metabolism, 'selected_kcat_estimates'):
			print('Skipping: sim_data does not have selected_kcat_estimates.')
			return

		ref_cells = self.ap.get_cells(
			variant=[ref_variant], generation=[0], only_successful=True)
		if len(ref_cells) == 0:
			print('Skipping: no cells found for reference variant.')
			return

		fba_reader = TableReader(
			os.path.join(ref_cells[0], 'simOut', 'FBAResults'))
		listener_rxn_ids = fba_reader.readAttribute('reactionIDs')
		listener_cat_ids = fba_reader.readAttribute('catalyst_ids')

		rxn_id_to_idx = {r: i for i, r in enumerate(listener_rxn_ids)}
		cat_id_to_idx = {c: i for i, c in enumerate(listener_cat_ids)}

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
		flux_sum = np.zeros(n_rxns)

		per_cat_matrix = {}
		for cat in constrained_cats:
			per_cat_matrix[cat] = _compute_saturation(
				self.ap, variants_by_cat[cat], rxn_to_pairs,
				constrained_rxn_ids, rxn_indexes, all_cat_indexes,
				cat_idx_local, cell_density, n_total_gens,
				ignore_first_n_gens, flux_sum)

		# Drop reactions where flux == 0 everywhere (always undefined).
		keep_rxn = flux_sum > 0
		if not np.any(keep_rxn):
			print('Skipping: all reactions have zero flux.')
			return
		kept_rxn_ids = [
			r for r, k in zip(constrained_rxn_ids, keep_rxn) if k]

		# Sort rows (globally) by pooled mean saturation rate across all
		# constrained categories, descending.
		pooled = np.concatenate(
			[per_cat_matrix[cat] for cat in constrained_cats], axis=1)
		pooled = pooled[keep_rxn]
		mean_sat = np.nanmean(pooled, axis=1)
		mean_sat_sort = np.where(np.isnan(mean_sat), -np.inf, mean_sat)
		order = np.argsort(mean_sat_sort)[::-1]
		kept_rxn_ids = [kept_rxn_ids[i] for i in order]

		for cat in constrained_cats:
			per_cat_matrix[cat] = per_cat_matrix[cat][keep_rxn][order]

		# ------------------------------------------------------------------ #
		# Export one CSV per category                                         #
		# ------------------------------------------------------------------ #
		for cat in constrained_cats:
			sat = per_cat_matrix[cat]
			variants = variants_by_cat[cat]
			csv_path = os.path.join(
				plotOutDir,
				f'kcat_saturation_rates_{category_label(cat)}.csv')
			with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
				writer = csv.writer(fh)
				writer.writerow(
					['reaction_id'] + [f'variant_{vi}' for vi in variants])
				for ri, rxn_id in enumerate(kept_rxn_ids):
					row = [rxn_id] + [
						f'{sat[ri, ci]:.4f}'
						if np.isfinite(sat[ri, ci]) else 'NaN'
						for ci in range(len(variants))]
					writer.writerow(row)
			print(f'Wrote {csv_path}')

		# ------------------------------------------------------------------ #
		# Plot: one subplot per constrained category (stacked vertically)   #
		# ------------------------------------------------------------------ #
		n_show = len(kept_rxn_ids)
		max_n_variants = max(
			len(variants_by_cat[c]) for c in constrained_cats)
		fig_height = max(
			6, len(constrained_cats) * (n_show * 0.3 + 3))
		fig_width = max(14, max_n_variants * 0.8 + 4)

		fig, axes = plt.subplots(
			len(constrained_cats), 1,
			figsize=(fig_width, fig_height), squeeze=False)

		cmap = plt.cm.YlOrRd.copy()
		cmap.set_bad('lightgray')

		for ax_idx, cat in enumerate(constrained_cats):
			ax = axes[ax_idx, 0]
			sat = per_cat_matrix[cat]
			variants = variants_by_cat[cat]

			im = ax.imshow(
				sat * 100,
				aspect='auto', origin='upper',
				cmap=cmap,
				vmin=0, vmax=100,
				interpolation='nearest')
			cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
			cbar.set_label(
				f'Saturation rate (% timesteps with flux > '
				f'{SATURATION_THRESHOLD} × bound)',
				fontsize=9)

			variant_labels = [_variant_label(vi) for vi in variants]
			ax.set_xticks(np.arange(len(variants)))
			ax.set_xticklabels(variant_labels, fontsize=7, rotation=90)
			ax.set_yticks(np.arange(n_show))
			ax.set_yticklabels(kept_rxn_ids, fontsize=6)

			ax.set_xlabel(
				f'Variant ({category_label(cat)})', fontsize=10)
			ax.set_ylabel(
				'Reaction (sorted by mean saturation rate)', fontsize=10)
			ax.set_title(
				f'Kcat bound saturation rate -- {category_label(cat)}\n'
				f'(last {min(8, n_total_gens)} gens; '
				f'exp 10^{EXPRESSION_FACTOR - 1}; gray = no data)',
				fontsize=11)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
