"""
Analysis script to estimate kcat values for below-line essential metabolic reactions.

This script:
1. Identifies catalysts (enzymes) associated with below-line proteome genes
2. Maps these catalysts to their corresponding metabolic reactions
3. Calculates kcat estimates by dividing reaction fluxes by catalyst concentrations
4. Outputs per-quantile TSV files of kcat estimates for essential reactions

Memory-efficient design: reads only the needed reaction/catalyst columns and
processes one cell at a time, accumulating per-pair kcat values, rather than
loading all cells into memory simultaneously.
"""

import csv
import datetime
import os
import pickle

import numpy as np
from scipy.ndimage import median_filter

from reconstruction.spreadsheets import CSV_DIALECT, JsonWriter
from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

IGNORE_FIRST_N_GENS = 4

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
	os.path.dirname(os.path.abspath(__file__))))))
below_line_directory = os.path.join(
	_REPO_ROOT, "reconstruction", "ecoli", "scripts",
	"new_gene_below_line_proteome_ids")
below_line_monomer_ids_filepath = os.path.join(below_line_directory, "below_line_monomer_ids_variant16.csv")
below_line_complex_ids_filepath = os.path.join(below_line_directory, "below_line_complex_ids_variant16.csv")
below_line_essential_monomer_ids_filepath = os.path.join(below_line_directory, "below_line_essential_monomer_ids_variant16.csv")
below_line_essential_complex_ids_filepath = os.path.join(below_line_directory, "below_line_essential_complex_ids_variant16.csv")

QUANTILES = {
	'median': 0.50,
	'p05': 0.05,
	'p10': 0.10,
	'p90': 0.90,
	'p95': 0.95,
	'p99': 0.99,
	'p999': 0.999,
	'max': 1.0,
}

SMOOTH_WINDOW = 10  # timesteps (seconds) for smoothed_max estimation


def _read_csv_ids(filepath):
	"""Read the first column of a CSV, skipping the header row."""
	ids = []
	with open(filepath, 'r') as f:
		reader = csv.reader(f)
		next(reader)  # skip header
		for row in reader:
			ids.append(row[0])
	return ids


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Load below-line IDs
		below_line_monomer_ids = _read_csv_ids(below_line_monomer_ids_filepath)
		below_line_complex_ids = _read_csv_ids(below_line_complex_ids_filepath)
		below_line_essential_monomer_ids = _read_csv_ids(below_line_essential_monomer_ids_filepath)
		below_line_essential_complex_ids = _read_csv_ids(below_line_essential_complex_ids_filepath)
		below_line_essential_ids = set(below_line_essential_monomer_ids + below_line_essential_complex_ids)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, ap.n_generation),
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# Read IDs and build mappings from the first cell
		fba_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'FBAResults'))
		listener_fba_reaction_ids = fba_reader.readAttribute('reactionIDs')
		listener_catalyst_ids = fba_reader.readAttribute('catalyst_ids')

		reaction_id_to_catalyst_ids_dict = sim_data.process.metabolism.reaction_catalysts
		catalyst_id_to_reaction_ids_dict = {}
		for reaction_id, cat_ids in reaction_id_to_catalyst_ids_dict.items():
			for cid in cat_ids:
				catalyst_id_to_reaction_ids_dict.setdefault(cid, []).append(reaction_id)

		# Build ordered list of (reaction_id, catalyst_id) pairs for essential reactions
		below_line_essential_catalyst_ids_unique = list(
			set(below_line_essential_ids) & set(listener_catalyst_ids))

		pair_reaction_ids = []
		pair_catalyst_ids = []
		for catalyst_id in below_line_essential_catalyst_ids_unique:
			if catalyst_id not in catalyst_id_to_reaction_ids_dict:
				print(f"Warning: Below line essential catalyst {catalyst_id} not in reaction_catalysts mapping.")
				continue
			for reaction_id in catalyst_id_to_reaction_ids_dict[catalyst_id]:
				pair_reaction_ids.append(reaction_id)
				pair_catalyst_ids.append(catalyst_id)

		# Warn on and record reaction IDs that have missing catalyst coverage, i.e. other catalysts not in this essential set
		reaction_id_to_idx = {r: i for i, r in enumerate(listener_fba_reaction_ids)}
		catalyst_id_to_idx = {c: i for i, c in enumerate(listener_catalyst_ids)}

		reactions_to_filter_out = set()
		unique_essential_reaction_ids = set(pair_reaction_ids)
		for reaction_id in unique_essential_reaction_ids:
			all_cats = reaction_id_to_catalyst_ids_dict.get(reaction_id, [])
			essential_cat_count = sum(
				1 for c in all_cats if c in below_line_essential_ids)
			if len(all_cats) > 1 and essential_cat_count < len(all_cats):
				reactions_to_filter_out.add(reaction_id)
		print(
			"Reactions with multiple catalysts where not all are below-line essential "
			f"(filtered out): {len(reactions_to_filter_out)}")

		# Resolve indexes, filtering out problematic reactions
		rxn_indexes = []
		cat_indexes = []
		valid_pair_reaction_ids = []
		valid_pair_catalyst_ids = []
		for rxn_id, cat_id in zip(pair_reaction_ids, pair_catalyst_ids):
			if rxn_id in reactions_to_filter_out:
				continue
			if rxn_id not in reaction_id_to_idx:
				print(f"Warning: Reaction {rxn_id} not found in FBA reaction IDs.")
				continue
			if cat_id not in catalyst_id_to_idx:
				print(f"Warning: Catalyst {cat_id} not found in catalyst IDs.")
				continue
			rxn_indexes.append(reaction_id_to_idx[rxn_id])
			cat_indexes.append(catalyst_id_to_idx[cat_id])
			valid_pair_reaction_ids.append(rxn_id)
			valid_pair_catalyst_ids.append(cat_id)

		rxn_indexes = np.array(rxn_indexes)
		cat_indexes = np.array(cat_indexes)
		n_pairs = len(rxn_indexes)

		if n_pairs == 0:
			print('No valid reaction-catalyst pairs found.')
			return

		cell_density = sim_data.constants.cell_density
		quantile_values = list(QUANTILES.values())

		# Accumulate per-cell quantiles for each pair: shape (n_cells_processed, n_pairs, n_quantiles).
		# Taking the median of per-cell quantiles at the end gives an accurate approximation
		# with O(n_cells * n_pairs * n_quantiles) memory instead of O(n_cells * T * n_pairs).
		per_cell_quantiles = []  # list of (n_pairs, n_quantiles) arrays, one per cell
		per_cell_smoothed_max = []  # list of (n_pairs,) arrays, one per cell

		for cell_path in cell_paths:
			sim_out = os.path.join(cell_path, 'simOut')
			try:
				counts_to_molar = TableReader(
					os.path.join(sim_out, 'EnzymeKinetics')
				).readColumn('countsToMolar', squeeze=False)[1:]  # remove first timestep

				cell_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('cellMass', squeeze=True)[1:]

				dry_mass = TableReader(
					os.path.join(sim_out, 'Mass')
				).readColumn('dryMass', squeeze=True)[1:]

				# Read only the needed catalyst columns
				catalyst_counts = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('catalyst_counts', indices=cat_indexes, squeeze=False)[1:]

				# Read only the needed reaction flux columns
				reaction_fluxes = TableReader(
					os.path.join(sim_out, 'FBAResults')
				).readColumn('reactionFluxes', indices=rxn_indexes, squeeze=False)[1:]

			except Exception as e:
				print(f"Ignored exception reading {sim_out}: {e!r}")
				continue

			# Compute conversion coefficient: dry_mass/cell_mass * density
			conversion_coeffs = (
				dry_mass / cell_mass
				* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
			)  # shape: (T,)

			# Compute fluxes in mmol/g DCW/h; shape: (T, n_pairs)
			fluxes = (
				(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
				* (reaction_fluxes / conversion_coeffs[:, np.newaxis])
			).asNumber(units.mmol / units.g / units.h)
			fluxes[np.isinf(fluxes)] = np.nan

			# Compute catalyst concentrations; shape: (T, n_pairs)
			# counts_to_molar is (T, 1) or (T,)
			if counts_to_molar.ndim > 1:
				conc = catalyst_counts * counts_to_molar
			else:
				conc = catalyst_counts * counts_to_molar[:, np.newaxis]

			# Compute kcat estimates for this cell; shape: (T, n_pairs)
			kcat = fluxes / conc
			kcat[np.isinf(kcat)] = np.nan

			# Compute all quantiles for this cell in one pass; shape: (n_pairs, n_quantiles)
			cell_q = np.nanquantile(kcat, quantile_values, axis=0).T
			per_cell_quantiles.append(cell_q)

			# Smoothed max: median filter over non-zero kcat, then take max.
			# Median filter eliminates spikes shorter than half the window.
			cell_sm = np.full(n_pairs, np.nan)
			for pair_idx in range(n_pairs):
				col = kcat[:, pair_idx].copy()
				valid = np.isfinite(col) & (col > 0)
				if not np.any(valid):
					continue
				col[~valid] = 0.0
				smoothed = median_filter(col, size=SMOOTH_WINDOW)
				# Only consider timesteps that originally had non-zero kcat
				smoothed[~valid] = np.nan
				cell_sm[pair_idx] = np.nanmax(smoothed)
			per_cell_smoothed_max.append(cell_sm)

		if not per_cell_quantiles:
			print('No cells successfully processed.')
			return

		# Aggregate smoothed_max across cells using nanmax
		smoothed_max_stacked = np.stack(per_cell_smoothed_max, axis=0)  # (n_cells, n_pairs)
		final_smoothed_max = np.nanmax(smoothed_max_stacked, axis=0)   # (n_pairs,)

		# Stack to (n_cells, n_pairs, n_quantiles) and take the median across cells
		# to get the final (n_pairs, n_quantiles) estimates.  Exception: the 'max'
		# column uses nanmax instead of nanmedian so it reflects the true global
		# maximum observed across all cells, ensuring the bound is never violated
		# by definition in the training data.
		stacked = np.stack(per_cell_quantiles, axis=0)          # (n_cells, n_pairs, n_quantiles)
		final_quantiles = np.nanmedian(stacked, axis=0)         # (n_pairs, n_quantiles)

		# For the 'max' label, use the true global maximum across all cells
		# rather than the median of per-cell maxima, so the resulting bound
		# is never violated by definition in the training data.
		max_col_idx = list(QUANTILES.keys()).index('max')
		final_quantiles[:, max_col_idx] = np.nanmax(stacked[:, :, max_col_idx], axis=0)

		# Write one TSV per quantile using tsv_writer so values are JSON-encoded
		# and compatible with JsonReader in knowledge_base_raw.
		# Rows with NaN estimates (no valid data across any cell) are skipped.
		n_gens_used = ap.n_generation - IGNORE_FIRST_N_GENS
		timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
		provenance_comment = (
			f'# Generated by models/ecoli/analysis/cohort/kcat_estimations.py'
			f' on wildtype simulations:'
			f' {ap.n_seed} seed(s),'
			f' {n_gens_used} generation(s) used (gens {IGNORE_FIRST_N_GENS}-{ap.n_generation - 1},'
			f' first {IGNORE_FIRST_N_GENS} skipped),'
			f' {len(cell_paths)} total cells.'
			f' Generated {timestamp}.\n'
		)

		fieldnames = ['reaction_id', 'catalyst_id', 'kcat_estimate']
		for col_idx, (label, _) in enumerate(QUANTILES.items()):
			out_path = os.path.join(plotOutDir, f'kcat_estimates_{label}.tsv')
			with open(out_path, 'w', encoding='utf-8', newline='') as fh:
				fh.write(provenance_comment)
				writer = JsonWriter(fh, fieldnames, dialect=CSV_DIALECT)
				writer.writeheader()
				for i in range(n_pairs):
					estimate = final_quantiles[i, col_idx]
					if np.isnan(estimate):
						continue
					writer.writerow({
						'reaction_id': valid_pair_reaction_ids[i],
						'catalyst_id': valid_pair_catalyst_ids[i],
						'kcat_estimate': estimate,
					})
			print(f"Wrote {out_path}")

		# Write smoothed_max TSV
		sm_out_path = os.path.join(plotOutDir, 'kcat_estimates_smoothed_max.tsv')
		with open(sm_out_path, 'w', encoding='utf-8', newline='') as fh:
			fh.write(provenance_comment)
			writer = JsonWriter(fh, fieldnames, dialect=CSV_DIALECT)
			writer.writeheader()
			for i in range(n_pairs):
				estimate = final_smoothed_max[i]
				if np.isnan(estimate):
					continue
				writer.writerow({
					'reaction_id': valid_pair_reaction_ids[i],
					'catalyst_id': valid_pair_catalyst_ids[i],
					'kcat_estimate': estimate,
				})
		print(f"Wrote {sm_out_path}")


if __name__ == '__main__':
	Plot().cli()
