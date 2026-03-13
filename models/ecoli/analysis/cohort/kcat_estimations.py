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
import os
import pickle

import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

IGNORE_FIRST_N_GENS = 4

below_line_directory = "reconstruction/ecoli/scripts/new_gene_below_line_proteome_ids/"
below_line_monomer_ids_filepath = below_line_directory + "below_line_monomer_ids_variant16.csv"
below_line_complex_ids_filepath = below_line_directory + "below_line_complex_ids_variant16.csv"
below_line_essential_monomer_ids_filepath = below_line_directory + "below_line_essential_monomer_ids_variant16.csv"
below_line_essential_complex_ids_filepath = below_line_directory + "below_line_essential_complex_ids_variant16.csv"

QUANTILES = {
	'median': 0.50,
	'p05': 0.05,
	'p10': 0.10,
	'p90': 0.90,
	'p95': 0.95,
	'p99': 0.99,
}


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

		# Accumulate per-pair kcat values across all cells, one cell at a time
		# kcat_values[i] is a list of scalar kcat estimates for pair i
		kcat_values = [[] for _ in range(n_pairs)]

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

			# Accumulate per-pair values
			for i in range(n_pairs):
				vals = kcat[:, i]
				vals = vals[~np.isnan(vals)]
				if len(vals) > 0:
					kcat_values[i].extend(vals.tolist())

		# Write one TSV per quantile
		tsv_header = ['reaction_id', 'catalyst_id', 'kcat_estimate']
		for label, q in QUANTILES.items():
			out_path = os.path.join(plotOutDir, f'kcat_estimates_{label}.tsv')
			with open(out_path, 'w', newline='') as f:
				writer = csv.writer(f, delimiter='\t')
				writer.writerow(tsv_header)
				for i in range(n_pairs):
					vals = kcat_values[i]
					if len(vals) == 0:
						estimate = np.nan
					else:
						estimate = np.quantile(vals, q)
					writer.writerow([
						valid_pair_reaction_ids[i],
						valid_pair_catalyst_ids[i],
						estimate,
					])
			print(f"Wrote {out_path}")


if __name__ == '__main__':
	Plot().cli()
