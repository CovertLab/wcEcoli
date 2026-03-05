"""
TODO: write this
"""

import csv
import json
import os
import pickle

import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import (read_stacked_bulk_molecules,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

IGNORE_FIRST_N_GENS = 4

# TODO: put my list(s) in here

# TODO: better way to handle this vs hardcoding reading in CSVs
below_line_directory = "reconstruction/ecoli/scripts/new_gene_below_line_proteome_ids/"
below_line_monomer_ids_filepath = below_line_directory + "below_line_monomer_ids_variant16.csv"
below_line_complex_ids_filepath = below_line_directory + "below_line_complex_ids_variant16.csv"
below_line_essential_monomer_ids_filepath = below_line_directory + "below_line_essential_monomer_ids_variant16.csv"
below_line_essential_complex_ids_filepath = below_line_directory + "below_line_essential_complex_ids_variant16.csv"

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# Get lists of monomers and complexes associated with proteome fraction decreases in high GFP sims
		# TODO: reduce redundancy
		below_line_monomer_ids = []
		with open(below_line_monomer_ids_filepath, 'r') as f:
			reader = csv.reader(f)
			for row in reader:
				below_line_monomer_ids.append(row[0])
		below_line_monomer_ids = below_line_monomer_ids[1:]
		below_line_complex_ids = []
		with open(below_line_complex_ids_filepath, 'r') as f:
			reader = csv.reader(f)
			for row in reader:
				below_line_complex_ids.append(row[0])
		below_line_complex_ids = below_line_complex_ids[1:]
		below_line_essential_monomer_ids = []
		with open(below_line_essential_monomer_ids_filepath, 'r') as f:
			reader = csv.reader(f)
			for row in reader:
				below_line_essential_monomer_ids.append(row[0])
		below_line_essential_monomer_ids = below_line_essential_monomer_ids[1:]
		below_line_essential_complex_ids = []
		with open(below_line_essential_complex_ids_filepath, 'r') as f:
			reader = csv.reader(f)
			for row in reader:
				below_line_essential_complex_ids.append(row[0])
		below_line_essential_complex_ids = below_line_essential_complex_ids[1:]
		below_line_ids = list(set(below_line_monomer_ids + below_line_complex_ids))
		below_line_essential_ids = list(set(below_line_essential_monomer_ids + below_line_essential_complex_ids))

		ap = AnalysisPaths(variantDir, cohort_plot=True)

		# Ignore first N generations
		cell_paths = ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, ap.n_generation),
			only_successful=True)

		if len(cell_paths) == 0:
			print('Skipping analysis -- not enough simulations run.')
			return

		# Read columns
		# remove_first=True because countsToMolar is 0 at first time step
		counts_to_molar = read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True, ignore_exception=True)

		# Load attributes for metabolic fluxes
		cell_density = sim_data.constants.cell_density
		listener_fba_reaction_ids = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults')).readAttribute('reactionIDs')
		listener_catalyst_ids = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults')).readAttribute('catalyst_ids')

		# Determine the catalysts associated with the below line monomers and complexes
		below_line_catalyst_ids = list(set(below_line_ids) & set(listener_catalyst_ids))
		below_line_essential_catalyst_ids = list(set(below_line_essential_ids) & set(listener_catalyst_ids))

		# Map catalyst IDs to reaction IDs
		# TODO: think about how to handle reactions with multiple catalysts
		reaction_id_to_catalyst_ids_dict = sim_data.process.metabolism.reaction_catalysts
		catalyst_id_to_reaction_ids_dict = {}
		for reaction_id, catalyst_ids in reaction_id_to_catalyst_ids_dict.items():
			for catalyst_id in catalyst_ids:
				if catalyst_id not in catalyst_id_to_reaction_ids_dict:
					catalyst_id_to_reaction_ids_dict[catalyst_id] = []
				catalyst_id_to_reaction_ids_dict[catalyst_id].append(reaction_id)
		below_line_reaction_ids = [reaction_id for catalyst_id in below_line_catalyst_ids for reaction_id in catalyst_id_to_reaction_ids_dict[catalyst_id]]
		below_line_associated_catalyst_ids = [catalyst_id for reaction_id in below_line_reaction_ids for catalyst_id in reaction_id_to_catalyst_ids_dict[reaction_id]]
		below_line_essential_reaction_ids = [reaction_id for catalyst_id in below_line_essential_catalyst_ids for reaction_id in catalyst_id_to_reaction_ids_dict[catalyst_id]]
		below_line_essential_associated_catalyst_ids = [catalyst_id for reaction_id in below_line_essential_reaction_ids for catalyst_id in reaction_id_to_catalyst_ids_dict[reaction_id]]

		below_line_pairs = zip(below_line_reaction_ids, below_line_associated_catalyst_ids)
		pair_counts = {}
		for pair in below_line_pairs:
			if pair not in pair_counts:
				pair_counts[pair] = 0
			pair_counts[pair] += 1
		unique_below_line_pairs = [pair for pair, count in pair_counts.items() if count == 1]
		below_line_reaction_ids, below_line_associated_catalyst_ids = zip(*unique_below_line_pairs) # Note the same reaction id can still appear multiple times if it is associated with multiple catalysts
		below_line_reaction_idx = np.array([np.where(np.array(listener_fba_reaction_ids) == reaction_id)[0][0] for reaction_id in below_line_reaction_ids])
		below_line_associated_catalyst_idx = np.array([np.where(np.array(listener_catalyst_ids) == catalyst_id)[0][0] for catalyst_id in below_line_associated_catalyst_ids])

		below_line_essential_pairs = zip(below_line_essential_reaction_ids, below_line_essential_associated_catalyst_ids)
		essential_pair_counts = {}
		for pair in below_line_essential_pairs:
			if pair not in essential_pair_counts:
				essential_pair_counts[pair] = 0
			essential_pair_counts[pair] += 1
		unique_below_line_essential_pairs = [pair for pair, count in essential_pair_counts.items() if count == 1]
		below_line_essential_reaction_ids, below_line_essential_associated_catalyst_ids = zip(*unique_below_line_essential_pairs)
		below_line_essential_reaction_idx = np.array([np.where(np.array(listener_fba_reaction_ids) == reaction_id)[0][0] for reaction_id in below_line_essential_reaction_ids])
		below_line_essential_associated_catalyst_idx = np.array([np.where(np.array(listener_catalyst_ids) == catalyst_id)[0][0] for catalyst_id in below_line_essential_associated_catalyst_ids])

		# Read columns
		catalyst_counts = read_stacked_columns(
			cell_paths, 'FBAResults', 'catalyst_counts', ignore_exception=True, remove_first = True)
		cell_mass = read_stacked_columns(
			cell_paths, 'Mass', 'cellMass', ignore_exception=True, remove_first = True)
		dry_mass = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass', ignore_exception=True, remove_first = True)
		conversion_coeffs = (
			dry_mass / cell_mass
			* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
		)

		# Calculate flux in units of mmol/g DCW/h
		fluxes = (
			(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
			* (read_stacked_columns(cell_paths, 'FBAResults', 'reactionFluxes', ignore_exception=True, remove_first = True) / conversion_coeffs)
			).asNumber(units.mmol / units.g / units.h)

		below_line_reaction_fluxes = fluxes[:, below_line_reaction_idx]
		below_line_essential_reaction_fluxes = fluxes[:, below_line_essential_reaction_idx]
		below_line_associated_catalyst_counts = catalyst_counts[:, below_line_associated_catalyst_idx]
		below_line_essential_associated_catalyst_counts = catalyst_counts[:, below_line_essential_associated_catalyst_idx]
		below_line_associated_catalyst_conc = below_line_associated_catalyst_counts * counts_to_molar
		below_line_essential_associated_catalyst_conc = below_line_essential_associated_catalyst_counts * counts_to_molar

		below_line_reaction_fluxes[np.isinf(below_line_reaction_fluxes)] = np.nan
		below_line_essential_reaction_fluxes[np.isinf(below_line_essential_reaction_fluxes)] = np.nan
		below_line_reaction_flux_avg = np.nanmean(below_line_reaction_fluxes, axis=0)
		below_line_essential_reaction_flux_avg = np.nanmean(below_line_essential_reaction_fluxes, axis=0)
		below_line_associated_catalyst_conc_avg = np.nanmean(below_line_associated_catalyst_conc, axis=0)
		below_line_essential_associated_catalyst_conc_avg = np.nanmean(below_line_essential_associated_catalyst_conc, axis=0)

		below_line_kcat_est = below_line_reaction_fluxes / below_line_associated_catalyst_conc
		below_line_essential_kcat_est = below_line_essential_reaction_fluxes / below_line_essential_associated_catalyst_conc

		below_line_kcat_est[np.isinf(below_line_kcat_est)] = np.nan
		below_line_essential_kcat_est[np.isinf(below_line_essential_kcat_est)] = np.nan
		below_line_kcat_est_avg = np.nanmean(below_line_kcat_est, axis=0)
		below_line_essential_kcat_est_avg = np.nanmean(below_line_essential_kcat_est, axis=0)
		below_line_kcat_est_std = np.nanstd(below_line_kcat_est, axis=0)
		below_line_essential_kcat_est_std = np.nanstd(below_line_essential_kcat_est, axis=0)

		below_line_kcat_est_avg_avg = np.nanmean(below_line_kcat_est_avg)
		below_line_kcat_est_avg_std = np.nanstd(below_line_kcat_est_avg)
		below_line_essential_kcat_est_avg_avg = np.nanmean(below_line_essential_kcat_est_avg)
		below_line_essential_kcat_est_avg_std = np.nanstd(below_line_essential_kcat_est_avg)

		# Pick 3 large ones
		vals = below_line_essential_kcat_est_avg
		finite_idx = np.flatnonzero(np.isfinite(vals))
		k = 3
		top_local = np.argpartition(vals[finite_idx], -k)[-k:]
		top_5_essential_idx = finite_idx[top_local]
		top_5_essential_reaction_ids = [below_line_essential_reaction_ids[idx] for idx in top_5_essential_idx]
		top_5_essential_catalyst_ids = [below_line_essential_associated_catalyst_ids[idx] for idx in top_5_essential_idx]
		top_5_kcat_means = below_line_essential_kcat_est_avg[top_5_essential_idx]
		top_5_kcat_stds = below_line_essential_kcat_est_std[top_5_essential_idx]

		# Pick 3 small ones
		vals = below_line_essential_kcat_est_avg
		finite_idx = np.flatnonzero(np.isfinite(vals))
		k = 3
		bottom_local = np.argpartition(vals[finite_idx], k)[:k]
		bottom_5_essential_idx = finite_idx[bottom_local]
		bottom_5_essential_reaction_ids = [below_line_essential_reaction_ids[idx] for idx in bottom_5_essential_idx]
		bottom_5_essential_catalyst_ids = [below_line_essential_associated_catalyst_ids[idx] for idx in bottom_5_essential_idx]
		bottom_5_kcat_means = below_line_essential_kcat_est_avg[bottom_5_essential_idx]
		bottom_5_kcat_stds = below_line_essential_kcat_est_std[bottom_5_essential_idx]

		# pick 3 in the median
		vals = below_line_essential_kcat_est_avg
		finite_idx = np.flatnonzero(np.isfinite(vals))
		k = 3
		median_local = np.argpartition(np.abs(vals[finite_idx] - np.nanmedian(vals)), k)[:k]
		median_5_essential_idx = finite_idx[median_local]
		median_5_essential_reaction_ids = [below_line_essential_reaction_ids[idx] for idx in median_5_essential_idx]
		median_5_essential_catalyst_ids = [below_line_essential_associated_catalyst_ids[idx] for idx in median_5_essential_idx]
		median_5_kcat_means = below_line_essential_kcat_est_avg[median_5_essential_idx]
		median_5_kcat_stds = below_line_essential_kcat_est_std[median_5_essential_idx]

		# save to csv
		output_filepath = os.path.join(plotOutDir, plotOutFileName)
		with open(output_filepath, 'w', newline='') as f:
			writer = csv.writer(f)
			writer.writerow(['reaction_id', 'catalyst_id', 'kcat_mean', 'kcat_std'])
			for reaction_id, catalyst_id, kcat_mean, kcat_std in zip(top_5_essential_reaction_ids, top_5_essential_catalyst_ids, top_5_kcat_means, top_5_kcat_stds):
				writer.writerow([reaction_id, catalyst_id, kcat_mean, kcat_std])
			for reaction_id, catalyst_id, kcat_mean, kcat_std in zip(bottom_5_essential_reaction_ids, bottom_5_essential_catalyst_ids, bottom_5_kcat_means, bottom_5_kcat_stds):
				writer.writerow([reaction_id, catalyst_id, kcat_mean, kcat_std])
			for reaction_id, catalyst_id, kcat_mean, kcat_std in zip(median_5_essential_reaction_ids, median_5_essential_catalyst_ids, median_5_kcat_means, median_5_kcat_stds):
				writer.writerow([reaction_id, catalyst_id, kcat_mean, kcat_std])



		# save the whole table
		output_filepath = os.path.join(plotOutDir, 'all_essential_' + plotOutFileName)
		with open(output_filepath, 'w', newline='') as f:
			writer = csv.writer(f)
			writer.writerow(['reaction_id', 'catalyst_id', 'kcat_mean', 'kcat_std'])
			for reaction_id, catalyst_id, kcat_mean, kcat_std in zip(below_line_essential_reaction_ids, below_line_essential_associated_catalyst_ids, below_line_essential_kcat_est_avg, below_line_essential_kcat_est_std):
				writer.writerow([reaction_id, catalyst_id, kcat_mean, kcat_std])





		import ipdb
		ipdb.set_trace()





if __name__ == '__main__':
	Plot().cli()
