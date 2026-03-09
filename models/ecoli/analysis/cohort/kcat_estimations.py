"""
Analysis script to estimate kcat values for below-line essential metabolic reactions.

This script:
1. Identifies catalysts (enzymes) associated with below-line proteome genes
2. Maps these catalysts to their corresponding metabolic reactions
3. Calculates kcat estimates by dividing reaction fluxes by catalyst concentrations
4. Outputs statistics on kcat estimates for essential reactions
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
		cell_density = sim_data.constants.cell_density
		counts_to_molar = read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True, ignore_exception=True)
		catalyst_counts = read_stacked_columns(
			cell_paths, 'FBAResults', 'catalyst_counts', ignore_exception=True, remove_first = True)
		catalyst_concentrations = catalyst_counts * counts_to_molar
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

		# Load attributes for metabolic fluxes
		listener_fba_reaction_ids = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults')).readAttribute('reactionIDs')
		listener_catalyst_ids = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'FBAResults')).readAttribute('catalyst_ids')

		# Determine the catalysts associated with the below line monomers and complexes
		below_line_catalyst_ids = list(set(below_line_ids) & set(listener_catalyst_ids))
		below_line_essential_catalyst_ids_unique = list(set(below_line_essential_ids) & set(listener_catalyst_ids))

		# Map catalyst IDs to reaction IDs
		# Note some reactions will have multiple catalysts
		reaction_id_to_catalyst_ids_dict = sim_data.process.metabolism.reaction_catalysts
		catalyst_id_to_reaction_ids_dict = {}
		for reaction_id, catalyst_ids in reaction_id_to_catalyst_ids_dict.items():
			for catalyst_id in catalyst_ids:
				if catalyst_id not in catalyst_id_to_reaction_ids_dict:
					catalyst_id_to_reaction_ids_dict[catalyst_id] = []
				catalyst_id_to_reaction_ids_dict[catalyst_id].append(reaction_id)

		# TODO: decide how to handle cases where a catalyst is associated with multiple reactions
		# maybe compute kcat estimates for all associated reactions and take the max? but save them all for now

		# Get reaction IDs associated with below line essential catalysts
		below_line_essential_catalyst_ids = [] # Will contain repeats if the catalyst is associated with multiple reactions, which is fine because we will compute kcat estimates for each reaction-catalyst pair
		below_line_essential_reaction_ids = []
		for catalyst_id in below_line_essential_catalyst_ids_unique:
			if catalyst_id in catalyst_id_to_reaction_ids_dict:
				below_line_essential_reaction_ids.extend(
					catalyst_id_to_reaction_ids_dict[catalyst_id])
				below_line_essential_catalyst_ids.extend(
					[catalyst_id] * len(catalyst_id_to_reaction_ids_dict[catalyst_id]))
			else:
				print(f"Warning: Below line essential catalyst ID {catalyst_id} not found in reaction_catalysts mapping.")

		# Determine the indexes for below line essential reactions in the flux array
		below_line_essential_reaction_indexes = []
		reaction_id_to_idx_dict = {rxn_id: idx for idx, rxn_id in enumerate(listener_fba_reaction_ids)}
		for reaction_id in below_line_essential_reaction_ids:
			if reaction_id in reaction_id_to_idx_dict:
				below_line_essential_reaction_indexes.append(reaction_id_to_idx_dict[reaction_id])
			else:
				print(f"Warning: Reaction ID {reaction_id} associated with below line essential catalysts not found in listener FBA reaction IDs.")
		below_line_essential_reaction_indexes = np.array(below_line_essential_reaction_indexes)

		# Determine the indexes for below line essential catalysts in the catalyst concentration array
		below_line_essential_catalyst_indexes = []
		catalyst_id_to_idx_dict = {catalyst_id: idx for idx, catalyst_id in enumerate(listener_catalyst_ids)}
		for catalyst_id in below_line_essential_catalyst_ids:
			if catalyst_id in catalyst_id_to_idx_dict:
				below_line_essential_catalyst_indexes.append(catalyst_id_to_idx_dict[catalyst_id])
			else:
				print(f"Warning: Catalyst ID {catalyst_id} associated with below line essential reactions not found in listener catalyst IDs.")
		below_line_essential_catalyst_indexes = np.array(below_line_essential_catalyst_indexes)

		# Calculate kcat estimates for below line essential reactions by dividing the fluxes by the associated catalyst concentrations
		below_line_essential_fluxes = fluxes[:, below_line_essential_reaction_indexes]
		below_line_essential_catalyst_concentrations = catalyst_concentrations[:, below_line_essential_catalyst_indexes]
		below_line_essential_kcat_estimates = below_line_essential_fluxes / below_line_essential_catalyst_concentrations

		# Save reaction id, catalyst id, average flux, averacge catalyst concentration, and average kcat estimate for below line essential reactions to a CSV
		output_file_avg = os.path.join(plotOutDir, 'below_line_essential_kcat_estimates_averages.csv')
		with open(output_file_avg, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['Reaction ID', 'Catalyst ID', 'Average Flux (mmol/g DCW/h)', 'Average Catalyst Concentration (mM)', 'Average kcat Estimate (1/h)'])
			for i in range(len(below_line_essential_reaction_indexes)):
				reaction_id = listener_fba_reaction_ids[below_line_essential_reaction_indexes[i]]
				catalyst_id = listener_catalyst_ids[below_line_essential_catalyst_indexes[i]]
				avg_flux = np.mean(below_line_essential_fluxes[:, i])
				avg_catalyst_conc = np.mean(below_line_essential_catalyst_concentrations[:, i])
				avg_kcat_estimate = np.mean(below_line_essential_kcat_estimates[:, i])
				writer.writerow([reaction_id, catalyst_id, avg_flux, avg_catalyst_conc, avg_kcat_estimate])

		# Save reaction id, catalyst id, median kcat estimate, and each 5% quantile kcat estimate for below line essential reactions to a CSV
		output_file_quantiles = os.path.join(plotOutDir, 'below_line_essential_kcat_estimates_quantiles.csv')
		with open(output_file_quantiles, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['Reaction ID', 'Catalyst ID', 'Median kcat Estimate (1/h)', '5% Quantile kcat Estimate (1/h)', '95% Quantile kcat Estimate (1/h)'])
			for i in range(len(below_line_essential_reaction_indexes)):
				reaction_id = listener_fba_reaction_ids[below_line_essential_reaction_indexes[i]]
				catalyst_id = listener_catalyst_ids[below_line_essential_catalyst_indexes[i]]
				median_kcat_estimate = np.median(below_line_essential_kcat_estimates[:, i])
				quantile_5_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.05)
				quantile_95_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.95)
				writer.writerow([reaction_id, catalyst_id, median_kcat_estimate, quantile_5_kcat_estimate, quantile_95_kcat_estimate])

		# Give a csv of the averages for k of the top kcat estimates, lowest kcat estimates, and closest to the median kcat estimates among the below line essential reactions
		# These are candidates for checking multigen plots
		k = 10
		output_file_top_k = os.path.join(plotOutDir, 'below_line_essential_kcat_estimates_top_k.csv')
		with open(output_file_top_k, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['Reaction ID', 'Catalyst ID', 'Average kcat Estimate (1/h)'])
			# Get indexes of top k highest average kcat estimates
			top_k_indexes = np.argsort(np.mean(below_line_essential_kcat_estimates, axis=0))[-k:]
			for i in top_k_indexes:
				reaction_id = listener_fba_reaction_ids[below_line_essential_reaction_indexes[i]]
				catalyst_id = listener_catalyst_ids[below_line_essential_catalyst_indexes[i]]
				avg_kcat_estimate = np.mean(below_line_essential_kcat_estimates[:, i])
				writer.writerow([reaction_id, catalyst_id, avg_kcat_estimate])
			# Get indexes of top k closest to the median kcat estimates
			median_kcat_estimates = np.median(below_line_essential_kcat_estimates, axis=0)
			closest_to_median_k_indexes = np.argsort(np.abs(np.mean(below_line_essential_kcat_estimates, axis=0) - median_kcat_estimates))[:k]
			for i in closest_to_median_k_indexes:
				reaction_id = listener_fba_reaction_ids[below_line_essential_reaction_indexes[i]]
				catalyst_id = listener_catalyst_ids[below_line_essential_catalyst_indexes[i]]
				avg_kcat_estimate = np.mean(below_line_essential_kcat_estimates[:, i])
				writer.writerow([reaction_id, catalyst_id, avg_kcat_estimate])
			# Get indexes of top k lowest average kcat estimates
			lowest_k_indexes = np.argsort(np.mean(below_line_essential_kcat_estimates, axis=0))[:k]
			for i in lowest_k_indexes:
				reaction_id = listener_fba_reaction_ids[below_line_essential_reaction_indexes[i]]
				catalyst_id = listener_catalyst_ids[below_line_essential_catalyst_indexes[i]]
				avg_kcat_estimate = np.mean(below_line_essential_kcat_estimates[:, i])
				writer.writerow([reaction_id, catalyst_id, avg_kcat_estimate])

		# Output a csv of the essential reaction ids where they have multiple catalysts, and how many of those catalysts are in the essential catalyst ids list
		unique_below_line_essential_reaction_ids = set(below_line_essential_reaction_ids)
		num_catalysts = []
		num_below_line_essential_catalysts = []
		reactions_to_filter_out = [] # reactions with multiple catalysts where not all of the catalysts are in the below line essential catalyst list, need to think about how to handle these
		for reaction_id in unique_below_line_essential_reaction_ids:
			reaction_catalyst_ids = reaction_id_to_catalyst_ids_dict.get(reaction_id, [])
			reaction_num_catalysts = len(reaction_catalyst_ids)
			num_catalysts.append(reaction_num_catalysts)
			reaction_num_below_line_essential_catalysts = sum([1 for this_catalyst_id in reaction_catalyst_ids if this_catalyst_id in below_line_essential_catalyst_ids_unique])
			num_below_line_essential_catalysts.append(reaction_num_below_line_essential_catalysts)
			if len(reaction_catalyst_ids) > 1 and reaction_num_below_line_essential_catalysts < reaction_num_catalysts:
				reactions_to_filter_out.append(reaction_id)
		print("Number of below line essential reactions with multiple catalysts where not all catalysts are in the below line essential catalyst list (need to filter these out for some analyses):", len(reactions_to_filter_out))
		output_file_multiple_catalysts = os.path.join(plotOutDir, 'below_line_essential_reactions_multiple_catalysts.csv')
		with open(output_file_multiple_catalysts, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['Reaction ID', 'Number of Catalysts', 'Number of Below Line Essential Catalysts'])
			for i, reaction_id in enumerate(unique_below_line_essential_reaction_ids):
				writer.writerow([reaction_id, num_catalysts[i], num_below_line_essential_catalysts[i]])

		# Save reaction id, catalyst id, median kcat estimate, and each 5% quantile kcat estimate for below line essential reactions to a CSV but filter out the reactions with multiple catalysts where not all of the catalysts are in the below line essential catalyst list
		output_file_quantiles_filtered = os.path.join(plotOutDir, 'below_line_essential_kcat_estimates_quantiles_filtered.csv')
		with open(output_file_quantiles_filtered, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerow(['Reaction ID', 'Catalyst ID', 'Median kcat Estimate (1/h)', '5% Quantile kcat Estimate (1/h)', '95% Quantile kcat Estimate (1/h)', '10% Quantile kcat Estimate (1/h)', '90% Quantile kcat Estimate (1/h)', '99% Quantile kcat Estimate (1/h)'])
			for i in range(len(below_line_essential_reaction_indexes)):
				reaction_id = listener_fba_reaction_ids[below_line_essential_reaction_indexes[i]]
				if reaction_id in reactions_to_filter_out:
					continue
				catalyst_id = listener_catalyst_ids[below_line_essential_catalyst_indexes[i]]
				median_kcat_estimate = np.median(below_line_essential_kcat_estimates[:, i])
				quantile_5_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.05)
				quantile_95_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.95)
				quantile_10_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.10)
				quantile_90_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.90)
				quantile_99_kcat_estimate = np.quantile(below_line_essential_kcat_estimates[:, i], 0.99)
				writer.writerow([reaction_id, catalyst_id, median_kcat_estimate, quantile_5_kcat_estimate, quantile_95_kcat_estimate, quantile_10_kcat_estimate, quantile_90_kcat_estimate, quantile_99_kcat_estimate])


if __name__ == '__main__':
	Plot().cli()
