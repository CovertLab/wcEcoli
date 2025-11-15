'''
Analyze results from metabolism_kinetic_objective_weight variant by comparing
the weighted and unweighted components of the objective
'''

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns

IGNORE_FIRST_N_GENS = 0

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		n_total_gens = self.ap.n_generation
		variant_indexes = self.ap.get_variants()

		homeostatic_weights = np.zeros(len(variant_indexes))
		kinetic_weights = np.zeros(len(variant_indexes))
		secretion_penalty_coeffs = np.zeros(len(variant_indexes))
		avg_homeostatic_obj_values = np.zeros(len(variant_indexes))
		avg_kinetic_obj_values = np.zeros(len(variant_indexes))
		avg_overall_obj_values = np.zeros(len(variant_indexes))

		# Loop through all variant indexes
		for i, variant_index in enumerate(variant_indexes):
			# Load sim data for this variant
			with open(self.ap.get_variant_kb(variant_index), 'rb') as f:
				variant_sim_data = pickle.load(f)

			# Get objective weights
			homeostatic_weight = 1.0
			kinetic_weight = variant_sim_data.process.metabolism.kinetic_objective_weight
			secretion_penalty = variant_sim_data.process.metabolism.secretion_penalty_coeff

			homeostatic_weights[i] = homeostatic_weight
			kinetic_weights[i] = kinetic_weight
			secretion_penalty_coeffs[i] = secretion_penalty

			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Get the homeostatic objective value
			# Weight is 1 for homeostatic so weighted is the same as unweighted
			homeostatic_obj_values = np.sum(read_stacked_columns(
				all_cells, 'FBAResults', 'homeostaticObjectiveValues',
				fun=lambda x: np.mean(x, axis = 0)), axis = 1).squeeze()

			# Note: these have the range weight applied but not the kinetic weight
			# so really this is unweighted
			kinetic_obj_values = np.sum(read_stacked_columns(
				all_cells, 'FBAResults', 'kineticObjectiveValues',
				fun=lambda x: np.mean(x, axis = 0)), axis = 1).squeeze()

			overall_obj_values = read_stacked_columns(
				all_cells, 'FBAResults', 'objectiveValue',
				fun=lambda x: np.mean(x)).squeeze()

			# Average across seeds and gens for this variant
			avg_homeostatic_obj_values[i] = np.mean(homeostatic_obj_values)
			avg_kinetic_obj_values[i] = np.mean(kinetic_obj_values)
			avg_overall_obj_values[i] = np.mean(overall_obj_values)

		# Create bar plots
		plt.figure(figsize=(8.5, 22))
		subplots = 9

		# TODO: percentage of sims that finished

		# Overall objective value
		ax = plt.subplot(subplots, 1, 1)
		plt.bar(
			[str(w) for w in kinetic_weights],
			avg_overall_obj_values,
		)
		plt.ylabel('Overall\nObjective Value', fontsize=8)
		plt.xlabel('Kinetic Objective Weight')

		# Weighted kinetic objective value
		ax = plt.subplot(subplots, 1, 2)
		plt.bar(
			[str(w) for w in kinetic_weights],
			avg_kinetic_obj_values * kinetic_weights,
		)
		plt.ylabel('Weighted Kinetic\nObjective Value', fontsize=8)
		plt.xlabel('Kinetic Objective Weight')

		# Unweighted kinetic objective value
		ax = plt.subplot(subplots, 1, 3)
		plt.bar(
			[str(w) for w in kinetic_weights],
			avg_kinetic_obj_values,
		)
		plt.ylabel('Unweighted Kinetic\nObjective Value', fontsize=8)
		plt.xlabel('Kinetic Objective Weight')

		# Homeostatic objective Value
		ax = plt.subplot(subplots, 1, 4)
		plt.bar(
			[str(w) for w in kinetic_weights],
			avg_homeostatic_obj_values,
		)
		plt.ylabel('Unweighted/Weighted Homeostatic\nObjective Value', fontsize=8)
		plt.xlabel('Kinetic Objective Weight')

		# Weighted leftover objective value
		ax = plt.subplot(subplots, 1, 5)
		plt.bar(
			[str(w) for w in kinetic_weights],
			(
				avg_overall_obj_values - avg_homeostatic_obj_values
				- (avg_kinetic_obj_values * kinetic_weights)
			),
		)
		plt.ylabel('Weighted Leftover\nObjective Value', fontsize=8)
		plt.xlabel('Kinetic Objective Weight')

		# Weighted leftover objective value divided by secretion penalty
		ax = plt.subplot(subplots, 1, 6)
		plt.bar(
			[str(w) for w in kinetic_weights],
			(
				avg_overall_obj_values - avg_homeostatic_obj_values
				- (avg_kinetic_obj_values * kinetic_weights)
			) / secretion_penalty_coeffs,
		)
		plt.ylabel('Weighted Leftover\nObjective Value\nDivided by Secretion Penalty Coeff', fontsize=8)
		plt.xlabel('Kinetic Objective Weight')

		# TODO: homeostatic objective / kinetic objective? maybe log scale?

		plt.subplots_adjust(hspace=0.6)
		exportFigure(plt, plotOutDir, '{}_obj'.format(plotOutFileName), metadata)
		plt.close('all')

		# import ipdb
		# ipdb.set_trace()

if __name__ == "__main__":
	Plot().cli()
