"""
Save the average flux for each reaction as a column in a csv file.
"""

import os
import pickle

import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.utils import units

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		n_total_gens = self.ap.n_generation
		selected_variant_indexes = self.ap.get_variants()
		n_variants = len(selected_variant_indexes)

		# Load attributes for metabolic fluxes
		cell_density = sim_data.constants.cell_density
		reaction_ids = sim_data.process.metabolism.base_reaction_ids

		reaction_ids = np.reshape(reaction_ids, (len(reaction_ids), 1))
		# Initialize a pandas dataframe to store the data
		all_avg_fluxes = np.zeros((len(reaction_ids), n_variants))

		# Loop through variant indexes
		for i, variant_index in enumerate(selected_variant_indexes):
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Read columns
			cell_mass = read_stacked_columns(
				all_cells, 'Mass', 'cellMass', ignore_exception=True)
			dry_mass = read_stacked_columns(
				all_cells, 'Mass', 'dryMass', ignore_exception=True)
			conversion_coeffs = (
					dry_mass / cell_mass
					* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
			)

			# TODO: double check that this means the fluxes are normalized by cell mass
			# Calculate flux in units of mmol/g DCW/h
			fluxes = ((COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
					* (read_stacked_columns(all_cells, 'FBAResults',
					'base_reaction_fluxes', ignore_exception=True) / conversion_coeffs)
			).asNumber(units.mmol / units.g / units.h)

			# Calculate derived flux values
			fluxes_avg = fluxes.mean(axis=0)
			all_avg_fluxes[:, i] = np.copy(fluxes_avg)

		# Save all_avg_fluxes as a csv file with row names as reaction_ids
		# and col names as selected_variant_indexes
		all_labeled_avg_fluxes = np.hstack((reaction_ids, all_avg_fluxes))
		header = np.array(['base_reaction_ids'] + selected_variant_indexes)
		all_labeled_avg_fluxes = np.vstack((header, all_labeled_avg_fluxes))

		np.savetxt(
			os.path.join(plotOutDir, 'all_labeled_avg_fluxes.csv'),
			all_labeled_avg_fluxes, delimiter=',', fmt='%s')

if __name__ == "__main__":
	Plot().cli()
