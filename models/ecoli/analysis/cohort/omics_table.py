"""
Generates tables of data to share with EcoCyc for display on the "modeling" tab.

TODO:
	other values
		weighted average for counts (time step weighted and cell cycle progress weighted)
		max/min
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

IGNORE_FIRST_N_GENS = 6

# TODO (ggsun): Add this to sim_data somewhere?
# Maps media names used in model to IDs used in EcoCyc
MEDIA_NAME_TO_ID = {
	'minimal': 'MIX0-57',
	'minimal_minus_oxygen': 'MIX0-57-ANAEROBIC',
	'minimal_plus_amino_acids': 'MIX0-850',
	'minimal_acetate': 'MIX0-58',
	'minimal_succinate': 'MIX0-844',
	}


def save_file(out_dir, filename, columns, values):
	output_file = os.path.join(out_dir, filename)
	print(f'Saving data to {output_file}')
	with open(output_file, 'w') as f:
		writer = csv.writer(f, delimiter='\t')

		# Header for columns
		writer.writerow(['# Column descriptions:'])

		for col, desc in columns.items():
			writer.writerow([f'# {col}', desc])
		writer.writerow(list(columns.keys()))

		# Data rows
		for i in np.arange(len(values[0])):
			writer.writerow([v[i] for v in values])


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		media_name = sim_data.conditions[sim_data.condition]['nutrients']
		media_id = MEDIA_NAME_TO_ID.get(media_name, media_name)

		ap = AnalysisPaths(variantDir, cohort_plot=True)

		# Ignore first N generations
		cell_paths = ap.get_cells(generation=[7], seed=[0])

		# Load attributes for enzymes
		all_enzyme_ids = []
		for enzymes in sim_data.process.metabolism.reaction_catalysts.values():
			all_enzyme_ids.extend(enzymes)
		all_enzyme_ids = list(set(all_enzyme_ids))

		# Read columns
		# remove_first=True because countsToMolar is 0 at first time step
		(enzyme_counts, ) = read_stacked_bulk_molecules(
			cell_paths, (all_enzyme_ids, ), remove_first=True,
			ignore_exception=True)

		columns = {
			'id': 'Object ID, according to EcoCyc',
			't=0': 'Integer',
			't=1': 'Integer',
			't=2': 'Integer',
			't=3': 'Integer',
			't=4': 'Integer',
			't=5': 'Integer',
			't=6': 'Integer',
			}
		values = [
			[enzyme_id[:-3] for enzyme_id in all_enzyme_ids],
			enzyme_counts[0, :], enzyme_counts[500, :],
			enzyme_counts[1000, :], enzyme_counts[1500, :],
			enzyme_counts[2000, :], enzyme_counts[2500, :],
			enzyme_counts[2906, :]
			]

		save_file(
			plotOutDir, f'wcm_enzymes_{media_id}.tsv', columns, values)

		# Load attributes for metabolic fluxes
		cell_density = sim_data.constants.cell_density
		reaction_ids = sim_data.process.metabolism.base_reaction_ids

		# Read columns
		cell_mass = read_stacked_columns(
			cell_paths, 'Mass', 'cellMass', ignore_exception=True)
		dry_mass = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass', ignore_exception=True)
		conversion_coeffs = (
			dry_mass / cell_mass
			* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
		)

		# Calculate flux in units of mmol/g DCW/h
		fluxes = (
			(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
			* (read_stacked_columns(cell_paths, 'FBAResults', 'base_reaction_fluxes', ignore_exception=True) / conversion_coeffs)
			).asNumber(units.mmol / units.g / units.h)


		columns = {
			'id': 'Object ID, according to EcoCyc',
			't=0': 'Flux',
			't=1': 'Flux',
			't=2': 'Flux',
			't=3': 'Flux',
			't=4': 'Flux',
			't=5': 'Flux',
			't=6': 'Flux',
			}
		values = [
			reaction_ids, fluxes[0, :], fluxes[500, :], fluxes[1000, :],
			fluxes[1500, :], fluxes[2000, :], fluxes[2500, :], fluxes[2906, :]
			]

		save_file(
			plotOutDir, f'wcm_fluxes_{media_id}.tsv', columns, values)



if __name__ == '__main__':
	Plot().cli()
