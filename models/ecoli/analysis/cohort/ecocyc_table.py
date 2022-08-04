"""
Save data to share with EcoCyc to display on the simulation tab page.

TODO:
	save a specific form of simulation metadata to share? (or existing metadata file?)
	save output as one file?
		stacked rows
		combine monomer values on the same row with a map of monomer to mRNA
	other values
		TPM for transcripts
		weighted average for counts (time step weighted and cell cycle progress weighted)
		exponential/stochastic/constant determination
		max/min
	other molecules
		values for complexes?
		charged/uncharged tRNA
		rRNA including in ribosomes/complexes
"""

import csv
import numpy as np
import os
import pickle

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


def save_file(out_dir, filename, ids,
	counts, concentrations, relative_counts, relative_masses_to_total_type_mass,
	relative_masses_to_total_cell_dry_mass, validation_counts,
	aerobic=True):
	"""
	TODO:
		generalize the data passed into this function to allow for arbitrary
		columns and labels for each column
	"""

	output_file = os.path.join(out_dir, filename)
	print(f'Saving data to {output_file}')
	with open(output_file, 'w') as f:
		writer = csv.writer(f, delimiter='\t')

		# Header for columns
		writer.writerow([f'# {"aerobic" if aerobic else "anaerobic"} condition'])
		writer.writerow(['# Column descriptions:'])
		columns = {
			'id': 'Object ID, according to EcoCyc',
			# 'change-symbol': 'One of: constant, exponential, or stochastic',  TODO
			'avg-count': 'A floating point number',
			'count, standard deviation': 'A floating point number',
			'avg-concentration': 'A floating point number in mM units',
			'concentration, standard deviation': 'A floating point number in mM units',
			'avg-relative-count-to-total-counts': 'A floating point number',
			'avg-relative-mass-to-total-molecule-type-mass': 'A floating point number',
			'avg-relative-mass-to-total-cell-dry-mass': 'A floating point number',
			'validation-count': 'A floating point number',
			}
		for col, desc in columns.items():
			writer.writerow([f'# {col}', desc])
		writer.writerow(list(columns.keys()))

		# Data rows
		for id_, count, conc, rel_count, rel_mass1, rel_mass2, val_count in zip(
				ids, counts.T, concentrations.T, relative_counts.T,
				relative_masses_to_total_type_mass.T, relative_masses_to_total_cell_dry_mass.T,
				validation_counts):
			writer.writerow([id_, count.mean(), count.std(), conc.mean(), conc.std(),
				rel_count.mean(), rel_mass1.mean(), rel_mass2.mean(), val_count])


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Load listener data
		## Tables
		mrna_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'mRNACounts'))
		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		mass_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'Mass'))

		## Attributes
		mrna_ids = mrna_reader.readAttribute('mRNA_ids')
		monomer_ids = monomer_reader.readAttribute('monomerIds')
		mass_unit =	mass_reader.readAttribute('cellDry_units')
		assert mass_unit == 'fg'

		## Columns
		## remove_first=True because countsToMolar is 0 at first time step
		mrna_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_counts', remove_first=True, ignore_exception=True)
		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts', remove_first=True, ignore_exception=True)
		counts_to_molar = read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar', remove_first=True, ignore_exception=True)
		dry_masses = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass', remove_first=True, ignore_exception=True)

		# Validation
		mrna_validation_counts = np.zeros(mrna_counts.shape[1])  # Zero for now
		protein_id_to_schmidt_counts = {
			item[0]: item[1] for item in validation_data.protein.schmidt2015Data
			}
		protein_validation_counts = np.array([
			protein_id_to_schmidt_counts.get(protein_id, np.nan) for protein_id in monomer_ids
			])

		# Derived mRNA values
		mrna_conc = mrna_counts * counts_to_molar
		mrna_counts_relative_to_total_mrna_counts = mrna_counts / mrna_counts.sum(1).reshape(-1, 1)
		mrna_mw = sim_data.getter.get_masses(mrna_ids).asNumber(units.fg / units.count)
		mrna_masses = mrna_counts * mrna_mw
		mrna_masses_relative_to_total_mrna_mass = mrna_masses / mrna_masses.sum(1).reshape(-1, 1)
		mrna_masses_relative_to_total_dcw = mrna_masses / dry_masses

		# Derived monomer values
		monomer_conc = monomer_counts * counts_to_molar
		monomer_counts_relative_to_total_monomer_counts = monomer_counts / monomer_counts.sum(1).reshape(-1, 1)
		monomer_mw = sim_data.getter.get_masses(monomer_ids).asNumber(units.fg / units.count)
		monomer_masses = monomer_counts * monomer_mw
		monomer_masses_relative_to_total_monomer_mass = monomer_masses / monomer_masses.sum(1).reshape(-1, 1)
		monomer_masses_relative_to_total_dcw = monomer_masses / dry_masses

		# Save data in tables
		mrna_ecocyc_ids = [rna[:-7] for rna in mrna_ids]  # strip _RNA[c]
		monomer_ecocyc_ids = [monomer[:-3] for monomer in monomer_ids]  # strip [*]
		media_id = 'MIX0-51'  # TODO: have a map of condition to EcoCyc media ID (temporarily hard-coded for minimal glc)

		save_file(plotOutDir, f'wcm-mrna-data-{media_id}.tsv',
			mrna_ecocyc_ids, mrna_counts, mrna_conc,
			mrna_counts_relative_to_total_mrna_counts,
			mrna_masses_relative_to_total_mrna_mass,
			mrna_masses_relative_to_total_dcw,
			mrna_validation_counts,
			)
		save_file(plotOutDir, f'wcm-monomer-data-{media_id}.tsv',
			monomer_ecocyc_ids, monomer_counts, monomer_conc,
			monomer_counts_relative_to_total_monomer_counts,
			monomer_masses_relative_to_total_monomer_mass,
			monomer_masses_relative_to_total_dcw,
			protein_validation_counts,
			)


if __name__ == '__main__':
	Plot().cli()
