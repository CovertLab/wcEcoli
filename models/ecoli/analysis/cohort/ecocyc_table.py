"""
Generates tables of data to share with EcoCyc for display on the "modeling" tab
of each gene page.

TODO:
	save a specific form of simulation metadata to share? (or existing metadata file?)
	save output as one file?
		stacked rows
		combine monomer values on the same row with a map of monomer to mRNA
	other values
		TPM for transcripts
		weighted average for counts (time step weighted and cell cycle progress weighted)
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

# TODO (ggsun): Add this to sim_data somewhere?
# Maps media names used in model to IDs used in EcoCyc
MEDIA_NAME_TO_ID = {
	'minimal': 'MIX0-51',
	'minimal_minus_oxygen': 'MIX0-51-ANAEROBIC',
	'minimal_plus_amino_acids': 'MIX0-847',
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

		# Ignore first two generations
		cell_paths = ap.get_cells(generation=np.arange(2, ap.n_generation))

		# Load tables and attributes for mRNAs
		mRNA_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'mRNACounts'))
		mass_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'Mass'))

		mRNA_ids = mRNA_reader.readAttribute('mRNA_cistron_ids')
		mrna_mw = sim_data.getter.get_masses(mRNA_ids).asNumber(units.fg/units.count)
		mass_unit =	mass_reader.readAttribute('cellDry_units')
		assert mass_unit == 'fg'

		# Read columns
		# remove_first=True because countsToMolar is 0 at first time step
		mRNA_counts = read_stacked_columns(
			cell_paths, 'mRNACounts', 'mRNA_cistron_counts',
			remove_first=True, ignore_exception=True)
		counts_to_molar = read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True, ignore_exception=True)
		dry_masses = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass',
			remove_first=True, ignore_exception=True)

		# Calculate derived mRNA values
		mRNA_counts_avg = mRNA_counts.mean(axis=0)
		mRNA_counts_std = mRNA_counts.std(axis=0)
		mRNA_conc = mRNA_counts * counts_to_molar
		mRNA_conc_avg = mRNA_conc.mean(axis=0)
		mRNA_conc_std = mRNA_conc.std(axis=0)
		mRNA_counts_relative_to_total_mRNA_counts = mRNA_counts_avg / mRNA_counts_avg.sum()
		mRNA_masses_avg = mRNA_counts_avg * mrna_mw
		mrna_masses_relative_to_total_mrna_mass = mRNA_masses_avg / mRNA_masses_avg.sum()
		mrna_masses_relative_to_total_dcw = mRNA_masses_avg / dry_masses.mean()

		# Save RNA data in table
		cistron_id_to_gene_id = {
			cistron['id']: cistron['gene_id']
			for cistron in sim_data.process.transcription.cistron_data
			}
		gene_ids = [cistron_id_to_gene_id[x] for x in mRNA_ids]

		columns = {
			'id': 'Object ID, according to EcoCyc',
			'rna-count-avg': 'A floating point number',
			'rna-count-std': 'A floating point number',
			'rna-concentration-avg': 'A floating point number in mM units',
			'rna-concentration-std': 'A floating point number in mM units',
			'relative-rna-count-to-total-rna-counts': 'A floating point number',
			'relative-rna-mass-to-total-rna-mass': 'A floating point number',
			'relative-rna-mass-to-total-cell-dry-mass': 'A floating point number',
			}
		values = [
			gene_ids, mRNA_counts_avg, mRNA_counts_std, mRNA_conc_avg,
			mRNA_conc_std, mRNA_counts_relative_to_total_mRNA_counts,
			mrna_masses_relative_to_total_mrna_mass,
			mrna_masses_relative_to_total_dcw,
			]

		save_file(
			plotOutDir, f'wcm-mrna-data-{media_id}.tsv', columns, values)

		# Load tables and attributes for proteins
		monomer_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
		monomer_ids = monomer_reader.readAttribute('monomerIds')

		# Read columns
		# remove_first=True because countsToMolar is 0 at first time step
		monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts',
			remove_first=True, ignore_exception=True)

		# Read validation data
		protein_id_to_schmidt_counts = {
			item[0]: item[1] for item in validation_data.protein.schmidt2015Data
			}
		protein_counts_val = np.array([
			protein_id_to_schmidt_counts.get(protein_id, np.nan) for protein_id in monomer_ids
			])

		# Calculate derived protein values
		monomer_counts_avg = monomer_counts.mean(axis=0)
		monomer_counts_std = monomer_counts.std(axis=0)
		monomer_conc = monomer_counts * counts_to_molar
		monomer_conc_avg = monomer_conc.mean(axis=0)
		monomer_conc_std = monomer_conc.std(axis=0)
		monomer_counts_relative_to_total_monomer_counts = monomer_counts_avg / monomer_counts_avg.sum()
		monomer_mass_avg = monomer_counts_avg * mrna_mw
		monomer_mass_relative_to_total_monomer_mass = monomer_mass_avg / monomer_mass_avg.sum()
		monomer_mass_relative_to_total_dcw = monomer_mass_avg / dry_masses.mean()

		# Save monomer data in table
		monomer_ecocyc_ids = [monomer[:-3] for monomer in monomer_ids]  # strip [*]

		columns = {
			'id': 'Object ID, according to EcoCyc',
			'protein-count-avg': 'A floating point number',
			'protein-count-std': 'A floating point number',
			'protein-concentration-avg': 'A floating point number in mM units',
			'protein-concentration-std': 'A floating point number in mM units',
			'relative-protein-count-to-protein-rna-counts': 'A floating point number',
			'relative-protein-mass-to-total-protein-mass': 'A floating point number',
			'relative-protein-mass-to-total-cell-dry-mass': 'A floating point number',
			'validation-count': 'A floating point number',
			}
		values = [
			monomer_ecocyc_ids, monomer_counts_avg, monomer_counts_std,
			monomer_conc_avg, monomer_conc_std,
			monomer_counts_relative_to_total_monomer_counts,
			monomer_mass_relative_to_total_monomer_mass,
			monomer_mass_relative_to_total_dcw,
			protein_counts_val,
			]

		save_file(
			plotOutDir, f'wcm-monomer-data-{media_id}.tsv', columns, values)


if __name__ == '__main__':
	Plot().cli()
