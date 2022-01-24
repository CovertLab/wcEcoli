"""
Save data to share with EcoCyc to display on the simulation tab page.

TODO:
	save a specific form of simulation metadata to share? (or existing metadata file?)
	other values
		TPM for transcripts
		weighted average for counts (time step weighted and cell cycle progress weighted)
		exponential/stochastic/constant determination
	other molecules
		values for complexes?
		charged/uncharged tRNA
		rRNA including in ribosomes/complexes
"""

import csv
import os
import pickle

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


def save_file(out_dir, filename, ids, counts, concentrations, relative_counts, relative_masses):
	"""
	TODO:
		generalize the data passed into this function to allow for arbitrary columns and labels for each column
		print header and column descriptions
	"""

	output_file = os.path.join(out_dir, filename)
	print(f'Saving data to {output_file}')
	with open(output_file, 'w') as f:
		writer = csv.writer(f, delimiter='\t')

		for id_, count, conc, rel_count, rel_mass in zip(ids, counts.T,
				concentrations.T, relative_counts.T, relative_masses.T):
			writer.writerow([id_, count.mean(), count.std(), conc.mean(), conc.std(),
				rel_count.mean(), rel_count.std(), rel_mass.mean(), rel_mass.std()])


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Load listener data
		## Tables
		mrna_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'mRNACounts'))
		monomer_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))

		## Attributes
		mrna_ids = mrna_reader.readAttribute('mRNA_ids')
		monomer_ids = monomer_reader.readAttribute('monomerIds')

		## Columns
		## remove_first=True because countsToMolar is 0 at first time step
		mrna_counts = read_stacked_columns(cell_paths, 'mRNACounts', 'mRNA_counts', remove_first=True)
		monomer_counts = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts', remove_first=True)
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar', remove_first=True)

		mrna_conc = mrna_counts * counts_to_molar
		mrna_relative_counts = mrna_counts / mrna_counts.sum(1).reshape(-1, 1)
		mrna_mw = sim_data.getter.get_masses(mrna_ids).asNumber(units.g / units.mol)
		mrna_masses = mrna_counts * mrna_mw
		mrna_relative_masses = mrna_masses / mrna_masses.sum(1).reshape(-1, 1)

		monomer_conc = monomer_counts * counts_to_molar
		monomer_relative_counts = monomer_counts / monomer_counts.sum(1).reshape(-1, 1)
		monomer_mw = sim_data.getter.get_masses(monomer_ids).asNumber(units.g / units.mol)
		monomer_masses = monomer_counts * monomer_mw
		monomer_relative_masses = monomer_masses / monomer_masses.sum(1).reshape(-1, 1)

		# Save data in tables
		mrna_ecocyc_ids = [rna[:-7] for rna in mrna_ids]  # strip _RNA[c]
		monomer_ecocyc_ids = [monomer[:-3] for monomer in monomer_ids]  # strip [*]
		media_id = 'MIX0-51'  # TODO: have a map of condition to EcoCyc media ID (temporarily hard-coded for minimal glc)
		save_file(plotOutDir, f'wcm-mrna-data-{media_id}.tsv', mrna_ecocyc_ids,
			mrna_counts, mrna_conc, mrna_relative_counts, mrna_relative_masses)
		save_file(plotOutDir, f'wcm-monomer-data-{media_id}.tsv', monomer_ecocyc_ids,
			monomer_counts, monomer_conc, monomer_relative_counts, monomer_relative_masses)


if __name__ == '__main__':
	Plot().cli()
