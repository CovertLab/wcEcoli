"""
Save data to share with EcoCyc to display on the simulation tab page.

TODO:
	specific form of metadata to share
	media condition labeled filename
	header describing columns
	other values
		TPM for transcripts
		weighted average for counts
		exponential/stochastic/constant?
	other molecules
		values for complexes?
		charged/uncharged tRNA
		rRNA including in ribosomes/complexes
"""

import os
import pickle

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import read_stacked_columns
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


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

		# TODO: save mrna_ and monomer_ data (mean/std): _counts, _conc, _relative_counts, relative_masses


if __name__ == '__main__':
	Plot().cli()
