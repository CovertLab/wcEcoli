"""
Template for cohort analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Load data
		## Simple stacking functions for data from all cells
		names = ['ATP[c]']  # Replace with desired list of names
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(counts,) = read_stacked_bulk_molecules(cell_paths, (names,))

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		#proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein],
					#			  int)
		#proteinCountsBulk = bulkMoleculeCounts[:, proteinIndexes]

		# other items from bulk molecules:
		bulkMolecules_time = monomer_counts_reader.readColumn("time")# equal to the right number of seconds? based on the length of the other attributes below
		bulkMolecules_simulationStep = monomer_counts_reader.readColumn("simulationStep")
		bulkMolecules_monomerCounts = monomer_counts_reader.readColumn("monomerCounts")



		# extract the numbers of interest
		monomerIds = monomer_counts_reader.readAttribute("monomerIds") # this is the one that matches the indexing  I used earlier to construct the listeners!

		# these likely are the actual counts to be degraded at the timestep
		bulkMolecules_pd_CR2_counts = monomer_counts_reader.readColumn("protein_deg_CR2_counts") # this is likely still zero
		bulkMolecules_pd_ES1_counts = monomer_counts_reader.readColumn("protein_deg_ES1_counts") # of the two, this is likely the one where counts is updated to the correct number

		# total counts from protein deg:
		bulkMolecules_pd_CR2__TC = monomer_counts_reader.readColumn("protein_deg_CR2__totalCount") # this is likely still zero
		bulkMolecules_pd_CR1__TC = monomer_counts_reader.readColumn("protein_deg_CR1__totalCount") # this is likely still zero


		# peptide elongation:
		bulkMolecules_pe_ES1_counts = monomer_counts_reader.readColumn("peptide_elongate_ES1_counts")

		# peptide elgonation total counts:
		bulkMolecules_pe_ES1__TC = monomer_counts_reader.readColumn("peptide_elongate_ES1__totalCount")



		# extract the data over all generations:
		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(free_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomerIds)

		bulkMolecules_pd_CR2_counts_all_gens = read_stacked_columns(cell_paths, 'MonomerCounts', 'protein_deg_CR2_counts')


		## Or iterate on each cell if additional processing is needed
		for sim_dir in cell_paths:
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))

			# Load data
			time = main_reader.readColumn('time')

			(counts,) = read_bulk_molecule_counts(simOutDir, (names,))

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
