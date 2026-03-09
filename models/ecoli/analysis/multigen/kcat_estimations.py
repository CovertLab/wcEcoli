"""
Plot estimated kcat values for below-the-line reactions over time, where kcat is
estimated as flux / catalyst concentration. This is a rough estimation that assumes
all of the catalyst is active and that the reaction is substrate-saturated, but
it can provide some insight into how kcat values may be changing over the course
of the simulation. The plot will show the estimated kcat values for each reaction
of interest over time, along with the average and standard deviation of the estimated
kcat values across all time points.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
# noinspection PyUnresolvedReferences
import numpy as np
from numpy import inf

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns, read_bulk_molecule_counts)
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

LINE_COLOR = (66/255, 170/255, 154/255)
IGNORE_FIRST_N_GENS = 4

reaction_catalyst_id_tuples = [
	("PREPHENATEDEHYDROG-RXN","CHORISMUTPREPHENDEHYDROG-CPLX[c]"),
	("MALONYL-COA-ACP-TRANSACYL-RXN","MALONYL-COA-ACP-TRANSACYL-MONOMER[c]"),
	("GLYOHMETRANS-RXN-SER/THF//GLY/METHYLENE-THF/WATER.33. (reverse)","GLYOHMETRANS-CPLX[c]"),
	("RXN-17009","1-ACYLGLYCEROL-3-P-ACYLTRANSFER-MONOMER[i]"),
	("PROTOHEMEFERROCHELAT-RXN[CCO-CYTOSOL]-PROTOHEME/PROTON//PROTOPORPHYRIN_IX/FE+2.54.","CPLX0-7810[c]"),
	("PROTOHEMEFERROCHELAT-RXN[CCO-CYTOSOL]-PROTOHEME/PROTON//PROTOPORPHYRIN_IX/FE+2.54.","G7266-MONOMER[c]"),
	("PROTOHEMEFERROCHELAT-RXN[CCO-CYTOSOL]-PROTOHEME/PROTON//PROTOPORPHYRIN_IX/FE+2.54.","PROTOHEME-FERROCHELAT-MONOMER[c]"),
	("RXN-9558 (reverse)","CPLX0-8006[c]"),
	("PSERTRANSAMPYR-RXN","PSERTRANSAM-CPLX[c]"),
	("RXN-22914","CPLX0-341[c]"),
	("UDPNACETYLGLUCOSAMENOLPYRTRANS-RXN","UDPNACETYLGLUCOSAMENOLPYRTRANS-MONOMER[c]"),
]

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# Ignore first N generations
		cell_paths = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS, self.ap.n_generation),
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

		# get indices of reactions and catalysts of interest
		below_line_reaction_idx = np.array([listener_fba_reaction_ids.index(rxn_id) for (rxn_id, catalyst_id) in reaction_catalyst_id_tuples])
		below_line_associated_catalyst_idx = np.array([listener_catalyst_ids.index(catalyst_id) for (rxn_id, catalyst_id) in reaction_catalyst_id_tuples])

		# Read columns
		below_line_associated_catalyst_counts = read_stacked_columns(
			cell_paths, 'FBAResults', 'catalyst_counts', ignore_exception=True, remove_first = True)[:, below_line_associated_catalyst_idx]
		cell_mass = read_stacked_columns(
			cell_paths, 'Mass', 'cellMass', ignore_exception=True, remove_first = True)
		dry_mass = read_stacked_columns(
			cell_paths, 'Mass', 'dryMass', ignore_exception=True, remove_first = True)
		conversion_coeffs = (
			dry_mass / cell_mass
			* cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)
		)

		# Calculate flux in units of mmol/g DCW/h
		below_line_reaction_fluxes = (
			(COUNTS_UNITS / MASS_UNITS / TIME_UNITS)
			* (read_stacked_columns(cell_paths, 'FBAResults', 'reactionFluxes', ignore_exception=True, remove_first = True) / conversion_coeffs)
			).asNumber(units.mmol / units.g / units.h)[:, below_line_reaction_idx]

		below_line_associated_catalyst_conc = below_line_associated_catalyst_counts * counts_to_molar
		below_line_kcat_est = below_line_reaction_fluxes / below_line_associated_catalyst_conc

		# create a folder called kcat_estimation_timeseries in the plotOutDir
		kcat_estimation_timeseries_dir = os.path.join(plotOutDir, 'kcat_estimation_timeseries')
		if not os.path.exists(kcat_estimation_timeseries_dir):
			os.makedirs(kcat_estimation_timeseries_dir)

		# Loop through each tuple and plot the fluxes, counts, and kcat estimations over time
		for i, (rxn_id, catalyst_id) in enumerate(reaction_catalyst_id_tuples):
			plt.figure(figsize=(10, 20))
			plt.subplot(3, 1, 1)
			plt.plot(below_line_reaction_fluxes[:, i], color=LINE_COLOR)
			plt.title(f'Flux of {rxn_id} over time')
			plt.ylabel('Flux (mmol/g DCW/h)')

			plt.subplot(3, 1, 2)
			plt.plot(below_line_associated_catalyst_conc[:, i], color=LINE_COLOR)
			plt.title(f'Concentration of {catalyst_id} over time')
			plt.ylabel('Concentration (mM)')

			plt.subplot(3, 1, 3)
			plt.plot(below_line_kcat_est[:, i], color=LINE_COLOR)
			avg_kcat_est = np.mean(below_line_kcat_est[:, i])
			std_kcat_est = np.std(below_line_kcat_est[:, i])
			plt.title(f'Estimated kcat over time (average: {avg_kcat_est:.2f}, std: {std_kcat_est:.2f})')
			plt.axhline(avg_kcat_est, color='red', linestyle='--', label='Average kcat est')
			plt.legend()
			plt.ylabel('Estimated kcat')
			plt.xlabel('Time step')

			# Save the figure
			plot_file_name = f'{rxn_id}_{catalyst_id}_kcat_estimation_timeseries'
			plot_file_name = plot_file_name.replace(' ', '_')
			plot_file_name = plot_file_name.replace('/', '|')

			print(f'Exporting figure to {os.path.join(kcat_estimation_timeseries_dir, plot_file_name)}')

			exportFigure(plt, kcat_estimation_timeseries_dir, plot_file_name, metadata)


if __name__ == '__main__':
	Plot().cli()
