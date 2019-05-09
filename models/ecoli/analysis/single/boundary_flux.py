"""
Plot boundary fluxes
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils import filepath

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

FLUX_UNITS = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS)
CONC_UNITS = COUNTS_UNITS / VOLUME_UNITS

BURN_IN_STEPS = 20 # remove initialization artifacts

def set_ticks(ax, time):
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.tick_params(which='both', direction='out', labelsize=10)
	ax.set_xticks([time.min() / 60., time.max() / 60.])

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		mass_reader = TableReader(os.path.join(simOutDir, "Mass"))
		enzyme_kinetics_reader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

		# main reader
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time")[BURN_IN_STEPS:] - initialTime

		# mass reader
		cell_mass = mass_reader.readColumn("cellMass")[BURN_IN_STEPS:]
		cell_density = sim_data.constants.cellDensity
		nAvogadro = sim_data.constants.nAvogadro

		cell_volume = np.array([mass * units.fg / cell_density for mass in cell_mass])
		counts_to_molar = np.array([1 / (nAvogadro * volume) for volume in cell_volume])

		# enzyme kinetics reader
		allTargetFluxes = FLUX_UNITS * enzyme_kinetics_reader.readColumn("targetFluxes")[BURN_IN_STEPS:]
		allActualFluxes = FLUX_UNITS * enzyme_kinetics_reader.readColumn("actualFluxes")[BURN_IN_STEPS:]
		kineticsConstrainedReactions = np.array(enzyme_kinetics_reader.readAttribute("kineticsConstrainedReactions"))
		boundaryConstrainedReactions = np.array(enzyme_kinetics_reader.readAttribute("boundaryConstrainedReactions"))

		allTargetFluxes = allTargetFluxes.asNumber(FLUX_UNITS)
		allActualFluxes = allActualFluxes.asNumber(FLUX_UNITS)

		# boundary target fluxes
		boundaryTargetFluxes = allTargetFluxes[:, len(kineticsConstrainedReactions):]
		boundaryActualFluxes = allActualFluxes[:, len(kineticsConstrainedReactions):]

		boundary_target_flux_dict = dict(zip(boundaryConstrainedReactions, boundaryTargetFluxes.T))
		boundary_actual_flux_dict = dict(zip(boundaryConstrainedReactions, boundaryActualFluxes.T))

		# get enzyme concentrations
		reaction_catalysts = sim_data.process.metabolism.reactionCatalysts
		enzyme_ids = [reaction_catalysts.get(rxn_id) for rxn_id in boundaryConstrainedReactions]
		enzyme_ids_list = list(set([enzymes for enzyme_set in enzyme_ids for enzymes in enzyme_set]))

		# bulk molecules reader
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		molecule_ids = bulkMolecules.readAttribute("objectNames")
		enzyme_indexes = np.array([molecule_ids.index(mol_id) for mol_id in enzyme_ids_list], np.int)
		enzyme_counts = bulkMolecules.readColumn("counts")[BURN_IN_STEPS:, enzyme_indexes]

		enzyme_concs = counts_to_molar * enzyme_counts.T
		enzyme_concs_mM = [[conc.asNumber(CONC_UNITS) for conc in conc_vec]
			for conc_vec in enzyme_concs]

		enzyme_concs_mM_dict = dict(zip(enzyme_ids_list, enzyme_concs_mM))

		## Plot
		rows = len(boundaryConstrainedReactions) + 1
		cols = 2

		n_char_of_reaction_id = 25 # number of characters of reaction_id string in title

		plt.figure(figsize=(5*cols, 2*rows))
		# for reaction_id, reaction_flux in boundary_actual_flux_dict.iteritems():
		for row_idx, (reaction_id, reaction_flux) in enumerate(boundary_actual_flux_dict.iteritems()):
			target_flux = boundary_target_flux_dict[reaction_id]

			# initialize flux column
			col = 1
			plot_index = row_idx * cols + col
			ax1 = plt.subplot(rows, cols, plot_index)

			# initialize enzyme concs column
			col = 2
			plot_index = row_idx * cols + col
			ax2 = plt.subplot(rows, cols, plot_index)

			# plot flux
			ax1.plot(time / 60., target_flux, color='Red', label='target flux')
			ax1.plot(time / 60., reaction_flux, color='Blue', label='actual flux')

			# plot enzyme concs
			reaction_enzymes = reaction_catalysts.get(reaction_id)
			for enzyme_id in reaction_enzymes:
				enzyme_concs = enzyme_concs_mM_dict[enzyme_id]
				ax2.plot(time / 60., enzyme_concs, label=enzyme_id)

			# add labels
			ax1.set_xlabel("time (min)", fontsize = 12)
			ax1.set_ylabel("flux (mmol/L/s)", fontsize = 12)
			ax1.set_title("%i. %s" % (row_idx, reaction_id[:n_char_of_reaction_id]), fontsize=14, y=1.15)
			set_ticks(ax1, time)
			ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)

			ax2.set_xlabel("Time (min)", fontsize = 12)
			ax2.set_ylabel("conc (mmol/L)", fontsize = 12)
			ax2.set_title("enzyme concentrations", fontsize=14, y=1.15)
			set_ticks(ax1, time)
			ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)


		plt.tight_layout()
		plt.subplots_adjust(hspace=2.0, wspace=2.0)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
