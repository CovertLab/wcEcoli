"""
Central carbon metabolism comparison to Toya et al for figure 3c
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		cell_paths = self.ap.get_cells()

		validation_data = self.read_pickle_file(validationDataFile)
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		sim_data = self.read_pickle_file(simDataFile)

		modelFluxes = {}
		toyaOrder = []

		for rxn in toyaReactions:
			modelFluxes[rxn] = []
			toyaOrder.append(rxn)

		mmol_per_g_per_h = units.mmol / units.g / units.h
		for simDir in cell_paths:
			simOutDir = os.path.join(simDir, "simOut")

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			coefficient = dryMass / cellMass * sim_data.constants.cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			base_reaction_ids = fbaResults.readAttribute("base_reaction_ids")
			base_reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fbaResults.readColumn("base_reaction_fluxes").T / coefficient).T
			base_reaction_id_to_index = {
				rxn_id: i for (i, rxn_id) in enumerate(base_reaction_ids)
			}

			for toyaReaction in toyaReactions:
				if toyaReaction in base_reaction_id_to_index:
					rxn_index = base_reaction_id_to_index[toyaReaction]
					fluxTimeCourse = base_reaction_fluxes[:, rxn_index]
					modelFluxes[toyaReaction].append(
						np.mean(fluxTimeCourse).asNumber(mmol_per_g_per_h))

		toyaVsReactionAve = []
		for rxn, toyaFlux in toyaFluxesDict.items():
			if rxn in modelFluxes:
				toyaVsReactionAve.append((np.mean(modelFluxes[rxn]), toyaFlux.asNumber(units.mmol / units.g / units.h), np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(units.mmol / units.g / units.h)))

		toyaVsReactionAve = np.array(toyaVsReactionAve)
		correlationCoefficient = np.corrcoef(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])[0,1]

		plt.figure(figsize = (8, 8))
		plt.title("Central Carbon Metabolism Flux, Pearson R = {:.2}".format(correlationCoefficient))
		plt.errorbar(toyaVsReactionAve[:,1], toyaVsReactionAve[:,0], xerr = toyaVsReactionAve[:,3], yerr = toyaVsReactionAve[:,2], fmt = "o", ecolor = "k")
		ylim = plt.ylim()
		plt.plot([ylim[0], ylim[1]], [ylim[0], ylim[1]], color = "k")
		plt.xlabel("Toya 2010 Reaction Flux [mmol/g/hr]")
		plt.ylabel("Mean WCM Reaction Flux [mmol/g/hr]")
		ax = plt.gca()
		ax.set_ylim(plt.xlim())
		whitePadSparklineAxis(ax)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
