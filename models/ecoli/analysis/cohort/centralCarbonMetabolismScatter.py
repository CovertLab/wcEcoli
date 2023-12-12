"""
Central carbon metabolism comparison to Toya et al
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from models.ecoli.analysis import cohortAnalysisPlot


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		allDir = self.ap.get_cells()

		validation_data = self.read_pickle_file(validationDataFile)
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		sim_data = self.read_pickle_file(simDataFile)
		cellDensity = sim_data.constants.cell_density

		modelFluxes = {}
		toyaOrder = []
		for rxn in toyaReactions:
			modelFluxes[rxn] = []
			toyaOrder.append(rxn)

		mmol_per_g_per_h = units.mmol / units.g / units.h
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			coefficient = dryMass / cellMass * cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

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
				toyaVsReactionAve.append(
					(np.mean(modelFluxes[rxn]),
					toyaFlux.asNumber(mmol_per_g_per_h),
					np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(mmol_per_g_per_h)))

		toyaVsReactionAve = np.array(toyaVsReactionAve)
		idx = np.abs(toyaVsReactionAve[:,0]) < 5 * np.abs(toyaVsReactionAve[:,1])
		rWithAll = pearsonr(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])
		rWithoutOutliers = pearsonr(toyaVsReactionAve[idx,0], toyaVsReactionAve[idx,1])

		plt.figure(figsize = (3.5, 3.5))
		ax = plt.axes()
		plt.title("Central Carbon Metabolism Flux, Pearson R = %.4f, p = %s\n(%.4f, %s without outliers)" % (rWithAll[0], rWithAll[1], rWithoutOutliers[0], rWithoutOutliers[1]), fontsize = 6)
		plt.errorbar(toyaVsReactionAve[:,1], toyaVsReactionAve[:,0], xerr = toyaVsReactionAve[:,3], yerr = toyaVsReactionAve[:,2], fmt = ".", ecolor = "k", alpha = 0.5, linewidth = 0.5)
		ylim = plt.ylim()
		plt.plot([ylim[0], ylim[1]], [ylim[0], ylim[1]], color = "k")
		plt.plot(toyaVsReactionAve[:,1], toyaVsReactionAve[:,0], "ob", markeredgewidth = 0.1, alpha = 0.9)
		plt.xlabel("Toya 2010 Reaction Flux [mmol/g/hr]")
		plt.ylabel("Mean WCM Reaction Flux [mmol/g/hr]")
		whitePadSparklineAxis(ax)

		# noinspection PyTypeChecker
		ax.set_xlim([-20, 30])
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		ax.set_yticks(list(range(int(ylim[0]), int(ylim[1]) + 1, 10)))
		ax.set_xticks(list(range(int(xlim[0]), int(xlim[1]) + 1, 10)))

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_title("")
		ax.set_xticklabels([])
		ax.set_yticklabels([])

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
