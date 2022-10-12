"""
Central carbon metabolism comparison to Toya et al for figure 3c
"""

from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle
import re

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from models.ecoli.analysis import multigenAnalysisPlot
import six
from six.moves import zip


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		cell_paths = self.ap.get_cells()

		validation_data = cPickle.load(open(validationDataFile, "rb"))
		toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
		toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
		toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
		toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

		sim_data = cPickle.load(open(simDataFile, 'rb'))

		modelFluxes = {}
		toyaOrder = []

		for rxn in toyaReactions:
			modelFluxes[rxn] = []
			toyaOrder.append(rxn)

		for simDir in cell_paths:
			simOutDir = os.path.join(simDir, "simOut")

			massListener = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = massListener.readColumn("cellMass")
			dryMass = massListener.readColumn("dryMass")
			coefficient = dryMass / cellMass * sim_data.constants.cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			compiled_reaction_ids = fbaResults.readAttribute("compiled_reaction_ids")
			compiled_reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fbaResults.readColumn("compiled_reaction_fluxes").T / coefficient).T

			for toyaReaction in toyaReactions:
				fluxTimeCourse = None

				for i, rxn in enumerate(compiled_reaction_ids):
					if re.findall(toyaReaction, rxn):
						if fluxTimeCourse is not None:
							fluxTimeCourse += compiled_reaction_fluxes[:, i]
						else:
							fluxTimeCourse = compiled_reaction_fluxes[:, i]

				if fluxTimeCourse is not None:
					modelFluxes[toyaReaction].append(np.mean(fluxTimeCourse).asNumber(units.mmol / units.g / units.h))

		toyaVsReactionAve = []
		for rxn, toyaFlux in six.viewitems(toyaFluxesDict):
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
