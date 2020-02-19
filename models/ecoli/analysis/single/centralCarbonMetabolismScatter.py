"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/2016
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle
import re

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units

from wholecell.analysis.plotting_tools import CMAP_COLORS_255

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		validation_data = cPickle.load(open(validationDataFile, "rb"))
		sim_data = cPickle.load(open(simDataFile, "rb"))

		cellDensity = sim_data.constants.cellDensity

		massListener = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massListener.readColumn("cellMass") * units.fg
		dryMass = massListener.readColumn("dryMass") * units.fg
		massListener.close()

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
		reactionFluxes = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * np.array(fbaResults.readColumn("reactionFluxes"))
		fbaResults.close()

		dryMassFracAverage = np.mean(dryMass / cellMass)

		#Toya Graph
		toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
		toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])
		toya_stdev = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]])
		toya_fluxes_dict = dict(zip(toya_reactions, toya_fluxes))
		toya_stdev_dict = dict(zip(toya_reactions, toya_stdev))

		toyaVsReactionAve = []
		toya_order = []
		for toyaReactionID, toyaFlux in toya_fluxes_dict.iteritems():
			fluxTimeCourse = []

			for rxn in reactionIDs:
				if toyaReactionID in rxn:
				#if re.findall(toyaReactionID, rxn):
					reverse = 1
					if "(reverse)" in rxn:
					#if re.findall("(reverse)", rxn):
						reverse = -1

					if len(fluxTimeCourse):
						fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
					else:
						fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

			if len(fluxTimeCourse):
				fluxAve = np.mean(fluxTimeCourse)
				fluxStdev = np.std(fluxTimeCourse.asNumber(FLUX_UNITS))
				toyaVsReactionAve.append((fluxAve.asNumber(FLUX_UNITS), toyaFlux.asNumber(FLUX_UNITS), fluxStdev, toya_stdev_dict[toyaReactionID].asNumber(FLUX_UNITS)))
				toya_order.append(toyaReactionID)

		toyaVsReactionAve = FLUX_UNITS * np.array(toyaVsReactionAve)
		correlationCoefficient = np.corrcoef(toyaVsReactionAve[:,0].asNumber(FLUX_UNITS), toyaVsReactionAve[:,1].asNumber(FLUX_UNITS))[0,1]

		plt.figure()

		plt.title("Central Carbon Metabolism Flux, Pearson R = {:.2}".format(correlationCoefficient))
		plt.errorbar(toyaVsReactionAve[:,1].asNumber(FLUX_UNITS), toyaVsReactionAve[:,0].asNumber(FLUX_UNITS), xerr = toyaVsReactionAve[:,3].asNumber(FLUX_UNITS),
			yerr = toyaVsReactionAve[:,2].asNumber(FLUX_UNITS), fmt = "o", ecolor = "k")
		plt.plot(toyaVsReactionAve[:, 1].asNumber(FLUX_UNITS), toyaVsReactionAve[:, 0].asNumber(FLUX_UNITS), "o")
		for x, y, label in zip(toyaVsReactionAve[:, 1].asNumber(FLUX_UNITS), toyaVsReactionAve[:, 0].asNumber(FLUX_UNITS), toya_order):
			plt.text(x, y, label)
		plt.ylabel("Mean WCM Reaction Flux {}".format(FLUX_UNITS.strUnit()))
		plt.xlabel("Toya 2010 Reaction Flux {}".format(FLUX_UNITS.strUnit()))

		exportFigure(plt, plotOutDir, plotOutFileName + "_toya", metadata)
		plt.close("all")

		#Long  Graph

		long_reactions = validation_data.reactionFlux.long2019fluxes2010fluxes["reactionID"]
		long_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in
											 validation_data.reactionFlux.long2019fluxes["reactionFlux"]])
		long_range = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in
											validation_data.reactionFlux.long2019fluxes["reactionFluxRange"]])
		long_fluxes_dict = dict(zip(long_reactions, long_fluxes))
		long_range_dict = dict(zip(long_reactions, long_range))

		longVsReactionAve = []
		long_order = []
		for longReactionID, longFlux in long_fluxes_dict.iteritems():
			fluxTimeCourse = []

			for rxn in reactionIDs:
				if longReactionID in rxn:
					# if re.findall(longReactionID, rxn):
					reverse = 1
					if "(reverse)" in rxn:
						# if re.findall("(reverse)", rxn):
						reverse = -1

					if len(fluxTimeCourse):
						fluxTimeCourse += reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]
					else:
						fluxTimeCourse = reverse * reactionFluxes[:, np.where(reactionIDs == rxn)]

			if len(fluxTimeCourse):
				fluxAve = np.mean(fluxTimeCourse)
				fluxRange = np.std(fluxTimeCourse.asNumber(FLUX_UNITS))
				longVsReactionAve.append((fluxAve.asNumber(FLUX_UNITS), longFlux.asNumber(FLUX_UNITS), fluxRange,
										  long_range_dict[longReactionID].asNumber(FLUX_UNITS)))
				long_order.append(longReactionID)

		longVsReactionAve = FLUX_UNITS * np.array(longVsReactionAve)
		correlationCoefficient = \
		np.corrcoef(longVsReactionAve[:, 0].asNumber(FLUX_UNITS), longVsReactionAve[:, 1].asNumber(FLUX_UNITS))[0, 1]

		plt.figure()

		plt.title("Central Carbon Metabolism Flux, Pearson R = {:.2}".format(correlationCoefficient))
		# plt.errorbar(longVsReactionAve[:,1].asNumber(FLUX_UNITS), longVsReactionAve[:,0].asNumber(FLUX_UNITS), xerr = longVsReactionAve[:,3].asNumber(FLUX_UNITS),
		# yerr = longVsReactionAve[:,2].asNumber(FLUX_UNITS), fmt = "o", ecolor = "k")
		plt.plot(longVsReactionAve[:, 1].asNumber(FLUX_UNITS), longVsReactionAve[:, 0].asNumber(FLUX_UNITS), "o")
		for x, y, label in zip(longVsReactionAve[:, 1].asNumber(FLUX_UNITS),
							   longVsReactionAve[:, 0].asNumber(FLUX_UNITS), long_order):
			plt.text(x, y, label)
		plt.ylabel("Mean WCM Reaction Flux {}".format(FLUX_UNITS.strUnit()))
		plt.xlabel("Long 2019 Reaction Flux {}".format(FLUX_UNITS.strUnit()))

		exportFigure(plt, plotOutDir, plotOutFileName + "_long", metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
