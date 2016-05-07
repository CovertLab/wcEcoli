#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/2016
"""

from __future__ import division

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from scipy.stats import pearsonr

import mpld3
from mpld3 import plugins, utils

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, VOLUME_UNITS, TIME_UNITS

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

CMAP_COLORS_255 = [
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_OVER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	validation_data = cPickle.load(open(validationDataFile, "rb"))
	sim_data = cPickle.load(open(simDataFile, "rb"))

	cellDensity = sim_data.constants.cellDensity

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

	massListener = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = massListener.readColumn("cellMass") * units.fg
	dryMass = massListener.readColumn("dryMass") * units.fg
	massListener.close()

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
	reactionFluxes = (COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) * np.array(fbaResults.readColumn("reactionFluxes"))
	fluxes_dict = dict(zip(reactionIDs, reactionFluxes))
	fbaResults.close()

	coefficients = dryMass / cellMass * cellDensity * (timeStepSec * units.s)

	dryMassFracAverage = np.mean(dryMass / cellMass)

	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = FLUX_UNITS * np.array([(dryMassFracAverage * cellDensity * x).asNumber(FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]])

	# Time courses of only reactions for which we have toya fluxes
	trimmedReactions = reactionFluxes[:,np.array([np.nonzero(reactionIDs == toya_reaction)[0][0] for toya_reaction in toya_reactions])]

	corrCoefTimecourse = []
	for fluxes in trimmedReactions:
		correlationCoefficient = np.corrcoef(fluxes.asNumber(FLUX_UNITS), toya_fluxes.asNumber(FLUX_UNITS))[0,1]
		corrCoefTimecourse.append(correlationCoefficient)

	meanCorr = np.mean(np.array(corrCoefTimecourse)[~np.isnan(corrCoefTimecourse)])

	fig = plt.figure()
	plt.plot(time / 60., corrCoefTimecourse)
	plt.axhline(y=meanCorr, color='r')
	plt.title("Measured vs. Simulated Central Carbon Fluxes")
	plt.text(.5*np.max(time / 60.),1.05*meanCorr, "Mean = {:.2}".format(meanCorr), horizontalalignment="center")
	plt.xlabel("Time (min)")
	plt.ylabel("Pearson R")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
	
