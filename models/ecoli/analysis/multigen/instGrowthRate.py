#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))
	ap = AnalysisPaths(variantDir, multi_gen_plot = True)

	fig, axesList = plt.subplots(5)

	for simDir in ap.get_cells():
		simOutDir = os.path.join(simDir, "simOut")
		
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
		massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
		growthRate = massDataFile.readColumn("instantaniousGrowthRate")

		rnaMass = massDataFile.readColumn("rnaMass")
		proteinMass = massDataFile.readColumn("proteinMass")
		cellMass = massDataFile.readColumn("cellMass")

		rnaGrowthRate = np.hstack((0, np.diff(rnaMass))) / rnaMass / timeStepSec
		proteinGrowthRate = np.hstack((0, np.diff(proteinMass))) / proteinMass / timeStepSec

		cellMassAdded = np.hstack((0, np.diff(cellMass)))
		rnaMassAdded = np.hstack((0, np.diff(rnaMass)))
		proteinMassAdded = np.hstack((0, np.diff(proteinMass)))

		fracProteinMassAdded = proteinMassAdded / cellMassAdded
		fracRnaMassAdded = rnaMassAdded / cellMassAdded

		width = 100
		if time.size < width:
			width = time.size
		# axesList[0].plot(time / 60., growthRate)
		# axesList[1].plot(time / 60., rnaGrowthRate / growthRate)
		# axesList[2].plot(time / 60., proteinGrowthRate / growthRate)

		axesList[0].plot(time / 60., np.convolve(growthRate, np.ones(width) / width, mode = "same"))
		axesList[1].plot(time / 60., np.convolve(rnaGrowthRate / growthRate, np.ones(width) / width, mode = "same"))
		axesList[2].plot(time / 60., np.convolve(proteinGrowthRate / growthRate, np.ones(width) / width, mode = "same"))
		axesList[3].plot(time / 60., np.convolve(fracProteinMassAdded, np.ones(width) / width, mode = "same"))
		axesList[4].plot(time / 60., np.convolve(fracRnaMassAdded, np.ones(width) / width, mode = "same"))

	expectedGrowthRate = np.log(2) / sim_data.doubling_time.asNumber(units.s)
	bound_offset = 0.7
	axesList[0].set_ylim([expectedGrowthRate * (1-bound_offset), expectedGrowthRate*(1+bound_offset)])
	axesList[1].set_ylim([1.0 * (1-bound_offset), 1.0*(1+bound_offset)])
	axesList[2].set_ylim([1.0 * (1-bound_offset), 1.0*(1+bound_offset)])
	axesList[3].set_ylim([0., 0.2])
	axesList[4].set_ylim([0., 0.2])

	axesList[4].set_xlabel("Time (min)")

	axesList[0].set_ylabel("Instantanious\ngrowth rate (1/s)")
	axesList[1].set_ylabel("RNA/total\ngrowth rates")
	axesList[2].set_ylabel("Protein/total\ngrowth rates")
	axesList[3].set_ylabel("Frac. protein\nmass added")
	axesList[4].set_ylabel("Frac. rna\nmass added")

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
