#!/usr/bin/env python
"""
Plots rRNA control loop parameters

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/17/2017
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

PLACE_HOLDER = -1

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load KB
	sim_data = cPickle.load(open(simDataFile, "rb"))

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load control loop data
	errorInElongationRate = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("errorInElongationRate")
	proportionalTerm = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("proportionalTerm")
	integralTerm = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("integralTerm")
	rRnaSynthRate_expected = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_expected")
	rRnaSynthRate_updated = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_updated")

	growth_rate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
	growthRate = (1 / units.s) * growthRate
	doublingTime = 1 / growthRate * np.log(2)

	effective_elongation_rate = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate")


	# nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel

	# sim_data.process.translation.ribosomeElongationRateDict[self._sim.processes["PolypeptideElongation"].currentNutrients]
	# expected_elongation_rate = 



	# Plot stuff
	plt.figure(figsize = (8.5, 11))

	fig, axesList = plt.subplots(5,1, sharex = True)

	axesList[0].plot(time / 60., errorInElongationRate)

	axesList[1].plot(time / 60., proportionalTerm)
	axesList[1].plot(time / 60., integralTerm)
	axesList[1].plot(time / 60., proportionalTerm + integralTerm)

	axesList[2].plot(time / 60., rRnaSynthRate_expected)
	axesList[2].plot(time / 60., rRnaSynthRate_updated)

	axesList[3].plot(time / 60., doublingTime.asNumber(units.min))

	axesList[4].plot(time / 60., effective_elongation_rate)

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
