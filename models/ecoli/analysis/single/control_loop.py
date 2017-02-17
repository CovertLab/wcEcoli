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

FONT_SIZE = 7

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
	bias = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_updated")

	growthRate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
	growthRate = (1 / units.s) * growthRate
	doublingTime = 1 / growthRate * np.log(2)

	effective_elongation_rate = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate")

	nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
	sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][0][1]
	media = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][0][1]
	expected_elongation_rate = sim_data.process.translation.ribosomeElongationRateDict[media]
	expected_doubling_time = sim_data.nutrientToDoublingTime[media]

	cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
	uniqueMoleculeCounts.close()
	ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass))

	# Plot stuff
	plt.figure(figsize = (8.5, 11))

	fig, axesList = plt.subplots(6,1, sharex = True)

	axesList[0].plot(time / 60., errorInElongationRate, label="raw")
	axesList[0].plot(time / 60., np.zeros(time.size), '--')
	axesList[0].plot(time / 60., errorInElongationRate + bias, label="bias")
	axesList[0].legend(fontsize=FONT_SIZE, loc=4,frameon=False)
	axesList[0].set_ylabel("Error " + r"$(e_{expected} - e_{actual})$", fontsize=FONT_SIZE)

	axesList[1].plot(time / 60., proportionalTerm, label = "proportional", alpha = 0.7)
	axesList[1].plot(time / 60., integralTerm, label = "integral", alpha = 0.7)
	axesList[1].plot(time / 60., proportionalTerm + integralTerm, label = "sum", alpha = 0.7)
	axesList[1].legend(fontsize=FONT_SIZE, loc=4,frameon=False)
	axesList[1].set_ylabel("Correction terms", fontsize=FONT_SIZE)

	axesList[2].plot(time / 60., rRnaSynthRate_expected, label = "expected")
	axesList[2].plot(time / 60., rRnaSynthRate_updated, label = "updated")
	axesList[2].legend(fontsize=FONT_SIZE, loc=4,frameon=False)
	axesList[2].set_ylabel("rRNA synthesis prob", fontsize=FONT_SIZE)

	axesList[3].plot(time / 60., doublingTime.asNumber(units.min))
	axesList[3].set_ylabel("Inst. doubling\ntime (min)", fontsize=FONT_SIZE)
	axesList[3].plot(time / 60., expected_doubling_time.asNumber(units.min) * np.ones(time.size), '--')

	axesList[4].plot(time / 60., effective_elongation_rate)
	axesList[4].set_ylabel("Eff. elng.\nrate", fontsize=FONT_SIZE)
	axesList[4].plot(time / 60., expected_elongation_rate.asNumber() * np.ones(time.size), '--')

	axesList[5].plot(time / 60., ribosomeConcentration.asNumber(units.mmol / units.L))
	axesList[5].set_ylabel("[Rib]", fontsize=FONT_SIZE)

	for a in axesList:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in a.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		# a.spines["right"].set_visible(False)
		a.spines["top"].set_visible(False)
		a.spines["bottom"].set_visible(False)
		a.xaxis.set_visible(False)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

	raise Exception()

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
