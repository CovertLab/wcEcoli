#!/usr/bin/env python
"""
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/8/2014
"""

from __future__ import division

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

GLUCOSE_ID = "GLC[p]"

FLUX_UNITS = units.mmol / units.L / units.s
MASS_UNITS = units.fg
GROWTH_UNITS = units.fg / units.s

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	kb = cPickle.load(open(kbFile, "rb"))

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStep")
	externalExchangeFluxes = fbaResults.readColumn("externalExchangeFluxes")

	externalMoleculeIDs = np.array(fbaResults.readAttribute("externalMoleculeIDs"))

	fbaResults.close()

	glucoseIdx = np.where(externalMoleculeIDs == GLUCOSE_ID)[0][0]
	glucoseFlux = -1. * FLUX_UNITS * externalExchangeFluxes[:, glucoseIdx] * kb.timeStepSec

	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = MASS_UNITS * mass.readColumn("cellMass")
	cellDryMass = MASS_UNITS * mass.readColumn("dryMass")
	growth = GROWTH_UNITS * mass.readColumn("growth") * kb.timeStepSec

	mass.close()

	cellDensity = kb.constants.cellDensity
	glucoseMW = kb.getter.getMass([GLUCOSE_ID])[0]

	glucoseMassFlux = glucoseFlux * glucoseMW * cellDryMass / cellDensity

	glucoseMassYield = growth / glucoseMassFlux
	glucoseMassYield = units.convertNoUnitToNumber(glucoseMassYield)

	fig, axesList = plt.subplots(3, sharex = True)
	fig.set_size_inches(8, 8)

	# Plot glucoes mas yield
	axesList[0].plot(time, glucoseMassYield)
	axesList[0].set_xlabel("Time (s)")
	axesList[0].set_ylabel("gDCW / g glucose")
	axesList[0].set_ylim([0., 1.])


	# Plot glucose exchange flux
	# fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
	# exFlux = fba_results.readColumn("externalExchangeFluxes")
	# exMolec = fba_results.readAttribute("externalMoleculeIDs")
	# glcFlux = exFlux[:,exMolec.index("GLC[p]")]

	# axesList[-1].plot(time / 60. / 60., -1. * glcFlux, label="Glucose exchange flux coefficient", color=COLORS[sim_idx])
	# axesList[-1].set_ylabel("External\nglucose\n(mmol/L/s)")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
