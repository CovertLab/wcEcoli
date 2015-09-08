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
CO2_ID = "CARBON-DIOXIDE[p]"
O2_ID = "OXYGEN-MOLECULE[p]"

FLUX_UNITS = units.mmol / units.L / units.s
MASS_UNITS = units.fg
COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L

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

	externalMoleculeIDs = fbaResults.readAttribute("externalMoleculeIDs")

	fbaResults.close()

	glucoseIdx = externalMoleculeIDs.index(GLUCOSE_ID)
	co2Index = externalMoleculeIDs.index(CO2_ID)
	o2Index = externalMoleculeIDs.index(O2_ID)
	glucoseFlux = -1. * FLUX_UNITS * externalExchangeFluxes[:, glucoseIdx] * kb.timeStepSec
	co2Flux = FLUX_UNITS * externalExchangeFluxes[:, co2Index] * kb.timeStepSec
	o2Flux = -1. * FLUX_UNITS * externalExchangeFluxes[:, o2Index] * kb.timeStepSec

	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = MASS_UNITS * mass.readColumn("cellMass")
	cellDryMass = MASS_UNITS * mass.readColumn("dryMass")
	dryMassGrowth = MASS_UNITS * mass.readColumn("growth") / (kb.timeStepSec * units.s)

	mass.close()

	cellDensity = kb.constants.cellDensity
	glucoseMW = kb.getter.getMass([GLUCOSE_ID])[0]
	cellVolume = cellMass / cellDensity

	glucoseMassFlux = glucoseFlux * glucoseMW * cellVolume

	glucoseMassYield = dryMassGrowth / glucoseMassFlux
	glucoseMassYield = units.convertNoUnitToNumber(glucoseMassYield)

	Co2GlucoseYield = (co2Flux / glucoseFlux) * (kb.getter.getMass(["CARBON-DIOXIDE"]) / kb.getter.getMass(["GLC"]))
	Co2GlucoseYield = units.convertNoUnitToNumber(Co2GlucoseYield)

	O2PerGlucose = o2Flux / glucoseFlux
	O2PerGlucose = units.convertNoUnitToNumber(O2PerGlucose)

	# Calculate constraints
	initWaterMass = kb.mass.avgCellWaterMassInit
	initDryMass = kb.mass.avgCellDryMassInit

	initCellMass = (
		initWaterMass
		+ initDryMass
		)

	coefficient = initDryMass / initCellMass * kb.constants.cellDensity * (kb.timeStepSec * units.s)

	externalMoleculeLevels = kb.process.metabolism.exchangeConstraints(
		externalMoleculeIDs,
		coefficient,
		COUNTS_UNITS / VOLUME_UNITS
		)

	fig, axesList = plt.subplots(5, sharex = True)
	fig.set_size_inches(8, 10)

	# Plot glucoes mas yield
	axesList[0].plot(time, glucoseMassYield)
	axesList[0].set_xlabel("Time (s)")
	axesList[0].set_ylabel("gDCW / g glucose")
	axesList[0].set_ylim([0., 1.])

	# Plot glucose exchange flux
	axesList[1].plot(time, (glucoseFlux / kb.timeStepSec).asNumber(FLUX_UNITS), label="Glucose exchange flux coefficient")
	axesList[1].plot(time, externalMoleculeLevels[glucoseIdx] * np.ones(time.shape), '--')
	axesList[1].set_ylabel("Glucose\nimport\n(mmol/L/s)")

	fluxesToPlot = np.where(np.abs(externalExchangeFluxes).max(axis=0) >= glucoseFlux.asNumber(FLUX_UNITS).max())[0]
	for idx in fluxesToPlot:
		moleculeId = externalMoleculeIDs[idx]
		flux = externalExchangeFluxes[:,idx]
		axesList[2].plot(time, flux, label=moleculeId)
	axesList[2].set_ylim([-6, 3])
	axesList[2].legend()
	axesList[2].set_ylabel('Exchange fluxes\n(mmol/L/s)')


	axesList[3].plot(time, Co2GlucoseYield)
	axesList[3].set_ylabel("g CO2 / g glucose")
	axesList[3].set_ylim([0., 1.])

	axesList[4].plot(time, O2PerGlucose)
	axesList[4].set_ylim([0., 5.])
	axesList[4].set_ylabel("mol O2 / mol glucose")

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
