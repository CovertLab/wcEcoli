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
	ap = AnalysisPaths(variantDir, cohort_plot = True)

	fig, axesList = plt.subplots(2,1)

	for simDir in ap.get_cells():
		simOutDir = os.path.join(simDir, "simOut")
		# Exchange flux
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		exFlux = fba_results.readColumn("externalExchangeFluxes")
		exMolec = fba_results.readAttribute("externalMoleculeIDs")
		moleculeIDs = ["GLC[p]", "OXYGEN-MOLECULE[p]"]
		
		# Plot
		for index, molecule in enumerate(["GLC[p]", "OXYGEN-MOLECULE[p]"]):
			if molecule not in exMolec:
				continue
			moleculeFlux = -1. * exFlux[:, exMolec.index(molecule)]
			axesList[index].plot(time / 60. / 60., moleculeFlux)

	axesList[0].set_ylim([0, 15])
	axesList[1].set_ylim([0, 50])


	# averageFlux = np.average(moleculeFlux)
	# yRange = np.min([np.abs(np.max(moleculeFlux) - averageFlux), np.abs(np.min(moleculeFlux) - averageFlux)])
	# ymin = np.round(averageFlux - yRange)
	# ymax = np.round(averageFlux + yRange)
	# axesList[index].set_ylim([ymin, ymax])



	axesList[0].set_ylabel("GLC[p] flux\n(mmol/gDCW/hr)")
	axesList[1].set_ylabel("O2[p] flux\n(mmol/gDCW/hr)")
	axesList[1].set_xlabel("Time (min)")

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
