#!/usr/bin/env python
"""
@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 30/03/2016
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	# Get all cells
	allDir = ap.get_cells()

	cleanNames = [
				"PCC",
				"PCC Log10",
				"Slope",
				"Slope Log10",
				]

	#plt.figure(figsize = (8.5, 11))
	fig, axesList = plt.subplots(len(cleanNames), sharex = True)

	currentMaxTime = 0

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))

		for idx, auxType in enumerate(cleanNames):

			PCC = fbaResults.readColumn("PCC")
			PCClog = fbaResults.readColumn("PCClog")

			slope = fbaResults.readColumn("slope")
			slopelog = fbaResults.readColumn("slopelog")
			if idx == 0: aux = PCC; aux1 = 0.8124
			if idx == 1: aux = PCClog; aux1 = 0.8185
			if idx == 2: aux = slope; aux1 = 1
			if idx == 3: aux = slopelog; aux1 = 1

			axesList[idx].semilogy(time / 60. / 60., aux)
			axesList[idx].semilogy(time / 60., np.ones(len(time)) * aux1, '--k')
			
			# set axes to size that shows all generations
			cellCycleTime = ((time[-1] - time[0]) / 60. / 60. )
			if cellCycleTime > currentMaxTime:
				currentMaxTime = cellCycleTime

			total_gens = 10; axesList[idx].set_xlim(0, currentMaxTime*int(total_gens)*1.1)
			# axesList[idx].set_xlim(0, currentMaxTime*int(metadata["total_gens"])*1.1)
			
			axesList[idx].set_ylabel(cleanNames[idx])

	for axes in axesList:
		axes.get_ylim()
		axes.set_yticks(list(axes.get_ylim()))

	axesList[0].set_title("Validation metabolic fluxes")
	axesList[len(cleanNames) - 1].set_xlabel("Time (hr)")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName,metadata)
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

	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	#main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
