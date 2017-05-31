#!/usr/bin/env python
"""
Plots fold change in cell cycle length from control for all simulations in all variants.

@author: Nicole Ferraro
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/31/2017
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	varDirs = os.listdir(seedOutDir)

	control_length = 0
	cellCycleLengths = []

	#Currently does not do anything with dead cell list, but collects if interested in future
	dead_cells = []

	for variant in varDirs:
		variantDir = os.path.join(seedOutDir, variant)
		ap = AnalysisPaths(variantDir, cohort_plot = True)
		allDir = ap.get_cells()
		for idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			try:
				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			except:
				dead_cells.append(variant)
				continue

			cellCycle = (time[-1] - time[0]) / 60. / 60.

			#Assumes first variant is the control condition
			if '000000' in variant:
				control_length = cellCycle
			else:
				#Make sure cell actually divided and did not just reach error or numerical instability
				if os.path.exists(os.path.join(simOutDir, "Daughter1")):
					cellCycleLengths.append(cellCycle)
				else:
					dead_cells.append(variant)
	
	cellCycleLengths.sort()
	cellCycleFC = cellCycleLengths / control_length

	plt.bar(range(0,len(cellCycleFC)), cellCycleFC)
	plt.xlabel('Individual simulations across variants, sorted')
	plt.ylabel('Fold change in cell cycle length over control')
	plt.title('Cell cycle fold change across variants')
	plt.gca().xaxis.set_major_locator(plt.NullLocator())
	plt.axhline(y=1.02, linewidth=4, color='r')
	plt.text(len(cellCycleFC)/4, 1.1, 'Fold change of 1.02', fontsize=12)

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

	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	#main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])