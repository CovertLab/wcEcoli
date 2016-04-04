#!/usr/bin/env python
"""
@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 30/03/2016
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
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

FONT = {
		'size'	:	12
		}

FLUX_UNITS = "M/s"

CMAP_COLORS_255 = [
	[103,0,31],
	[178,24,43],
	[214,96,77],
	[244,165,130],
	[253,219,199],
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]
CMAP_UNDER = [1, 0.2, 0.75]
CMAP_OVER = [0, 1, 0.75]

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	PCC = fbaResults.readColumn("PCC")
	PCClog = fbaResults.readColumn("PCClog")

	slope = fbaResults.readColumn("slope")
	slopelog = fbaResults.readColumn("slopelog")

	fbaResults.close()

	# compute correlation between observed and predicted metabolic fluxes	
	plt.figure(figsize = (8.5, 8.5))
	matplotlib.rc('font', **FONT)
	max_yticks = 5

	ax = plt.subplot(4,1,1)
	plt.plot(time / 60., PCC)
	plt.plot(time / 60., np.ones(len(time)) * np.average(PCC), '--k')
	PCCfba = 0.8124
	plt.plot(time / 60., np.ones(len(time)) * PCCfba, '-r')
	plt.ylabel("PCC")
	plt.title("Comparison observed vs. predicted metabolic fluxes\n\n PCC_avg = %.2f, PCC_sd = %.2f" % (np.average(PCC), np.std(PCC)), fontsize = 12)
	yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

	ax = plt.subplot(4,1,2)
	plt.plot(time / 60., PCClog)
	plt.plot(time / 60., np.ones(len(time)) * np.average(PCClog), '--k')
	PCCfbaLog = 0.8185
	plt.plot(time / 60., np.ones(len(time)) * PCCfbaLog, '-r')
	plt.ylabel("PCC, Log10 (.)")	
	plt.title("PCC_avg = %.2f, PCC_sd = %.2f" % (np.average(PCClog), np.std(PCClog)), fontsize = 12)
	yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

	ax = plt.subplot(4,1,3)
	plt.plot(time / 60., slope)
	plt.plot(time / 60., np.ones(len(time)) * np.average(slope), '--k')
	plt.plot(time / 60., np.ones(len(time)), '-r')
	plt.ylabel("Slope")
	plt.title("Slope_avg = %.2f, Slope_sd = %.2f" % (np.average(slope), np.std(slope)), fontsize = 12)
	yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

	ax = plt.subplot(4,1,4)
	plt.plot(time / 60., slopelog)
	plt.plot(time / 60., np.ones(len(time)) * np.average(slopelog), '--k')
	plt.plot(time / 60., np.ones(len(time)), '-r')
	plt.xlabel("Time (min)")
	plt.ylabel("Slope Log (.)")
	plt.title("SlopeLog10_avg = %.2f, SlopeLog10_sd = %.2f" % (np.average(slopelog), np.std(slopelog)), fontsize = 12)
	yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

	plt.subplots_adjust(hspace = 0.4, wspace = 0.2)

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

