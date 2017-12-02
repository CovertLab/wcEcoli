#!/usr/bin/env python
"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/1/2017
"""

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
from wholecell.utils.sparkline import whitePadSparklineAxis


def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"
	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load from sim data
	sim_data = cPickle.load(open(simDataFile, "rb"))
	aaIds = sim_data.moleculeGroups.aaIDs
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

	# Load from growth limits listener
	growthLimitsDataFile = TableReader(os.path.join(simOutDir, "GrowthLimits"))

	aaRequestedByTranslationSupply = growthLimitsDataFile.readColumn("aaRequestedByTranslationSupply")
	aaRequestedBySynthetaseKineticCapacity = growthLimitsDataFile.readColumn("aaRequestedBySynthetaseKineticCapacity")

	# Convention: pos means synthetase kinetic capacity is sufficient
	diff = aaRequestedBySynthetaseKineticCapacity - aaRequestedByTranslationSupply

	# Plot
	fig, axesList = plt.subplots(3, 1, figsize = (10, 10))
	ax1, ax2, ax3 = axesList

	limiting = ["ASN[c]", "TRP[c]", "SER[c]", "GLN[c]", "PHE[c]", "LYS[c]", "GLT[c]", "LEU[c]", "MET[c]", "L-ASPARTATE[c]"]

	from wholecell.analysis.analysis_tools import exportFigure
	for i in xrange(diff.shape[1]):
		if aaIds[i] == "L-SELENOCYSTEINE[c]":
			continue
		if aaIds[i] in limiting:
			tag = "*"
		else:
			tag = ""
		if np.all(diff[20:, i] > 0):
			ax1.plot(time / 60., diff[:, i], label = "%s%s" % (aaIds[i], tag))
		elif np.all(diff[20:, i] < 0):
			ax2.plot(time / 60., diff[:, i], label = "%s%s" % (aaIds[i], tag))
		else:
			ax3.plot(time / 60., diff[:, i], label = "%s%s" % (aaIds[i], tag))
		# import ipdb; ipdb.set_trace()

	for ax in axesList:
		ax.tick_params(right = "off", top = "off", which = "both", direction = "out")
		ax.set_xlim([time[0] / 60., time[-1] / 60.])
		# ax.set_ylim(np.round(diff.min(), -3), np.round(diff.max(), -2))

		plt.axes(ax)
		plt.legend(loc = 1, fontsize = 6)


	ax3.set_xlabel("Time (min)")
	ax1.set_title("Supply - Demand (amino acid counts)")
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.85, right = 0.95)

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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
