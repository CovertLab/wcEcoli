"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/1/2017
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
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
		fig, axesList = plt.subplots(3, 1, figsize = (8.5, 11))
		ax1, ax2, ax3 = axesList
		limiting = ["ASN[c]", "TRP[c]", "SER[c]", "GLN[c]", "PHE[c]", "LYS[c]", "GLT[c]", "LEU[c]", "MET[c]", "L-ASPARTATE[c]"]

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

		for ax in axesList:
			ax.tick_params(right = "off", top = "off", which = "both", direction = "out")
			ax.set_xlim([time[0] / 60., time[-1] / 60.])
			plt.axes(ax)
			plt.legend(loc = 1, fontsize = 6)

		ax3.set_xlabel("Time (min)")
		ax1.set_title("Supply - Demand (amino acid counts)")
		plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.85, right = 0.95)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
