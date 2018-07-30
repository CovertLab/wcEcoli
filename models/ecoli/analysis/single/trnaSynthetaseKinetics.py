"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/3/2017
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from wholecell.utils import units


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load kinetic data from validation
		validation_data = cPickle.load(open(validationDataFile, "rb"))
		kinetic_data = validation_data.trna.synthetaseKinetics
		synthetaseIds = kinetic_data.keys()

		# Order synthetases by their max kcat (in prep for plot)
		substrateList = ["tRNA", "AA", "ATP"]
		substrateColor = ["b", "g", "r"]
		maxKcats = []
		for synthetaseId in synthetaseIds:
			kcats = [kinetic_data[synthetaseId]["kcat %s" % x] for x in substrateList]
			maxKcats.append(max(kcats).asNumber())

		order = np.argsort(maxKcats)
		synthetaseIds = list(np.array(synthetaseIds)[order])

		# Load from sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		aaIds = sim_data.moleculeGroups.aaIDs
		timestep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

		# Load data from sim
		cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass") * units.fg
		cellVolume = cellMass / cellDensity

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		synthetaseIndexes = np.array([moleculeIds.index(x) for x in synthetaseIds], np.int)
		synthetaseCounts = bulkMolecules.readColumn("counts")[:, synthetaseIndexes]
		bulkMolecules.close()

		# Compute synthetase concentrations
		synthetaseConcentrations = (synthetaseCounts / nAvogadro) / cellVolume[:, None]

		# Load number of incorporations of each amino acid from sim
		aasUsed = TableReader(os.path.join(simOutDir, "GrowthLimits")).readColumn("aasUsed")
		selenocysteine_index = aaIds.index("L-SELENOCYSTEINE[c]")
		aaIds.remove("L-SELENOCYSTEINE[c]")
		aasUsed = np.delete(aasUsed, selenocysteine_index, 1)

		# Plot
		width = 0.25
		fig, axesList = plt.subplots(3, 1, figsize = (11, 8.5))
		ax1, ax2, ax3 = axesList
		x_ticklabels = []

		for i, synthetaseId in enumerate(synthetaseIds):
			aaId = kinetic_data[synthetaseId]["amino acid"]
			aaIndex = aaIds.index(aaId)
			x_ticklabels.append(aaId[:-3])
			E = np.average([x.asNumber() for x in synthetaseConcentrations[:, i]])
			ax2.bar(i + 1.5 * width, E * 1e6, width * 3, facecolor = "k", alpha = 0.5)

			# Calculate minimum kcat required to meet demands:
			# kcat required = (# aa used) / (E * timestep * volume * avogadro's number)
			denominator = [x.asNumber() for x in synthetaseConcentrations[:, i] * timestep * cellVolume * nAvogadro]
			kcat_required = np.average(aasUsed[:, aaIndex] / denominator)
			ax1.plot([i, i + width * 3], [kcat_required, kcat_required], color = "k", linewidth = 2)

			# Calculate synthetase concentration required to meet demands:
			# synthetase concentration required = (# aa used) / (kcat * timestep * volume * avogadro's number)
			kcat_max = np.max([kinetic_data[synthetaseId]["kcat %s" % x] for x in substrateList])
			denominator = timestep * units.s * cellVolume * nAvogadro * kcat_max
			E_required = np.average([x.asNumber() for x in aasUsed[:, aaIndex] / denominator]) * 1e6
			ax2.plot([i, i + width * 3], [E_required, E_required], color = "k", linewidth = 2)

			# Calculate average demand per timestep
			demand = np.average(aasUsed[:, aaIndex])
			ax3.bar(i + 1.5 * width, demand / 1000., width * 3, facecolor = "k", alpha = 0.5)

			for j, substrate in enumerate(substrateList):
				kcat = kinetic_data[synthetaseId]["kcat %s" % substrate]
				if kcat.asNumber() == 0:
					continue
				ax1.bar(i + j * width, kcat.asNumber(), width, facecolor = substrateColor[j])

		ax1.set_title("kcat (1/s)")
		ax2.set_title("average [synthetase] (uM)")
		ax3.set_title("average demand per timestep (1000s amino acids)")
		blue_patch = mpatches.Patch(color = "blue", label = "tRNA")
		green_patch = mpatches.Patch(color = "green", label = "amino acid")
		red_patch = mpatches.Patch(color = "red", label = "ATP")
		black_dashed_line = mlines.Line2D([], [], color = "k", linewidth = 2,
			label = "kcat required")

		plt.axes(ax1)
		plt.legend(
			handles = [blue_patch, green_patch, red_patch, black_dashed_line],
			loc = "upper left",
			bbox_to_anchor = [0, 1],
			)

		plt.axes(ax2)
		black_dashed_line = mlines.Line2D([], [], color = "k", linewidth = 2,
			label = "synthetase required (assumes max kcat value from literature)")
		plt.legend(
			handles = [black_dashed_line],
			loc = "upper right",
			)

		for ax in axesList:
			ax.set_xticks(np.arange(20) + 1.5 * width)
			ax.set_xticklabels(x_ticklabels, fontsize = 6)
			ax.tick_params(right = "off", top = "off", which = "both", direction = "out")
		plt.subplots_adjust(hspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
