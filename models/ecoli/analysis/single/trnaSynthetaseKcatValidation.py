"""
Plots tRNA synthetase kinetics validation.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/2018
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from matplotlib.patches import Rectangle
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
		trnaValidation = validation_data.trna.synthetaseKineticsValidation
		aaToSynthetase = validation_data.trna.aaToSynthetase
		kinetic_data = validation_data.trna.synthetaseKinetics
		synthetaseIds = kinetic_data.keys()

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		timestep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass") * units.fg
		cellVolume = cellMass / cellDensity
		aaOrdered = sim_data.amino_acid_1_to_3_ordered
		aaIds = sim_data.moleculeGroups.aaIDs

		# Load counts of synthetases
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		synthetaseIndexes = np.array([moleculeIds.index(x) for x in synthetaseIds], np.int)
		synthetaseCounts = bulkMolecules.readColumn("counts")[:, synthetaseIndexes]
		bulkMolecules.close()

		# Load number of amino acid incorporations
		aasUsed = TableReader(os.path.join(simOutDir, "GrowthLimits")).readColumn("aasUsed")

		# Pop selenocysteine from aaOrdered
		aaOrdered.pop("U")

		# Plot kinetic values from literature
		width = 0.8
		fig, ax = plt.subplots(1, 1, figsize = (11, 8.5))

		for i, aa in enumerate(aaOrdered.keys()):
			aaId = aaOrdered[aa]
			synthetaseId = aaToSynthetase[aaId]
			kcats = np.array([x["kcat"].asNumber() for x in trnaValidation[synthetaseId[:-3]]])
			temp37 = np.array([x["Temperature"] == 37 for x in trnaValidation[synthetaseId[:-3]]])

			# Remove 0 entries
			mask = [x != 0 for x in kcats]
			if False in mask:
				kcats = kcats[mask]
				temp37 = temp37[mask]
			kcat_min = min(kcats)
			kcat_max = max(kcats)

			# Draw range of kcat values from literature
			p = Rectangle((i - (0.5 * width), kcat_min), width,
				kcat_max - kcat_min, alpha = 0.25, fc = "b")
			ax.add_patch(p)

			# Draw each kcat value
			for j in np.where(temp37 == False)[0]:
				ax.hlines(kcats[j], i - (0.5 * width), i + (0.5 * width), lw = 0.5, color = "r")
			
			for j in np.where(temp37)[0]:
				ax.hlines(kcats[j], i - (0.5 * width), i + (0.5 * width), lw = 0.5, color = "k")

			# Indicate kcat value previous used
			kcat = max([kinetic_data[synthetaseId][x] for x in ["kcat AA", "kcat ATP", "kcat tRNA"]]).asNumber()
			ax.scatter(i, kcat, facecolor = "none", edgecolor = "b", s = 10)

			# Plot required kcat computed from sim
			aaIndex = aaIds.index(aaId)
			aaDemand = aasUsed[:, aaIndex]
			demand = aaDemand / nAvogadro / cellVolume / (timestep * units.s)

			synthetaseIndex = synthetaseIds.index(synthetaseId)
			synthetaseSupply = synthetaseCounts[:, synthetaseIndex]
			supply = synthetaseSupply / nAvogadro / cellVolume

			ratio = [x.asNumber() for x in (demand / supply)]
			kcatRequiredMax = max(ratio)
			kcatRequiredAvg = np.average(ratio)
			ax.scatter(i, kcatRequiredMax, c = "b")
			ax.scatter(i, kcatRequiredAvg, facecolor = "none", edgecolor = "b")

		ax.set_xticks(range(len(aaOrdered)))
		ax.set_xticklabels(aaOrdered.keys())
		ax.set_xlabel("Amino Acid")
		ax.set_ylabel("Synthetase kcat (1/s)")
		ax.set_yscale("log")

		patch = mpatches.Patch(color = "b", alpha = 0.25,
			label = "Range of literature values")
		kLine = mlines.Line2D([], [], color = "k", label = "Temperature 37")
		rLine = mlines.Line2D([], [], color = "r",
			label = "Temperature <37 or not reported")
		bFilled = mlines.Line2D([], [], marker = "o", color = "b",
			linestyle = "None", label = "Max required kcat")
		bEdge = mlines.Line2D([], [], marker = "o", color = "b",
			linestyle = "None", markerfacecolor = "None",
			label = "Avg required kcat")
		bEdgeS = mlines.Line2D([], [], marker = "o", markersize = 3,
			color = "b", linestyle = "None", markerfacecolor = "None",
			label = "Previously modeled kcat")
		plt.legend(handles = [patch, kLine, rLine, bFilled, bEdge, bEdgeS],
			prop = {'size': 8})

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
