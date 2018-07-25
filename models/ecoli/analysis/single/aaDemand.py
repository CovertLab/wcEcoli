"""
Plot distribution of amino acid demand compared to tRNA supply

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/30/2018
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


FONTSIZE = 8

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get amino acid usage
		growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))
		aasUsed = growthLimits.readColumn("aasUsed")
		growthLimits.close()

		# Get tRNA counts
		sim_data = cPickle.load(open(simDataFile))
		isTRna = sim_data.process.transcription.rnaData["isTRna"]
		rnaIds = sim_data.process.transcription.rnaData["id"][isTRna]
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
		rnaCounts = bulkMolecules.readColumn("counts")[:, rnaIndexes]
		bulkMolecules.close()

		# Map tRNAs to amino acids
		trnaToAaMapping = sim_data.relation.TRnaToAaMapping
		rnaCounts_summedByAa = np.dot(rnaCounts, trnaToAaMapping)

		# Plot
		fig, axesList = plt.subplots(6, 4, figsize = (8.5, 11))
		axToRemove = axesList.flatten()[-3:]
		axesList = axesList.flatten()[:-3]
		aaIds = sim_data.moleculeGroups.aaIDs
		for i, ax in enumerate(axesList):
			# Normalize data (to get fraction of timesteps)
			data, bins = np.histogram(aasUsed[:, i] / rnaCounts_summedByAa[:, i], bins = 20)
			data = data.astype(np.float32) / data.sum()

			ax.bar(bins[:-1], data, width = (bins[1] - bins[0]))
			ax.set_title(aaIds[i], fontsize = FONTSIZE)
			ax.set_ylabel("Fraction of simulation", fontsize = FONTSIZE)
			ax.set_yticks([0, 0.5, 1])
			ax.tick_params(axis = "both", labelsize = FONTSIZE)

			# Report mean average
			ax.text(0.5, 0.8,
				"mean = %0.2f" % np.mean(aasUsed[:, i] / rnaCounts_summedByAa[:, i]),
				transform = ax.transAxes,
				ha = "center",
				fontsize = FONTSIZE,
				)

		# Remove extra axes
		for ax in axToRemove:
			ax.axis("off")
		plt.suptitle("Ratio of amino acids used to its tRNA count")
		plt.subplots_adjust(wspace = 0.5, hspace = 0.5, bottom = 0.07, top = 0.92)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
