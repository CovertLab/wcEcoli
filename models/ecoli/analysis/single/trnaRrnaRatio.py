"""
Plot ratio of trna to 16S rRNA

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/31/2018
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units


FONTSIZE = 8

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get counts of tRNAs and 16S rRNA
		sim_data = cPickle.load(open(simDataFile))
		isTRna = sim_data.process.transcription.rnaData["isTRna"]
		trnaIds = sim_data.process.transcription.rnaData["id"][isTRna]
		rrnaIds = ["RRSA-RRNA[c]", "RRSB-RRNA[c]", "RRSC-RRNA[c]", "RRSD-RRNA[c]", "RRSE-RRNA[c]", "RRSG-RRNA[c]", "RRSH-RRNA[c]"]

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		trnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in trnaIds], np.int)
		rrnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rrnaIds], np.int)
		trnaCounts = bulkMolecules.readColumn("counts")[:, trnaIndexes]
		rrnaCounts = bulkMolecules.readColumn("counts")[:, rrnaIndexes]
		bulkMolecules.close()

		rrnaCounts_sum = rrnaCounts.sum(axis = 1, dtype = np.float)

		# Get known tRNA to 16S rRNA ratio
		validation_data = cPickle.load(open(validationDataFile))
		trnaTo16SRatio = validation_data.trna.trnaTo16SRatio

		# Identify closest growth rate
		growthRates = [0.4, 0.7, 1.07, 1.6, 2.5] # units of doublings/hr
		sim_growthRate = 1./sim_data.doubling_time.asNumber(units.h)
		closestGrowthRate = growthRates[np.argmin([abs(x - sim_growthRate) for x in growthRates])]

		# Plot
		fig, ax = plt.subplots(1, 1, figsize = (8.5, 11))

		for i, trnaId in enumerate(trnaIds):
			trnaToRranRatio = trnaCounts[:, i] / rrnaCounts_sum
			expectedRatio = trnaTo16SRatio[trnaId]["%s per hr" % closestGrowthRate]
			ax.scatter(np.mean(trnaToRranRatio), expectedRatio, c = "b")
			ax.scatter(np.mean(trnaToRranRatio), trnaTo16SRatio[trnaId]["%s per hr" % 1.07], c = "r")

		max_axis = max(ax.get_xlim()[1], ax.get_ylim()[1])
		ax.set_xlim([0, max_axis])
		ax.set_ylim([0, max_axis])
		ax.set_xlabel("Simulation: ratio of tRNA to 16S rRNA")
		ax.set_ylabel("Expected: ratio of tRNA to 16S rRNA")
		plt.subplots_adjust(bottom = 0.2, top = 0.8)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
