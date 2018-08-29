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
		aaIds = sim_data.moleculeGroups.aaIDs

		# Remove selenocysteine
		selenocysteineId = "L-SELENOCYSTEINE[c]"
		selenocysteineIndex = aaIds.index(selenocysteineId)
		aaIds.remove(selenocysteineId)
		aasUsed = np.delete(aasUsed, selenocysteineIndex, 1)
		rnaCounts_summedByAa = np.delete(rnaCounts_summedByAa, selenocysteineIndex, 1)

		# Plot
		xlabels = []
		fig, ax = plt.subplots(1, 1, figsize = (11, 8.5))
		for i, aaId in enumerate(aaIds):
			ratio = aasUsed[:, i] / rnaCounts_summedByAa[:, i]
			std = np.std(ratio)
			print std
			ax.errorbar(i, np.average(ratio), yerr=std, marker='o', mfc='b', mec='b', ecolor='b', capsize = 2)

			if len(aaId) < 7:
				xlabels.append(aaId[:-3])
			elif aaId == 'L-ALPHA-ALANINE[c]':
				xlabels.append("ALA")
			elif aaId == 'L-ASPARTATE[c]':
				xlabels.append("ASP")
			else:
				print "unexpected amino acid during labeling"

		ax.set_xticks(range(len(aaIds)))
		ax.set_xticklabels(xlabels)
		ax.axhline(y = 1, c = 'b')
		ax.set_ylabel("amino acids used / number of tRNAs")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
