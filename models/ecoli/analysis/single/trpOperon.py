"""
Plot trp operon

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/9/2018
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get sim_data
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Get mRNA synthesis probabilities
		trpGeneNames = ["trpL", "trpE", "trpD", "trpC", "trpB", "trpA"]
		trpRnaIds = ["EG11274", "EG11028", "EG11027", "EG11026", "EG11025", "EG11024"]
		rnaSynthProbListener = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		rnaIds = rnaSynthProbListener.readAttribute('rnaIds')
		trpRnaIdxs = [rnaIds.index("%s_RNA[c]" % x) for x in trpRnaIds]
		rnaSynthProb = rnaSynthProbListener.readColumn('rnaSynthProb')[:, trpRnaIdxs]
		time = rnaSynthProbListener.readColumn('time')
		rnaSynthProbListener.close()

		# Get mRNA counts
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		trpRnaIdxs = [moleculeIds.index("%s_RNA[c]" % x) for x in trpRnaIds]
		rnaCounts = bulkMolecules.readColumn("counts")[:, trpRnaIdxs]
		bulkMolecules.close()

		# Plot
		rows = 12
		cols = 1
		fig, axesList = plt.subplots(rows, cols, figsize = (8.5, 11), sharex = True)

		for iGene, iAx in enumerate(xrange(0, rows, 2)):
			axSynthProb = axesList[iAx]
			axSynthProb.plot(time, rnaSynthProb[:, iGene])
			axSynthProb.set_ylabel(trpGeneNames[iGene])
			ymin = rnaSynthProb[:, iGene].min()
			ymax = rnaSynthProb[:, iGene].max()
			if ymin == ymax:
				axSynthProb.set_yticks([ymin])
			else:
				axSynthProb.set_yticks([ymin, ymax])
			axSynthProb.text(-0.3, 0.5, "synth\nprob", transform = axSynthProb.transAxes, va = "center", ha = "right")
			fractionBound = (rnaSynthProb[:, iGene] == ymin).sum(dtype = float) / len(time)
			axSynthProb.text(1.025, 0.5, "%0.4f" % fractionBound, transform = axSynthProb.transAxes)

			if iGene == 0:
				axSynthProb.text(1.025, 1.5, "fraction\nbound", transform = axSynthProb.transAxes)

			axRnaCount = axesList[iAx + 1]
			axRnaCount.plot(time, rnaCounts[:, iGene])
			ymin = rnaCounts[:, iGene].min()
			ymax = rnaCounts[:, iGene].max()
			if ymin == ymax:
				axRnaCount.set_yticks([ymin])
			else:
				axRnaCount.set_yticks([ymin, ymax])
			axRnaCount.text(-0.3, 0.5, "mRNA\ncount", transform = axRnaCount.transAxes, va = "center", ha = "right")

		axesList[0].set_xlim([time[0], time[-1]])
		plt.subplots_adjust(hspace = 0.8, left = 0.3, right = 0.9)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
