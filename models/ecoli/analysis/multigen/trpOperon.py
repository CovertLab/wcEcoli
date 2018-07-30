"""
Plot expression of trp operon (designed for "000002_add_aa" shift)

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2018

todo:
- utilize SHIFT_TIME as time (seconds) rather than an index
- move this plot to a variant analysis script, so that we can check that this
	analysis only occurs for the upshift.
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
from models.ecoli.analysis import multigenAnalysisPlot


SHIFT_TIME = 11000
FONTSIZE = 8
TEXT_POS_X = 0.9
TEXT_POS_Y = 0.8

def plotComplexes(ax, plotComplexCounts):
	ax.hist(plotComplexCounts, alpha = 0.5, color = "red")
	ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2,
		"%0.2f" % np.mean(plotComplexCounts),
		transform = ax.transAxes,
		fontsize = FONTSIZE,
		ha = "right",
		fontdict = {"color": "red"},
		)
	return


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load directories
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDirs = ap.get_cells()

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity

		# Define molecule IDs
		trpGeneNames = ["trpL", "trpE", "trpD", "trpC", "trpB", "trpA"]
		trpRnaIds = ["EG11274", "EG11028", "EG11027", "EG11026", "EG11025", "EG11024"]
		trpMonomerIds = [
			"EG11274-MONOMER",
			"ANTHRANSYNCOMPI-MONOMER",
			"ANTHRANSYNCOMPII-MONOMER",
			"PRAI-IGPS",
			"TRYPSYN-BPROTEIN",
			"TRYPSYN-APROTEIN",
			]
		trpComplexIds = [
			"ANTHRANSYN-CPLX", # 2 trpE + 2 trpD
			"CPLX0-2401", # 2 trpB
			"TRYPSYN", # 2 trpA + 1 CPLX0-2401
			]

		# Select first simOutDir to get indices
		simOutDir = os.path.join(allDirs[0], "simOut")

		# Get mRNA synthesis probabilities
		rnaSynthProbListener = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		rnaIds = rnaSynthProbListener.readAttribute('rnaIds')
		trpRnaSynthProbIdxs = [rnaIds.index("%s_RNA[c]" % x) for x in trpRnaIds]
		time = rnaSynthProbListener.readColumn('time')
		rnaSynthProbListener.close()

		# Get mRNA counts
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		trpRnaCountsIdxs = [moleculeIds.index("%s_RNA[c]" % x) for x in trpRnaIds]
		trpMonomerCountsIdxs = [moleculeIds.index("%s[c]" % x) for x in trpMonomerIds]
		trpComplexCountsIdxs = [moleculeIds.index("%s[c]" % x) for x in trpComplexIds]
		bulkMolecules.close()

		# Plot
		nRows = 12
		nCols = 5
		fig = plt.figure(figsize = (11, 8.5))
		axesList = [plt.subplot2grid((nRows, nCols), (x, 0)) for x in xrange(nRows)]
		axesList_preshift = [plt.subplot2grid((nRows, nCols), (x, 1), rowspan = 2) for x in xrange(0, nRows, 2)]
		axesList_postshift = [plt.subplot2grid((nRows, nCols), (x, 2), rowspan = 2, sharex = axesList_preshift[i], sharey = axesList_preshift[i]) for i, x in enumerate(xrange(0, nRows, 2))]
		axesList_preshift_protein = [plt.subplot2grid((nRows, nCols), (x, 3), rowspan = 2) for x in xrange(0, nRows, 2)]
		axesList_postshift_protein = [plt.subplot2grid((nRows, nCols), (x, 4), rowspan = 2, sharex = axesList_preshift_protein[i], sharey = axesList_preshift_protein[i]) for i, x in enumerate(xrange(0, nRows, 2))]

		rnaCounts = []
		monomerCounts = []
		complexCounts = []
		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")

			# Read from RnaSynthProb listener
			rnaSynthProbListener = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			rnaSynthProb = rnaSynthProbListener.readColumn('rnaSynthProb')[:, trpRnaSynthProbIdxs]
			time = rnaSynthProbListener.readColumn('time')
			rnaSynthProbListener.close()

			# Read from BulkMolecules listener
			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")
			bulkMolecules.close()
			rnaCounts_thisGen = moleculeCounts[:, trpRnaCountsIdxs]
			monomerCounts_thisGen = moleculeCounts[:, trpMonomerCountsIdxs]
			complexCounts_thisGen = moleculeCounts[:, trpComplexCountsIdxs]

			# Store values
			rnaCounts.append(rnaCounts_thisGen)
			monomerCounts.append(monomerCounts_thisGen)
			complexCounts.append(complexCounts_thisGen)

			# Plot rna synthesis prob and rna counts
			for i, ax in enumerate(axesList):
				if i%2 == 0:
					ax.plot(time, rnaSynthProb[:, i//2])
					ax.set_ylabel(trpGeneNames[i//2])
				else:
					ax.plot(time, rnaCounts_thisGen[:, i//2])

		# Convert to numpy arrays
		rnaCounts = np.array(rnaCounts)[0]
		monomerCounts = np.array(monomerCounts)[0]
		complexCounts = np.array(complexCounts)[0]

		# Plot preshift and postshift rna counts
		startList = [0, SHIFT_TIME]
		endList = [SHIFT_TIME, None]

		for j, aList in enumerate([axesList_preshift, axesList_postshift]):
			start = startList[j]
			end = endList[j]

			# Plot rna counts
			for i, ax in enumerate(aList):
				plotCounts = rnaCounts[start:end, i]
				ax.hist(plotCounts)
				ax.text(TEXT_POS_X, TEXT_POS_Y,
					"%0.2f" % np.mean(plotCounts),
					transform = ax.transAxes,
					fontsize = FONTSIZE,
					ha = "right",
					)

		# Plot preshift and postshit protein counts
		for j, aList in enumerate([axesList_preshift_protein, axesList_postshift_protein]):
			start = startList[j]
			end = endList[j]

			# Plot protein counts
			for i, ax in enumerate(aList):
				plotCounts = monomerCounts[start:end, i]
				ax.hist(plotCounts)
				ax.text(TEXT_POS_X, TEXT_POS_Y,
					"%0.2f" % np.mean(plotCounts),
					transform = ax.transAxes,
					fontsize = FONTSIZE,
					ha = "right",
					)

				# Plot complexes
				if i in [1, 2]: # trpE, trpD
					plotComplexCounts = 2 * complexCounts[start:end, 0]
					plotComplexes(ax, plotComplexCounts)
				if i == 4: # trpB
					plotComplexCounts = 2 * np.sum((complexCounts[start:end, 1], complexCounts[start:end, 2]), axis = 0)
					plotComplexes(ax, plotComplexCounts)
				if i == 5: # trpA
					plotComplexCounts = 2 * complexCounts[start:end, 2]
					plotComplexes(ax, plotComplexCounts)

		for ax in [axesList_preshift[0], axesList_preshift_protein[0]]:
			ax.set_title("Pre-shift")

		for ax in [axesList_postshift[0], axesList_postshift_protein[0]]:
			ax.set_title("Post-shift")

		for ax in [axesList_preshift[-1], axesList_postshift[-1]]:
			ax.set_xlabel("mRNA count")

		for ax in [axesList_preshift_protein[-1], axesList_postshift_protein[-1]]:
			ax.set_xlabel("protein count")

		for aList in [axesList, axesList_preshift, axesList_postshift, axesList_preshift_protein, axesList_postshift_protein]:
			for ax in aList:
				ax.tick_params(axis = "both", labelsize = FONTSIZE)
				if len(ax.get_xticks()) > 6:
					xticks = ax.get_xticks()[::2]
					ax.set_xticks(xticks)

		for ax in axesList:
			ax.set_xlim([0, time[-1]])
			ax.set_yticks(ax.get_ylim())
			ax.set_xticklabels([])
		plt.subplots_adjust(hspace = 0.8, wspace = 0.2, left = 0.1, right = 0.975, bottom = 0.075, top = 0.95)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
