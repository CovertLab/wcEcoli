#!/usr/bin/env python
"""
Plot expression of trp operon (designed for "000002_add_aa" shift)

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2018
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

SHIFT_TIME = 11000
FONTSIZE = 8
TEXT_POS_X = 0.9
TEXT_POS_Y = 0.8

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
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

	# Define mRNAs
	trpGeneNames = ["trpL", "trpE", "trpD", "trpC", "trpB", "trpA"]
	trpRnaIds = ["EG11274", "EG11028", "EG11027", "EG11026", "EG11025", "EG11024"]
	trpMonomerIds = ["EG11274-MONOMER", "ANTHRANSYNCOMPI-MONOMER", "ANTHRANSYNCOMPII-MONOMER", "PRAI-IGPS", "TRYPSYN-BPROTEIN", "TRYPSYN-APROTEIN"]
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

	FIRST_READ = True
	rnaCounts = None
	monomerCounts = None
	complexCounts = None

	for simDir in allDirs:
		simOutDir = os.path.join(simDir, "simOut")

		# Read from RnaSynthProb listener
		rnaSynthProbListener = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		rnaSynthProb = rnaSynthProbListener.readColumn('rnaSynthProb')[:, trpRnaSynthProbIdxs]
		time = rnaSynthProbListener.readColumn('time')
		rnaSynthProbListener.close()

		# Read from BulkMolecules listener
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		rnaCounts_thisGen = bulkMolecules.readColumn("counts")[:, trpRnaCountsIdxs]
		monomerCounts_thisGen = bulkMolecules.readColumn("counts")[:, trpMonomerCountsIdxs]
		complexCounts_thisGen = bulkMolecules.readColumn("counts")[:, trpComplexCountsIdxs]
		bulkMolecules.close()

		# Store values
		if FIRST_READ:
			FIRST_READ = False
			rnaCounts = rnaCounts_thisGen
			monomerCounts = monomerCounts_thisGen
			complexCounts = complexCounts_thisGen

		else:
			rnaCounts = np.vstack([rnaCounts, rnaCounts_thisGen])
			monomerCounts = np.vstack([monomerCounts, monomerCounts_thisGen])
			complexCounts = np.vstack([complexCounts, complexCounts_thisGen])

		# Plot
		for i, ax in enumerate(axesList):
			if i%2 == 0:
				ax.plot(time, rnaSynthProb[:, i//2])
				ax.set_ylabel(trpGeneNames[i//2])
			else:
				ax.plot(time, rnaCounts_thisGen[:, i//2])

	for ax in axesList:
		ax.set_xlim([0, time[-1]])
		ax.set_yticks(ax.get_ylim())
		ax.set_xticklabels([])

	for i, ax in enumerate(axesList_preshift):
		if rnaCounts[:SHIFT_TIME, i].max() == 0:
			ax.hist(rnaCounts[:SHIFT_TIME, i])
		else:
			ax.hist(rnaCounts[:SHIFT_TIME, i], bins = range(rnaCounts[:SHIFT_TIME, i].max()))
		ax.set_yticks(ax.get_ylim())
		ax.text(TEXT_POS_X, TEXT_POS_Y, "%0.2f" % np.mean(rnaCounts[:SHIFT_TIME, i]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right")

	for i, ax in enumerate(axesList_postshift):
		if rnaCounts[SHIFT_TIME:, i].max() == 0:
			ax.hist(rnaCounts[SHIFT_TIME:, i])
		else:
			ax.hist(rnaCounts[SHIFT_TIME:, i], bins = range(rnaCounts[SHIFT_TIME:, i].max()))
		ax.set_yticks(ax.get_ylim())
		ax.text(TEXT_POS_X, TEXT_POS_Y, "%0.2f" % np.mean(rnaCounts[SHIFT_TIME:, i]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right")

	for i, ax in enumerate(axesList_preshift_protein):
		ax.hist(monomerCounts[:SHIFT_TIME, i])
		ax.text(TEXT_POS_X, TEXT_POS_Y, "%0.2f" % np.mean(monomerCounts[:SHIFT_TIME, i]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right")
		if i in [1, 2]: # trpE, trpD
			ax.hist(2 * complexCounts[:SHIFT_TIME, 0], alpha = 0.5, color = "red")
			ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2, "%0.2f" % np.mean(2 * complexCounts[:SHIFT_TIME, 0]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right", fontdict = {"color": "red"})
		if i == 4: # trpB
			trpB_complexed = 2 * np.sum((complexCounts[:SHIFT_TIME, 1], complexCounts[:SHIFT_TIME, 2]), axis = 0)
			ax.hist(trpB_complexed, alpha = 0.5, color = "red")
			ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2, "%0.2f" % np.mean(trpB_complexed), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right", fontdict = {"color": "red"})
		if i == 5: # trpA
			ax.hist(2 * complexCounts[:SHIFT_TIME, 2], alpha = 0.5, color = "red")
			ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2, "%0.2f" % np.mean(2 * complexCounts[:SHIFT_TIME, 2]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right", fontdict = {"color": "red"})
		ax.set_yticks(ax.get_ylim())


	for i, ax in enumerate(axesList_postshift_protein):
		ax.hist(monomerCounts[SHIFT_TIME:, i])
		ax.text(TEXT_POS_X, TEXT_POS_Y, "%0.2f" % np.mean(monomerCounts[SHIFT_TIME:, i]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right")
		if i in [1, 2]: # trpE, trpD
			ax.hist(2 * complexCounts[SHIFT_TIME:, 0], alpha = 0.5, color = "red")
			ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2, "%0.2f" % np.mean(2 * complexCounts[SHIFT_TIME:, 0]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right", fontdict = {"color": "red"})
		if i == 4: # trpB
			trpB_complexed = 2 * np.sum((complexCounts[SHIFT_TIME:, 1], complexCounts[SHIFT_TIME:, 2]), axis = 0)
			ax.hist(trpB_complexed, alpha = 0.5, color = "red")
			ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2, "%0.2f" % np.mean(trpB_complexed), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right", fontdict = {"color": "red"})
		if i == 5: # trpA
			ax.hist(2 * complexCounts[SHIFT_TIME:, 2], alpha = 0.5, color = "red")
			ax.text(TEXT_POS_X, TEXT_POS_Y - 0.2, "%0.2f" % np.mean(2 * complexCounts[SHIFT_TIME:, 2]), transform = ax.transAxes, fontsize = FONTSIZE, horizontalalignment = "right", fontdict = {"color": "red"})
		ax.set_yticks(ax.get_ylim())

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

	plt.subplots_adjust(hspace = 0.8, wspace = 0.2, left = 0.075, right = 0.975, bottom = 0.075, top = 0.95)
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
