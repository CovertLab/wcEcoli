#!/usr/bin/env python
"""
Plots fraction of mRNAs transcribed (out of all genes to be transcribed) for all seeds.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/8/2016
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get IDs of mRNAs
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	mRnaIds = np.where(isMRna)[0]

	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]
	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIds])
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# basalExpression = sim_data.process.transcription.rnaExpression["basal"]
	# degRate = sim_data.process.transcription.rnaData["degRate"]

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells()
	numMRnas = mRnaNames.shape[0]
	numCells = all_cells.shape[0]
	

	# Get number of mRNA transcripts produced
	transcribedBoolean = []
	for n, simDir in enumerate(all_cells):
		print n

		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNames])
		mRnaCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]
		bulkMolecules.close()

		rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
	 	countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
	 	countMRnaDegraded = countRnaDegraded[:, mRnaIds]

		
		# Calculate number of transcripts produced
		mRnaDegraded_skip0 = countMRnaDegraded[1:, :] # remove first timestep
		mRnaCounts_skip0 = mRnaCounts[1:, :]
		mRnaCounts_skiplast = mRnaCounts[:-1, :]

		mRnaProducedCounts = mRnaCounts_skip0 - (mRnaCounts_skiplast - mRnaDegraded_skip0)
		mRnaProducedCountsSumOverTime = mRnaProducedCounts.sum(axis = 0)


		# Record T/F if at least 1 transcript was produced
		mRnasTranscribed = np.array([x != 0 for x in mRnaProducedCountsSumOverTime])
		transcribedBoolean.append(mRnasTranscribed)

	transcribedBoolean = np.array(transcribedBoolean)
	transcribedBooleanSumOverCells = transcribedBoolean.sum(axis = 0)
	transcribedFreq = transcribedBooleanSumOverCells / float(numCells)

	# Plot
	fig = plt.figure(figsize = (14, 10))
	ax = plt.subplot(1, 1, 1)

	ax.scatter(np.log10(mRnaSynthProb), transcribedFreq, facecolors = "none", edgecolors = "b")
	ax.set_title("Correlation of synthesis probability and frequency of producing at least 1 transcript\nn = %s cells" % numCells, fontsize = 12)
	ax.set_xlabel("log_10(synthesis probability)")
	ax.set_ylabel("Frequency of producing at least 1 transcript")
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)
	ax.set_ylim([-0.1, 1.1])

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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
