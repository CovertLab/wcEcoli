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
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]
	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIds])

	expression = sim_data.process.transcription.rnaExpression["basal"]
	mRnaExpression = np.array([expression[x] for x in mRnaIds])

	degRate = sim_data.process.transcription.rnaData["degRate"]
	mRnaDegRate = np.array([degRate[x].asNumber(1 / units.s) for x in mRnaIds])

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells()
	numMRnas = mRnaNames.shape[0]
	numCells = all_cells.shape[0]
	numerical_zero = 0.1
	

	# Get number of mRNA transcripts produced
	fig = plt.figure(figsize = (14, 10))
	ax = plt.subplot(1, 1, 1)
	
	for n, simDir in enumerate(all_cells):
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

		mRnaProduced = mRnaCounts_skip0 - (mRnaCounts_skiplast - mRnaDegraded_skip0)
		mRnaProducedSumOverTime = mRnaProduced.sum(axis = 0)

		# Replace true zero with numerical zero
		zerosIndex = np.where(mRnaProducedSumOverTime == 0)[0]
		mRnaProducedSumOverTime[zerosIndex] = numerical_zero
		log10_mRnasProduceSumOverTime = np.log10(mRnaProducedSumOverTime)
		log10_mRnasProduceSumOverTime[zerosIndex] = np.log10(numerical_zero)

		ax.scatter(np.log10(mRnaSynthProb), log10_mRnasProduceSumOverTime, facecolors = "none", edgecolors = "b")


	# Plot
	ax.set_title("Correlation of synthesis probability and number of transcripts produced\nn = %s cells" % numCells, fontsize = 12)
	ax.set_xlabel("log_10(synthesis probability)")
	ax.set_ylabel("log_10(Transcripts produced)")
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)

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
