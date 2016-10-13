#!/usr/bin/env python
"""
Plots fraction of mRNAs transcribed (out of all genes to be transcribed) for all seeds.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/29/2016
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

	# import ipdb; ipdb.set_trace()

	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get IDs of mRNAs
	sim_data = cPickle.load(open(simDataFile, "rb"))

	# rnaIds = sim_data.process.transcription.rnaData["id"]
	# isMRna = sim_data.process.transcription.rnaData["isMRna"]
	# mRnaIds = np.where(isMRna)[0]

	mRnaNames = sim_data.process.transcription.rnaData["id"][sim_data.relation.rnaIndexToMonomerMapping]
	proteinNames = sim_data.process.translation.monomerData["id"]

	mRnaBasalExpression = sim_data.process.transcription.rnaExpression["basal"][sim_data.relation.rnaIndexToMonomerMapping]
	mRnaSynthProb = sim_data.process.transcription.rnaSynthProb["basal"][sim_data.relation.rnaIndexToMonomerMapping]
	mRnaDegRate = sim_data.process.transcription.rnaData["degRate"][sim_data.relation.rnaIndexToMonomerMapping]
	

	# mRnaBasalExpression = np.array([basalExpression[x] for x in mRnaIds])
	# mRnaSynthProb = np.array([synthProb[x] for x in mRnaIds])
	# mRnaDegRate = np.array([degRate[x].asNumber(1 / units.s) for x in mRnaIds])
	# mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	# Sort in order
	# descendingOrderIndexing = np.argsort(mRnaBasalExpression)[::-1]
	# mRnaNamesSorted = mRnaNames[descendingOrderIndexing]
	# mRnaBasalExpressionSorted = mRnaBasalExpression[descendingOrderIndexing]
	# mRnaSynthProbSorted = mRnaSynthProb[descendingOrderIndexing]
	# mRnaDegRateSorted = mRnaDegRate[descendingOrderIndexing]

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells()

	# Get number of mRNAs transcribed
	transcribedFreq = []
	translatedFreq = []
	PresentFreqAvg = []
	ProteinAvg = []

	# import ipdb; ipdb.set_trace()

	# post-processing NPZ file
	# outfile = os.path.join(plotOutDir) + 'FrequenciesRNAandProtein.npz'
	# npzfile = np.load(outfile)
	# DATA =  np.column_stack((npzfile['ProteinId'], npzfile['mRnaId'], npzfile['degRate'], npzfile['synthProb'], npzfile['ProteinFreq'], npzfile['mRnaFreq'], npzfile['PresentProteinFreq'], npzfile['ProteinAvg'], npzfile['expression']))
	# np.savetxt(os.path.join(plotOutDir) + 'FrequenciesRNAandProtein.txt', DATA, delimiter = "\t", fmt = "%s")

	for n, simDir in enumerate(all_cells):
		print n
		
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")

		mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNames])
		mRnaMoleculeCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]

		ProteinIndexes = np.array([moleculeIds.index(x) for x in proteinNames])
		ProteinMoleculeCounts = bulkMolecules.readColumn("counts")[:, ProteinIndexes]

		bulkMolecules.close()

		# rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		# countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
		# countMRnaDegraded = countRnaDegraded[:, mRnaIds]

		mRnaMoleculeCountsSumOverTime = mRnaMoleculeCounts.sum(axis = 0)
		mRnasTranscribed = np.array([x != 0 for x in mRnaMoleculeCountsSumOverTime])

		ProteinMoleculeCountsSumOverTime = ProteinMoleculeCounts.sum(axis = 0)
		ProteinTranslated = np.array([x != 0 for x in ProteinMoleculeCountsSumOverTime])


		PresentFreq = []
		ProteinAvgOverSim = []

		for idx in xrange(0, len(ProteinMoleculeCountsSumOverTime)):
			# print ProteinMoleculeCounts[:,idx]
			ProteinPresentOverTime = np.array([x != 0 for x in ProteinMoleculeCounts[:,idx]])

			ProteinMoleculeCountsGene = ProteinMoleculeCounts[:,idx]
			
			ProteinAvgGene = np.mean(ProteinMoleculeCountsGene[ProteinMoleculeCountsGene > 0])
			if len(ProteinMoleculeCountsGene[ProteinMoleculeCountsGene > 0]) == 0:
				ProteinAvgGene = 0
				
			ProteinAvgOverSim.append(ProteinAvgGene)

			# import ipdb; ipdb.set_trace()

			# print sum(ProteinPresentOverTime) / float(len(ProteinPresentOverTime))
			PresentFreq.append(sum(ProteinPresentOverTime) / float(len(ProteinPresentOverTime)))

		PresentFreqAvg = PresentFreqAvg + PresentFreq
		ProteinAvg = ProteinAvg + ProteinAvgOverSim

		# import ipdb; ipdb.set_trace()

		transcribedFreq.append(mRnasTranscribed)
		translatedFreq.append(ProteinTranslated)


	transcribedFreq = np.array(transcribedFreq)
	transcribedFreqSumOverSeeds = transcribedFreq.sum(axis = 0)

	translatedFreq = np.array(translatedFreq)
	translatedFreqSumOverSeeds = translatedFreq.sum(axis = 0)


	# Save data
	numMRnas = mRnaNames.shape[0]
	numCells = all_cells.shape[0]

	mRnaFreq = transcribedFreqSumOverSeeds / float(numCells)
	ProteinFreq = translatedFreqSumOverSeeds / float(numCells)
	ProteinFreq = translatedFreqSumOverSeeds / float(numCells)

	ProteinFreq = translatedFreqSumOverSeeds / float(numCells)

	PresentFreqAvg = np.array(PresentFreqAvg)
	PresentFreqAvg = PresentFreqAvg / float(numCells)

	ProteinAvg = np.array(ProteinAvg)
	ProteinAvg = ProteinAvg / float(numCells)


	# import ipdb; ipdb.set_trace()
	np.savez(os.path.join(plotOutDir) + 'FrequenciesRNAandProtein.npz', 
		mRnaFreq = mRnaFreq, 
		ProteinFreq = ProteinFreq,
		PresentProteinFreq = PresentFreqAvg,
		mRnaId = mRnaNames, 
		ProteinId = proteinNames,
		expression = mRnaBasalExpression, 
		synthProb = mRnaSynthProb, 
		degRate = mRnaDegRate,
		numCells = numCells,
		ProteinAvg = ProteinAvg,
	)

	# Plot
	fig = plt.figure(figsize = (14, 10))

	ax = plt.subplot(3, 1, 1)

	# ax.scatter(np.arange(numMRnas), transcribedFreqSumOverSeeds / float(numCells), facecolors = "none", edgecolors = "b")
	ax.scatter(np.log10(mRnaSynthProb), transcribedFreqSumOverSeeds / float(numCells), facecolors = "none", edgecolors = "b")
	ax.set_title("Correlation of synthesis probability and frequency of observing at least 1 transcript\nn = %s cells" % numCells, fontsize = 12)
	ax.set_xlabel("Log10(synthesis probability)")
	ax.set_ylabel("Prob(>=1 transcript)")
	ax.set_ylim([-.05, 1.05])
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)

	ax = plt.subplot(3, 1, 2)

	ax.scatter(np.log10(mRnaSynthProb), translatedFreqSumOverSeeds / float(numCells), facecolors = "none", edgecolors = "b")
	ax.set_title("Correlation of synthesis probability and frequency of observing at least 1 protein\nn = %s cells" % numCells, fontsize = 12)
	ax.set_xlabel("Log10(synthesis probability)")
	ax.set_ylabel("Prob(>=1 transcript)")
	ax.set_ylim([-.05, 1.05])
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)

	ax = plt.subplot(3, 1, 3)

	ax.scatter(transcribedFreqSumOverSeeds / float(numCells), translatedFreqSumOverSeeds / float(numCells), facecolors = "none", edgecolors = "b")
	ax.set_title("Correlation between frequency of observing at least 1 transcript and 1 protein\nn = %s cells" % numCells, fontsize = 12)
	ax.set_xlabel("Prob(>=1 transcript)")
	ax.set_ylabel("Prob(>=1 protein)")
	ax.set_xlim([-.05, 1.05])
	ax.set_ylim([-.05, 1.05])
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)

	plt.subplots_adjust(hspace = 0.9, wspace = 0.5)

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
