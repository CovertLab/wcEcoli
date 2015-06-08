#!/usr/bin/env python

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(seedOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir)
	kb = cPickle.load(open(kbFile, "rb"))
	cellDensity = kb.constants.cellDensity
	nAvogadro = kb.constants.nAvogadro

	# Get all cells
	allDir = ap.getAll()

	massNames = [
				"dryMass",
				#"proteinMass",
				#"tRnaMass",
				"rRnaMass",
				'mRnaMass',
				#"dnaMass"
				]

	cleanNames = [
				"Dry\nmass",
				#"Protein\nmass",
				#"tRNA\nmass",
				"rRNA\nmass",
				"mRNA\nmass",
				#"DNA\nmass"
				]

	#plt.figure(figsize = (8.5, 11))
	nAxes = len(massNames) + 3
	fig, axesList = plt.subplots(nAxes, sharex = True)

	for simDir in allDir:
		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = mass.readColumn("cellMass")

		# Plot masses
		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx])

			axesList[idx].plot(time / 60. / 60., massToPlot, linewidth = 2)
			axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")


		# Plot doubling time of each cell
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		doublingTime = np.ones(time.shape) * (time.max() - initialTime) / 60.
		axesList[-3].plot(time / 60. / 60., doublingTime, linewidth = 2)
		axesList[-3].set_ylabel("Doubling\ntime\n(min)")

		# Plot ppGpp
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		id_idx = bulkMolecules.readAttribute("objectNames").index("PPGPP[c]")
		ppGppCounts = bulkMolecules.readColumn("counts")[:, id_idx]

		cellMass = units.fg * cellMass
		cellVolume = cellMass / cellDensity
		ppGppConc = ((1 / nAvogadro) * (1 / cellVolume) * ppGppCounts).asNumber(units.umol / units.L)

		axesList[-2].plot(time / 60. / 60., ppGppConc, linewidth = 2)
		axesList[-2].set_ylabel("[ppGpp]\n(uM)")

		# Plot ratio of initiated transcripts
		initiatedTranscripts = TableReader(os.path.join(simOutDir, "InitiatedTranscripts"))
		ratioStableToToalInitalized = initiatedTranscripts.readColumn("ratioStableToToalInitalized")
		N = 20
		runningMeanRatio = np.zeros(ratioStableToToalInitalized.size)
		for i in range(ratioStableToToalInitalized.size):
			runningMeanRatio[i] = np.mean(ratioStableToToalInitalized[i:i+N])

		axesList[-1].plot(time / 60. / 60., ratioStableToToalInitalized)
		axesList[-1].plot(time / 60. / 60., runningMeanRatio, linestyle = '--', color = 'k')
		axesList[-1].set_ylabel('$r_s / r_t$')

	for axes in axesList[:-1]:
		axes.get_ylim()
		axes.set_yticks(list(axes.get_ylim()))

	axesList[0].set_title("Cell mass fractions")
	axesList[-1].set_xlabel("Time (hr)")

	plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
