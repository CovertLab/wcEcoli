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

COLORS = ['b','g','r','c','m','y','b']

def main(seedOutDir, plotOutDir, plotOutFileName, kbFile, metadata=None):

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

	nAxes = len(massNames) + 5
	fig, axesList = plt.subplots(nAxes, sharex = True)
	fig.set_figwidth(10)
	fig.set_figheight(10)

	for sim_idx, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = mass.readColumn("cellMass")

		# Plot masses
		for idx, massType in enumerate(massNames):
			massToPlot = mass.readColumn(massNames[idx])

			axesList[idx].plot(time / 60. / 60., massToPlot, linewidth = 2, color=COLORS[sim_idx])
			axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")


		# Plot doubling time of each cell
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		doublingTime = np.ones(time.shape) * (time.max() - initialTime) / 60.
		axesList[-5].plot(time / 60. / 60., doublingTime, linewidth = 2, color=COLORS[sim_idx])
		axesList[-5].set_ylabel("Doubling\ntime\n(min)")

		# Plot ppGpp
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		id_idx = bulkMolecules.readAttribute("objectNames").index("GUANOSINE-5DP-3DP[c]")
		ppGppCounts = bulkMolecules.readColumn("counts")[:, id_idx]

		cellMass = units.fg * cellMass
		cellVolume = cellMass / cellDensity
		ppGppConc = ((1 / nAvogadro) * (1 / cellVolume) * ppGppCounts).asNumber(units.umol / units.L)

		axesList[-4].plot(time / 60. / 60., ppGppConc, linewidth = 2, color=COLORS[sim_idx])
		axesList[-4].set_ylabel("[ppGpp]\n(uM)")

		# Plot ratio of initiated transcripts
		initiatedTranscripts = TableReader(os.path.join(simOutDir, "InitiatedTranscripts"))
		ratioStableToToalInitalized = initiatedTranscripts.readColumn("ratioStableToToalInitalized")
		N = 20
		runningMeanRatio = np.zeros(ratioStableToToalInitalized.size)
		for i in range(ratioStableToToalInitalized.size):
			runningMeanRatio[i] = np.mean(ratioStableToToalInitalized[i:i+N])

		axesList[-3].plot(time / 60. / 60., ratioStableToToalInitalized, color=COLORS[sim_idx])
		axesList[-3].plot(time / 60. / 60., runningMeanRatio, linestyle = '--', color = 'k')
		axesList[-3].set_ylabel('$r_s / r_t$')

		# Plot ribosome capacity
		ribosomeSubunitIds = []
		ribosomeSubunitIds.extend(kb.moleculeGroups.s50_fullComplex)
		ribosomeSubunitIds.extend(kb.moleculeGroups.s30_fullComplex)
		ribosomeSubunitIds.extend(kb.moleculeGroups.s50_proteinComplexes)
		ribosomeSubunitIds.extend(kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitIds'])
		ribosomeSubunitIds.extend(kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitIds'])

		moleculeIds = bulkMolecules.readAttribute("objectNames")
		ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)
		ribosomeSubunitCounts = bulkMolecules.readColumn("counts")[:, ribosomeSubunitIndexes]
		bulkMolecules.close()

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
		uniqueMoleculeCounts.close()

		massFile = TableReader(os.path.join(simOutDir, "RibosomeData"))
		actualElongations = massFile.readColumn("actualElongations")
		massFile.close()

		elongationRate = float(kb.constants.ribosomeElongationRate.asNumber(units.aa / units.s))
		totalRibosome = (activeRibosome + np.min(ribosomeSubunitCounts))
		totalRibosomeCapacity = totalRibosome * elongationRate

		axesList[-2].plot(time / 60. / 60., totalRibosomeCapacity, label="Theoretical total ribosome capacity", color=COLORS[sim_idx])
		axesList[-2].plot(time / 60. / 60., actualElongations, label="Actual elongations", color='k')
		axesList[-2].set_ylabel("AA\npoly-\nmerized")
		#axesList[-1].legend(ncol=2)

		# Plot glucose exchange flux
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		exFlux = fba_results.readColumn("externalExchangeFluxes")
		exMolec = fba_results.readAttribute("externalMoleculeIDs")
		glcFlux = exFlux[:,exMolec.index("GLC[p]")]

		axesList[-1].plot(time / 60. / 60., -1. * glcFlux, label="Glucose exchange flux coefficient", color=COLORS[sim_idx])
		axesList[-1].set_ylabel("External\nglucose\n(mmol/L/s)")

	for axes in axesList:
		possibleTicks = np.hstack((np.array(axes.get_ylim()), 0))
		possibleTicks = possibleTicks[possibleTicks.argsort()]
		if possibleTicks[0] == 0.:
			np.delete(possibleTicks, 0)
		axes.set_yticks(possibleTicks.tolist())

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
