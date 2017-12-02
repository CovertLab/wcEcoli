#!/usr/bin/env python
"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2017
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

FRACTION_EFTU_AATRNA = 0.84
FACTOR_RIBOSOME_AATRNA = 5.5

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	return
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"
	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load sim data
	sim_data = cPickle.load(open(simDataFile, "rb"))

	# Load validation data
	# Load kinetic data from validation
	validation_data = cPickle.load(open(validationDataFile, "rb"))
	kinetic_data = validation_data.trna.synthetaseKinetics
	synthetaseIds = kinetic_data.keys()

	# Load time data
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

	# Load counts of EF-Tu, ribosome subunits
	EFTu1Id = "EG11036-MONOMER[i]"
	EFTu2Id = "EG11037-MONOMER[c]"
	ribosome30SId = sim_data.moleculeGroups.s30_fullComplex[0]
	ribosome50SId = sim_data.moleculeGroups.s50_fullComplex[0]

	bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMoleculesReader.readAttribute("objectNames")
	moleculeIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in [EFTu1Id, EFTu2Id, ribosome30SId, ribosome50SId]], np.int)
	synthetaseIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in synthetaseIds])
	moleculeCounts = bulkMoleculesReader.readColumn("counts")[:, moleculeIndexes]
	synthetaseCounts = bulkMoleculesReader.readColumn("counts")[:, synthetaseIndexes]
	bulkMoleculesReader.close()

	# Load counts of active ribosomes
	uniqueMoleculeCountsReader = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
	ribosomeIndex = uniqueMoleculeCountsReader.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	activeRibosomeCounts = uniqueMoleculeCountsReader.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
	uniqueMoleculeCountsReader.close()

	# Load counts of amino acids incorporated into polypeptides
	aaIds = sim_data.moleculeGroups.aaIDs
	selenocysteine_index = aaIds.index("L-SELENOCYSTEINE[c]")
	aasUsed = TableReader(os.path.join(simOutDir, "GrowthLimits")).readColumn("aasUsed")
	aasUsed = np.delete(aasUsed, selenocysteine_index, 1)

	# Plot
	fig, axesList = plt.subplots(4, 1, figsize = (10, 10))
	ax1, ax2, ax3, ax4 = axesList

	ax1.plot(time, moleculeCounts[:, 0], color = "b", linestyle = "--")
	ax1.plot(time, moleculeCounts[:, 1], color = "b", linestyle = "--")
	ax1.plot(time, np.sum(moleculeCounts[:, :2], axis = 1), color = "b")
	ax1.text(time[-1], moleculeCounts[-1, 0], "EF-Tu 1")
	ax1.text(time[-1], moleculeCounts[-1, 1], "EF-Tu 2")
	ax1.text(time[-1], np.sum(moleculeCounts[-1, :2]), "EF-Tu (total)")

	ax2.plot(time, FRACTION_EFTU_AATRNA * np.sum(moleculeCounts[:, :2], axis = 1), color = "b")
	ax2.plot(time, FACTOR_RIBOSOME_AATRNA * activeRibosomeCounts, color = "g")
	ax2.text(time[-1], FRACTION_EFTU_AATRNA * np.sum(moleculeCounts[-1, :2]), "%s EF-Tu total" % FRACTION_EFTU_AATRNA)
	ax2.text(time[-1], FACTOR_RIBOSOME_AATRNA * activeRibosomeCounts[-1], "%s Ribosomes (active)" % FACTOR_RIBOSOME_AATRNA)

	aminoacylTrnaCounts = cPickle.load(open(os.path.join(plotOutDir, "aminoacylTrnaCounts.cPickle"), "rb"))
	totalAminoacylTrnaCounts = []
	for aa in aminoacylTrnaCounts.keys():
		totalAminoacylTrnaCounts.append(aminoacylTrnaCounts[aa])
	totalAminoacylTrnaCounts = np.array(totalAminoacylTrnaCounts)
	ax2.plot(time, np.sum(totalAminoacylTrnaCounts, axis = 0), color = "r")
	ax2.text(time[-1], totalAminoacylTrnaCounts[:, -1].sum(), "aminoacylated-tRNA (kinetic estimation)")

	ax2.plot(time, aasUsed.sum(axis = 1), "r")
	ax2.text(time[-1], aasUsed[-1, :].sum(), "aminoacylated-tRNA (demand estimation)")

	ax3.plot(time, activeRibosomeCounts, color = "g")
	ax3.plot(time, moleculeCounts[:, 2], color = "g", linestyle = "--")
	ax3.plot(time, moleculeCounts[:, 3], color = "g", linestyle = "--")
	ax3.text(time[-1], moleculeCounts[-1, 2], "30S subunit")
	ax3.text(time[-1], moleculeCounts[-1, 3], "50S subunit")
	ax3.text(time[-1], activeRibosomeCounts[-1], "Ribosomes (active)")

	ax4.plot(time, np.array(synthetaseCounts.sum(axis = 1), np.float) / np.sum(moleculeCounts[:, :2], axis = 1))
	ax4.set_title("total synthetase / total EF-Tu")

	for ax in axesList:
		whitePadSparklineAxis(ax)
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.9, right = 0.85)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

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