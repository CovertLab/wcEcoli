"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/17/2017
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils.sparkline import whitePadSparklineAxis


FRACTION_EFTU_AATRNA = 0.84
FACTOR_RIBOSOME_AATRNA = 5.5

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Load validation data
		validation_data = cPickle.load(open(validationDataFile, "rb"))
		kinetic_data = validation_data.trna.synthetaseKinetics
		synthetaseIds = kinetic_data.keys()

		# Load time data
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		# Load counts of EF-Tu, ribosome subunits
		EFTu1Id = "EG11036-MONOMER[i]"
		EFTu2Id = "EG11037-MONOMER[c]"
		ribosome30SIds = sim_data.moleculeGroups.s30_proteins
		ribosome50SIds = sim_data.moleculeGroups.s50_proteins

		bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMoleculesReader.readAttribute("objectNames")
		EFTu1Index = moleculeIds.index(EFTu1Id)
		EFTu2Index = moleculeIds.index(EFTu2Id)
		ribosome30SIndices = np.array([moleculeIds.index(x) for x in ribosome30SIds], np.int)
		ribosome50SIndices = np.array([moleculeIds.index(x) for x in ribosome50SIds], np.int)
		synthetaseIndices = np.array([moleculeIds.index(x) for x in synthetaseIds], np.int)
		moleculeCounts = bulkMoleculesReader.readColumn("counts")
		bulkMoleculesReader.close()
		EFTu1Counts = moleculeCounts[:, EFTu1Index]
		EFTu2Counts = moleculeCounts[:, EFTu2Index]
		ribosome30SCounts = moleculeCounts[:, ribosome30SIndices]
		ribosome50SCounts = moleculeCounts[:, ribosome50SIndices]
		synthetaseCounts = moleculeCounts[:, synthetaseIndices]

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
		fig, axesList = plt.subplots(4, 1, figsize = (8.5, 11))
		ax1, ax2, ax3, ax4 = axesList

		ax1.plot(time, EFTu1Counts, color = "b", linestyle = "--")
		ax1.plot(time, EFTu2Counts, color = "b", linestyle = "--")
		ax1.plot(time, EFTu1Counts + EFTu2Counts, color = "b")
		ax1.text(time[-1], EFTu1Counts[-1], "EF-Tu 1")
		ax1.text(time[-1], EFTu2Counts[-1], "EF-Tu 2")
		ax1.text(time[-1], EFTu1Counts[-1] + EFTu2Counts[-1], "EF-Tu (total)")

		ax2.plot(time, FRACTION_EFTU_AATRNA * (EFTu1Counts + EFTu2Counts), color = "b")
		ax2.plot(time, FACTOR_RIBOSOME_AATRNA * activeRibosomeCounts, color = "g")
		ax2.text(time[-1], FRACTION_EFTU_AATRNA * (EFTu1Counts + EFTu2Counts)[-1], "%s EF-Tu total" % FRACTION_EFTU_AATRNA)
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
		ax3.plot(time, ribosome30SCounts.sum(axis = 1), color = "g", linestyle = "--")
		ax3.plot(time, ribosome50SCounts.sum(axis = 1), color = "g", linestyle = "--")
		ax3.text(time[-1], ribosome30SCounts[-1, :].sum(), "30S subunit")
		ax3.text(time[-1], ribosome50SCounts[-1, :].sum(), "50S subunit")
		ax3.text(time[-1], activeRibosomeCounts[-1], "Ribosomes (active)")

		ax4.plot(time, np.array(synthetaseCounts.sum(axis = 1), np.float) / (EFTu1Counts + EFTu2Counts))
		ax4.set_title("total synthetase / total EF-Tu")

		for i, ax in enumerate(axesList):
			if i < 3:
				ax.set_ylabel("counts")
			whitePadSparklineAxis(ax)
		plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.15, bottom = 0.1, top = 0.9, right = 0.7)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
