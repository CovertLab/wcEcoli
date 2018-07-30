"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/3/2017
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units


SKIP_FIRST_N_TIMESTEPS = 10
LABEL = False

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load kinetic data from validation
		validation_data = cPickle.load(open(validationDataFile, "rb"))
		kinetic_data = validation_data.trna.synthetaseKinetics
		synthetaseIds = kinetic_data.keys()

		# Load from sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		aaIds = sim_data.moleculeGroups.aaIDs

		# Load data from sim
		timestep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
		cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass") * units.fg
		cellVolume = cellMass / cellDensity

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		synthetaseIndexes = np.array([moleculeIds.index(x) for x in synthetaseIds], np.int)
		synthetaseCounts = bulkMolecules.readColumn("counts")[:, synthetaseIndexes]
		bulkMolecules.close()

		# Compute synthetase concentrations
		synthetaseConcentrations = (synthetaseCounts / nAvogadro) / cellVolume[:, None]

		# Load number of incorporations of each amino acid from sim
		aasUsed = TableReader(os.path.join(simOutDir, "GrowthLimits")).readColumn("aasUsed")
		selenocysteine_index = aaIds.index("L-SELENOCYSTEINE[c]")
		aaIds.remove("L-SELENOCYSTEINE[c]")
		aasUsed = np.delete(aasUsed, selenocysteine_index, 1)

		# Plot
		fig, axesList = plt.subplots(1, 3, figsize = (11, 8.5), sharey = True)
		substrateList = ["tRNA", "AA", "ATP"]
		ax1, ax2, ax3 = axesList
		for i, synthetaseId in enumerate(synthetaseIds):
			aaId = kinetic_data[synthetaseId]["amino acid"]
			aaIndex = aaIds.index(aaId)
			E = synthetaseConcentrations[:, i]
			for j, substrate in enumerate(substrateList):
				kcat = kinetic_data[synthetaseId]["kcat %s" % substrate]
				if kcat.asNumber() == 0:
					continue

				vmax = E*kcat
				charged_tRNAs_max = np.floor(np.array(vmax * timestep*units.s * cellVolume * nAvogadro, dtype = float))
				y = np.average(charged_tRNAs_max[SKIP_FIRST_N_TIMESTEPS:] / aasUsed[SKIP_FIRST_N_TIMESTEPS:, aaIndex])
				x1 = kcat.asNumber(1/units.s)
				x2 = np.average([x.asNumber() for x in E]) * 1e6
				x3 = (np.average(aasUsed[:, aaIndex]) * 1e-3)

				for ax, x in zip(axesList, [x1, x2, x3]):
					ax.scatter(x, y, edgecolor = "none")

					if LABEL:
						ax.text(x, y, aaId, fontsize = 4)

		ax1.set_xlabel("kcat (1/s)")
		ax2.set_xlabel("[synthetase] (uM)")
		ax3.set_xlabel("demand (1000s amino acids)")
		ax1.set_ylabel("Ratio of (supply / demand)")
		for ax in axesList:
			ax.tick_params(right = "off", top = "off", which = "both", direction = "out")
		plt.subplots_adjust(wspace = 0.5, bottom = 0.3, top = 0.7)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
