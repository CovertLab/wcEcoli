"""
Plot polypeptide elongation rate vs tRNA synthetase kinetic capacity

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/30/2018
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load sim data
		sim_data = cPickle.load(open(simDataFile))
		aaIds = sim_data.moleculeGroups.aaIDs
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass") * units.fg
		cellVolume = cellMass / cellDensity

		# Get timesteps
		timestep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

		# Get number of amino acids used per timestep
		growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))
		aasUsed = growthLimits.readColumn("aasUsed")
		growthLimits.close()

		# Convert amino acid counts to concentrations
		aasUsedConcentrations = aasUsed / nAvogadro / cellVolume[:, None]

		# Compute per-amino acid elongation rates
		elongationRates = aasUsedConcentrations / (timestep[:, None] * units.s)

		# Get tRNA synthetase kcats
		validation_data = cPickle.load(open(validationDataFile, "rb"))
		kinetic_data = validation_data.trna.synthetaseKinetics

		# Get counts of tRNA synthetases
		synthetaseIds = kinetic_data.keys()
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		synthetaseIndices = np.array([moleculeIds.index(x) for x in synthetaseIds], np.int)
		synthetaseCounts = bulkMolecules.readColumn("counts")[:, synthetaseIndices]
		bulkMolecules.close()

		# Convert tRNA synthetase counts to concentrations
		synthetaseConcentrations = synthetaseCounts / nAvogadro / cellVolume[:, None]

		# Plot
		fig, ax = plt.subplots(1, 1, figsize = (8.5, 11))
		ax.set_xlabel("log10 (kcat * [E])\n([trna-aa]/s)")
		ax.set_ylabel("log10 (polypeptide elongation rate)\n([aa]/s)")

		for i, synthetaseId in enumerate(synthetaseIds):
			aaId = kinetic_data[synthetaseId]["amino acid"]

			# Compute vmax
			kcats = np.array([kinetic_data[synthetaseId][x].asNumber() for x in ["kcat AA", "kcat tRNA", "kcat ATP"]])
			synthetaseConcentrations_molar = [x.asNumber() for x in synthetaseConcentrations[:, i]]
			vmaxs = kcats * np.average(synthetaseConcentrations_molar)

			# Get elongation rate
			aaIndex = aaIds.index(aaId)
			elongationRate = np.average([x.asNumber() for x in elongationRates[:, aaIndex]])

			# Plot log values
			xlog = np.log10(max(vmaxs))
			ylog = np.log10(elongationRate)
			ax.scatter(xlog, ylog)
			ax.text(xlog, ylog, aaId[:-3], ha = "center")

		# Adjust axes ranges
		axmin = min(ax.get_xlim()[0], ax.get_ylim()[0])
		axmax = max(ax.get_xlim()[1], ax.get_ylim()[1])
		ax.set_xlim([axmin, axmax])
		ax.set_ylim([axmin, axmax])

		# Plot diagonal
		ax.plot([axmin, axmax], [axmin, axmax])
		plt.subplots_adjust(bottom = 0.2, top = 0.8)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
