#!/usr/bin/env python
"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/3/2017
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

SKIP_FIRST_N_TIMESTEPS = 10
LABEL = False

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	return
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"
	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load kinetic data from validation
	validation_data = cPickle.load(open(validationDataFile, "rb"))
	kinetic_data = validation_data.trna.synthetaseKinetics
	trna_pairings = validation_data.trna.trnaPairings
	aa_to_trna = validation_data.trna.aaToTrna
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
	synthetaseIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in synthetaseIds], np.int)
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
	fig, axesList = plt.subplots(1, 3, figsize = (14, 5))
	substrateList = ["tRNA", "AA", "ATP"]
	substrateColor = ["b", "g", "r"]
	ax1, ax2, ax3 = axesList

	# Determined by tRnaSynthetaseCapacity, hard-coded here
	low_supply_aas = ["ASN[c]", "L-ASPARTATE[c]", "GLT[c]", "GLN[c]", "LEU[c]", "LYS[c]", "MET[c]", "PHE[c]", "SER[c]", "TRP[c]"]

	for i, synthetaseId in enumerate(synthetaseIds):
		aa_id = kinetic_data[synthetaseId]["amino acid"]
		_aaIndex = aaIds.index(aa_id)
		E = synthetaseConcentrations[:, synthetaseIds.index(synthetaseId)]
		# E = np.average([x.asNumber() for x in synthetaseConcentrations[:, synthetaseIds.index(synthetaseId)]])
		for j, substrate in enumerate(substrateList):
			kcat = kinetic_data[synthetaseId]["kcat %s" % substrate]
			if kcat.asNumber() == 0:
				continue

			vmax = E*kcat
			charged_tRNAs_max = np.floor(np.array(vmax * timestep*units.s * cellVolume * nAvogadro, dtype = float))
			y = np.average(charged_tRNAs_max[SKIP_FIRST_N_TIMESTEPS:] / aasUsed[SKIP_FIRST_N_TIMESTEPS:, _aaIndex])
			x1 = kcat.asNumber(1/units.s)
			x2 = np.average([x.asNumber() for x in E]) * 1e6
			x3 = (np.average(aasUsed[:, _aaIndex]) * 1e-3)

			if aa_id in low_supply_aas:
				ax1.scatter(x1, y, facecolor = "orange", edgecolor = "none", s = 20)
				ax2.scatter(x2, y, facecolor = "orange", edgecolor = "none", s = 20)
				ax3.scatter(x3, y, facecolor = "orange", edgecolor = "none", s = 20)

				if LABEL:
					ax1.text(x1, y, aa_id, fontsize = 4)
					ax2.text(x2, y, aa_id, fontsize = 4)
					ax3.text(x3, y, aa_id, fontsize = 4)

			else:
				ax1.scatter(x1, y, facecolor = "k", alpha = 0.5, edgecolor = "none", s = 20)
				ax2.scatter(x2, y, facecolor = "k", alpha = 0.5, edgecolor = "none", s = 20)
				ax3.scatter(x3, y, facecolor = "k", alpha = 0.5, edgecolor = "none", s = 20)

				if LABEL:
					ax1.text(x1, y, aa_id, fontsize = 4)
					ax2.text(x2, y, aa_id, fontsize = 4)
					ax3.text(x3, y, aa_id, fontsize = 4)


	ax1.set_xlabel("kcat (1/s)")
	ax2.set_xlabel("[synthetase] (uM)")
	ax3.set_xlabel("demand (1000s amino acids)")
	ax1.set_ylabel("Ratio (supply / demand)")
	ax1.set_xlim([0, 100])
	ax2.set_xlim([0.2, 3.6])
	ax3.set_xlim([0, 32])
	for ax in axesList:
		ax.tick_params(right = "off")
		ax.tick_params(top = "off")
		ax.tick_params(which = "both", direction = "out")
		ax.set_ylim([0, 12])

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.85, right = 0.95)
	from wholecell.analysis.analysis_tools import exportFigure
	if LABEL:
		exportFigure(plt, plotOutDir, plotOutFileName + "__labeled", metadata)
	else:
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
