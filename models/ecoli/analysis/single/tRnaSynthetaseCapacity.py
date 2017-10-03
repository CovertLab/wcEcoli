#!/usr/bin/env python
"""
Analyzes kinetic capacity of tRNA synthetases.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/26/2017
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
PLOT_VMAX = False

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
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

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))
	nAvogadro = sim_data.constants.nAvogadro
	cellDensity = sim_data.constants.cellDensity

	# Load data from sim
	cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass") * units.fg
	cellVolume = cellMass / cellDensity
	timestep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

	isTRna = sim_data.process.transcription.rnaData["isTRna"]
	trnaIds = sim_data.process.transcription.rnaData["id"][isTRna].tolist()
	aaIds = sim_data.moleculeGroups.aaIDs
	atpId = "ATP[c]"

	# Remove selenocysteine
	selenocysteine_index = aaIds.index("L-SELENOCYSTEINE[c]")
	aaIds.remove("L-SELENOCYSTEINE[c]")
	trnaIds.remove("selC-tRNA[c]")

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	synthetaseIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in synthetaseIds], np.int)
	trnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in trnaIds], np.int)
	aaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in aaIds], np.int)
	atpIndex = np.array(moleculeIds.index(atpId), np.int)

	synthetaseCounts = bulkMolecules.readColumn("counts")[:, synthetaseIndexes]
	trnaCounts = bulkMolecules.readColumn("counts")[:, trnaIndexes]
	aaCounts = bulkMolecules.readColumn("counts")[:, aaIndexes]
	atpCounts = bulkMolecules.readColumn("counts")[:, atpIndex]
	bulkMolecules.close()

	# Load number of incorporations of each amino acid from sim
	aasUsed = TableReader(os.path.join(simOutDir, "GrowthLimits")).readColumn("aasUsed")
	aasUsed = np.delete(aasUsed, selenocysteine_index, 1)

	# Compute concentrations
	synthetaseConcentrations = (synthetaseCounts / nAvogadro) / cellVolume[:, None]
	trnaConcentrations = (trnaCounts / nAvogadro) / cellVolume[:, None]
	aaConcentrations = (aaCounts / nAvogadro) / cellVolume[:, None]
	atpConcentrations = (atpCounts / nAvogadro) / cellVolume


	# Plot
	fig, axesList = plt.subplots(4, 5, figsize = (14, 10))
	axesList = axesList.reshape(-1)

	substrateList = ["tRNA", "AA", "ATP"]
	substrateColor = ["b", "g", "r"]

	for synthetaseId in synthetaseIds:
		print synthetaseId
		_aaId = kinetic_data[synthetaseId]["amino acid"]
		_aaIndex = aaIds.index(_aaId)
		_trnaIds = aa_to_trna[_aaId]
		_trnaIndexes = [trnaIds.index(x) for x in _trnaIds]

		ax = axesList[_aaIndex]

		E = synthetaseConcentrations[:, synthetaseIds.index(synthetaseId)]
		AA = aaConcentrations[aaIds.index(_aaId)]
		T = trnaConcentrations[:, _trnaIndexes].sum(axis = 1)

		for i, substrate in enumerate(substrateList):
			kcat = kinetic_data[synthetaseId]["kcat %s" % substrate]
			kM = (kinetic_data[synthetaseId]["kM %s" % substrate])
			if kcat.asNumber() == 0 or kM.asNumber() == 0:
				continue

			if PLOT_VMAX:
				vmax = E*kcat
				charged_tRNAs_max = np.floor(np.array(vmax * timestep*units.s * cellVolume * nAvogadro, dtype = float))
				ax.plot(time[SKIP_FIRST_N_TIMESTEPS:], charged_tRNAs_max[SKIP_FIRST_N_TIMESTEPS:] / aasUsed[SKIP_FIRST_N_TIMESTEPS:, _aaIndex], substrateColor[i])
			else:
				v = E*T / (T + kM.asUnit(units.mol/units.L)) * kcat
				charged_tRNAs = np.floor(np.array(v * timestep*units.s * cellVolume * nAvogadro, dtype = float))
				ax.plot(time[SKIP_FIRST_N_TIMESTEPS:], charged_tRNAs[SKIP_FIRST_N_TIMESTEPS:] / aasUsed[SKIP_FIRST_N_TIMESTEPS:, _aaIndex], substrateColor[i])

	for i, ax in enumerate(axesList):
		ax.set_xlim([time[SKIP_FIRST_N_TIMESTEPS], time[-1]])
		whitePadSparklineAxis(ax)
		if ax.get_ylim()[1] > 1:
			ax.set_yticks([0, 1, ax.get_ylim()[1]])
		else:
			ax.set_yticks([0, ax.get_ylim()[0], ax.get_ylim()[1], 1])
		ax.set_title(aaIds[i], fontsize = 12)
	
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.85, right = 0.95)
	from wholecell.analysis.analysis_tools import exportFigure

	if SKIP_FIRST_N_TIMESTEPS != 0:
		text = "skipping first %s time steps" % SKIP_FIRST_N_TIMESTEPS
		name_tag = "skip_%s_timesteps" % SKIP_FIRST_N_TIMESTEPS
	if PLOT_VMAX:
		fig.suptitle("reaction rate = vmax only\n(blue = tRNA)  (green = amino acid)  (red = ATP)\ny-axis = ratio of # charged-tRNAs to # incorporations into polypeptides\nx-axis = time (s)\n%s" % text)
		exportFigure(plt, plotOutDir, plotOutFileName + "__vmax_only__%s" % name_tag, metadata)
	else:
		fig.suptitle("reaction rate = vmax * ([substrate] / [substrate] + kM)\n(blue = tRNA)  (green = amino acid)  (red = ATP)\ny-axis = ratio of # charged-tRNAs to # incorporations into polypeptides\nx-axis = time (s)\n%s" % text)
		exportFigure(plt, plotOutDir, plotOutFileName + "__%s" % name_tag, metadata)

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
