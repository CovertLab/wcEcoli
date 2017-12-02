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
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis

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
	timestep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

	# Load data from sim
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
	x_pos = np.arange(20)
	width = 0.25
	x_ticklabels = []

	fig, axesList = plt.subplots(3, 1, figsize = (14, 10))
	ax1 = axesList[0]
	ax2 = axesList[1]
	ax3 = axesList[2]
	substrateList = ["tRNA", "AA", "ATP"]
	substrateColor = ["b", "g", "r"]

	low_supply_aas = ["ASN[c]", "L-ASPARTATE[c]", "GLT[c]", "GLN[c]", "LEU[c]", "LYS[c]", "MET[c]", "PHE[c]", "SER[c]", "TRP[c]"]

	# order synthetases by their max kcat
	maxKcats = []
	for synthetaseId in synthetaseIds:
		kcats = [kinetic_data[synthetaseId]["kcat %s" % substrate] for substrate in substrateList]
		maxKcats.append(max(kcats).asNumber())

	synthetaseOrder = np.argsort(maxKcats)
	synthetaseIds = np.array(synthetaseIds)[synthetaseOrder]
	synthetaseIds = synthetaseIds.tolist()
	for i, synthetaseId in enumerate(synthetaseIds):
		aa_id = kinetic_data[synthetaseId]["amino acid"]
		_aaIndex = aaIds.index(aa_id)
		E = np.average([x.asNumber() for x in synthetaseConcentrations[:, synthetaseIds.index(synthetaseId)]])
		ax2.bar(x_pos[i], E * 1e6, width * 3, edgecolor = "none", facecolor = "k", alpha = 0.5)
		demand = np.average(aasUsed[:, _aaIndex])

		# calculate minimum kcat required to meet demands
		# kcat required = (# aa used) / (E * timestep * volume * avogadro's number)
		denominator = synthetaseConcentrations[:, synthetaseIds.index(synthetaseId)] * timestep * cellVolume * nAvogadro
		kcat_required = np.average([x.asNumber() for x in aasUsed[:, _aaIndex] / denominator])
		ax1.plot([x_pos[i], x_pos[i] + 3 * width], [kcat_required, kcat_required], color = "k", linestyle = "--", linewidth = 2)

		# calculate synthetase concentration required to meet demands
		# synthetase concentration required = (# aa used) / (kcat * timestep * volume * avogadro's number)
		kcat_max = np.max([kinetic_data[synthetaseId]["kcat %s" % substrate] for substrate in substrateList])
		denominator = timestep * cellVolume * nAvogadro * kcat_max
		E_required = np.average([x.asNumber() for x in aasUsed[:, _aaIndex] / denominator])
		ax2.plot([x_pos[i], x_pos[i] + width * 3], [E_required * 1e6, E_required * 1e6], color = "k", linestyle = "--", linewidth = 2)

		ax3.bar(x_pos[i], demand / 1000., width * 3, edgecolor = "none", facecolor = "k", alpha = 0.5)
		if aa_id in low_supply_aas:
			x_ticklabels.append(aa_id + "*")
		else:
			x_ticklabels.append(aa_id)

		for j, substrate in enumerate(substrateList):
			kcat = kinetic_data[synthetaseId]["kcat %s" % substrate]
			if kcat.asNumber() == 0:
				continue
			ax1.bar(x_pos[i] + j * width, kcat.asNumber(), width, edgecolor = "none", facecolor = substrateColor[j])

	ax1.set_title("kcat (1/s)")
	ax2.set_title("average [synthetase] (uM)")
	ax3.set_title("average demand per timestep (1000s amino acids)")
	blue_patch = mpatches.Patch(color = "blue", label = "tRNA")
	green_patch = mpatches.Patch(color = "green", label = "amino acid")
	red_patch = mpatches.Patch(color = "red", label = "ATP")
	black_dashed_line = mlines.Line2D([], [], color = "k", linestyle = "--", linewidth = 2, label = "kcat required")
	plt.axes(ax1)
	plt.legend(handles = [blue_patch, green_patch, red_patch, black_dashed_line],
		loc = "upper left",
		bbox_to_anchor = [0, 1])
	plt.axes(ax2)
	black_dashed_line = mlines.Line2D([], [], color = "k", linestyle = "--", linewidth = 2, label = "synthetase required (assumes max kcat value from literature)")
	plt.legend(handles = [black_dashed_line],
		loc = "upper right")

	for ax in axesList:
		ax.set_xticks(x_pos + 1.5 * width)
		ax.set_xticklabels(x_ticklabels, fontsize = 6)
		ax.tick_params(right = "off")
		ax.tick_params(top = "off")
		ax.tick_params(which = "both", direction = "out")

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.1, top = 0.85, right = 0.95)
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
