#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))
	ap = AnalysisPaths(variantDir, cohort_plot = True)

	n_cells = ap.get_cells().size
	growth_rates = np.zeros(0)
	doubling_times = np.zeros(0)
	ribConcs = np.zeros(0)

	for simDir in ap.get_cells():
		simOutDir = os.path.join(simDir, "simOut")

		growthRate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
		growthRate = (1/units.s) * growthRate
		growth_rates = np.hstack((growth_rates, np.nanmean(growthRate.asNumber(1/units.min))))
		doubling_time = (1 / growthRate * np.log(2) ).asNumber(units.min)
		doubling_times = np.hstack((doubling_times, np.nanmean(doubling_time)))

		massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
		rnaMass = massDataFile.readColumn("rnaMass")
		proteinMass = massDataFile.readColumn("proteinMass")
		cellMass = massDataFile.readColumn("cellMass")
		massDataFile.close()

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
		uniqueMoleculeCounts.close()

		ribosomeCounts = activeRibosome
		ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass))
		ribConcs = np.hstack((ribConcs, np.nanmean(ribosomeConcentration.asNumber(units.umol/units.L))))

	fig, axesList = plt.subplots(2)

	axesList[0].plot(ribConcs, doubling_times, 'x')
	axesList[1].set_xlabel("Average [Ribosomes] (umol/L)")
	axesList[1].set_ylabel("Average growth rate (1/min)")

	axesList[1].plot(ribConcs, growth_rates, 'x')
	axesList[0].set_ylabel("Average doubling time (min)")

	z = np.polyfit(ribConcs, growth_rates, 1)
	p = np.poly1d(z)
	axesList[1].plot(ribConcs, p(ribConcs), '--')
	text_x = np.mean(axesList[1].get_xlim())
	text_y = np.mean(axesList[1].get_ylim()) + np.mean(axesList[1].get_ylim())*0.1
	axesList[1].text(text_x, text_y, r"$mu$={}$\times$$rib$ + {}".format(z[0],z[1]))

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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
