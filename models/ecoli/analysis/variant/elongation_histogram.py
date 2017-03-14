#!/usr/bin/env python
"""
Plots rRNA control loop parameters

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/17/2017
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

FONT_SIZE = 7

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	allCells = ap.get_cells()

	# sim_data = cPickle.load(open(ap.get_variant_kb(0), "rb"))

	# fig = plt.figure(figsize = (8.5, 11))


	for var_idx in range(ap.n_variant):
		# sim_data = cPickle.load(open(ap.get_variant_kb(var_idx), "rb"))
		allCells = ap.get_cells(variant=[var_idx])

		net_variant_histogram = np.zeros(21, dtype=np.int)
		net_variant_noterm_histogram = np.zeros(21, dtype=np.int)

		for simDir in allCells:
			simOutDir = os.path.join(simDir, "simOut")

			# Load time
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

			# Load control loop data
			actualElongationHistogram = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("actualElongationHist")
			elongationsNonTerminatingHist = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("elongationsNonTerminatingHist")
			effective_elongation_rate = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate")

			net_variant_histogram += actualElongationHistogram.sum(axis=0)
			net_variant_noterm_histogram += elongationsNonTerminatingHist.sum(axis=0)

		ax0 = plt.subplot2grid((ap.n_variant,1), (var_idx,0))
		ax0.bar(np.arange(0,21) , net_variant_histogram / net_variant_histogram.sum())

		hist_avg = np.dot(net_variant_histogram / net_variant_histogram.sum().astype(np.float), np.arange(0,21))
		ax0.axvline(hist_avg, linewidth=2)
		ax0.set_xlabel("Number of elongations ribosome accomplished during timestep")




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
