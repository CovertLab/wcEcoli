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

from wholecell.utils.sparkline import whitePadSparklineAxis

from scipy.stats import linregress
def mm2inch(value):
	return value * 0.0393701

FONT_SIZE=9
trim = 0.05

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return r_value**2

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	
	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	if ap.n_generation == 1:
		print "Need more data to create addedMass"
		return

	fig = plt.figure()
	fig.set_figwidth(8)
	fig.set_figheight(4)
	ax0 = plt.subplot2grid((1,1), (0,0))

	label = [r"$\tau=$"+"44 min", r"$\tau=$"+"100 min", r"$\tau=$"+"22 min"]
	color = ["red", "blue", "green"]

	for varIdx in range(ap.n_variant):

		initial_masses = np.zeros(0)
		final_masses = np.zeros(0)

		if varIdx == 0:
			gen = [1,2,3]
		elif varIdx == 1:
			gen = [1,2,3]
		elif varIdx == 2:
			gen = [1,2,3]

		all_cells = ap.get_cells(generation=gen, variant=[varIdx])

		for simDir in all_cells:
			simOutDir = os.path.join(simDir, "simOut")
			mass = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = mass.readColumn("dryMass")

			initial_masses = np.hstack((initial_masses, cellMass[0]))
			final_masses = np.hstack((final_masses, cellMass[-1]))

		print final_masses.size
		added_masses = final_masses - initial_masses

		scaled_initial_masses = initial_masses / initial_masses.mean()
		scaled_added_masses = added_masses / added_masses.mean()

		n, bins, patches = ax0.hist(initial_masses, bins=np.arange(110,1000, 10),color=color[varIdx], alpha = 0.8, histtype="stepfilled", label=label[varIdx])
		ax0.set_xlim([110, 1000])

		# ax0.set_title("n = {}".format(n_cells))
		ax0.set_xlabel("Initial mass (fg)", fontsize=FONT_SIZE)
		ax0.set_ylabel("Frequency", fontsize=FONT_SIZE)

		plt.subplots_adjust(left = 0.2, bottom = 0.2, wspace= 0.6)

		whitePadSparklineAxis(ax0)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 

	ax0.legend(fontsize=FONT_SIZE)
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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
