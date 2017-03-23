#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.patches as patches


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from wholecell.utils.sparkline import whitePadSparklineAxis


PLACE_HOLDER = -1

FONT_SIZE=10
trim = 0.03
WIDTH = 100

# def sparklineAxis(axis):
# 	axis.spines['top'].set_visible(False)
# 	axis.spines['bottom'].set_visible(False)
# 	axis.xaxis.set_ticks_position('none')
# 	axis.tick_params(which = 'both', direction = 'out')


def mm2inch(value):
	return value * 0.0393701

def remove_artifacts(a):
	median = np.median(a)
	std = np.nanstd(a)
	a[a < median - 1.5*std] = np.nan
	a[a > median + 1.5*std] = np.nan

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	
	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	fig = plt.figure()
	fig.set_figwidth(8.5)
	fig.set_figheight(5)

	title_list = ["Induced variability in growth rate", "No induced variablity in growth rate"]

	for varIdx in range(ap.n_variant):

		all_cells = ap.get_cells(variant=[varIdx], seed=[2])

		import cPickle
		simDataFile = ap.get_variant_kb(all_cells[0])
		sim_data = cPickle.load(open(simDataFile, "rb"))


		ax0 = plt.subplot2grid((2,ap.n_variant), (0,varIdx))
		ax1 = plt.subplot2grid((2,ap.n_variant), (1,varIdx), sharex=ax0)


		for idx, simDir in enumerate(all_cells):
			color = "black"
			alpha = 0.8
			if idx % 2:
				color = "blue"
				blue = 0.8

			simOutDir = os.path.join(simDir, "simOut")

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			
			## Cell mass
			mass = TableReader(os.path.join(simOutDir, "Mass"))
			cellMass = mass.readColumn("cellMass")
			ax0.plot(time / 60., cellMass, color = color, alpha = alpha, linewidth=1)

			# Instantanious growth rate
			growthRate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
			growthRate = growthRate * 60.
			remove_artifacts(growthRate)
			ax1.plot(time / 60., growthRate, color = color, alpha = alpha, linewidth=1)
 

		ax0.set_xlim([0., time.max() / 60.])
		if varIdx == 0:
			ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax0.yaxis.set_visible(False)
		#ax0.set_ylim([400, 6000.])
		ax0.set_title(title_list[varIdx], fontsize=FONT_SIZE)
		# ax0.xaxis.set_visible(False)
		#ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

		# nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
		# try:
		# 	T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0] / 60.
		# except IndexError:
		# 	print "No shift detected"
		# 	T_ADD_AA = 1e9

		axes_list = [ax0, ax1]
		# for a in axes_list:
		# 	shift_time = T_ADD_AA
		# 	width = a.get_xlim()[1] - shift_time
		# 	height = a.get_ylim()[1] - a.get_ylim()[0]
		# 	a.add_patch(
		# 		patches.Rectangle(
		# 			(shift_time, a.get_ylim()[0]),   # (x,y)
		# 			width,          # width
		# 			height,          # height
		# 			alpha = 0.25,
		# 			color = "gray",
		# 			linewidth = 0.
		# 		)
		# 	)

		for a in axes_list:
			for tick in a.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE) 
			for tick in a.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE) 

		if varIdx == 0:
			ax1.set_ylabel("Instantanious\ngrowth rate\n(fg/fg-min)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax1.yaxis.set_visible(False)
		ax1.xaxis.set_visible(False)
		#ax1.set_ylim([0.004, 0.03])


		whitePadSparklineAxis(ax0, False)
		whitePadSparklineAxis(ax1)

	plt.subplots_adjust(left = 0.18, wspace=0.3, right=0.92)

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
