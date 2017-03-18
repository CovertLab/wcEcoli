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


# def sparklineAxis(axis):
# 	axis.spines['top'].set_visible(False)
# 	axis.spines['bottom'].set_visible(False)
# 	axis.xaxis.set_ticks_position('none')
# 	axis.tick_params(which = 'both', direction = 'out')


def mm2inch(value):
	return value * 0.0393701

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	
	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	fig = plt.figure()
	fig.set_figwidth(5)
	fig.set_figheight(4)

	for varIdx in [0]:

		all_cells = ap.get_cells(variant=[varIdx], seed=[2], generation=[0,1,2,3,4])

		import cPickle
		simDataFile = ap.get_variant_kb(all_cells[0])
		sim_data = cPickle.load(open(simDataFile, "rb"))

		oriC = sim_data.constants.oriCCenter.asNumber()
		terC = sim_data.constants.terCCenter.asNumber()
		genomeLength = len(sim_data.process.replication.genome_sequence)

		ax0 = plt.subplot2grid((3,1), (0,varIdx))
		ax2 = plt.subplot2grid((3,1), (1,varIdx), sharex=ax0)
		ax4 = plt.subplot2grid((3,1), (2,varIdx), sharex=ax0)

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
			massPerOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")
			idxInit = np.where(massPerOric >= 1)[0]
			ax0.plot(time / 60., cellMass, color = color, alpha = alpha, linewidth=1)
			ax0.plot(time[idxInit] / 60., cellMass[idxInit],  markersize=3, linewidth=0, marker="o", color = "red", markeredgewidth=0)

			## Number of oriC
			numberOfOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("numberOfOric")
			ax2.plot(time / 60., numberOfOric, color = color, alpha = alpha, linewidth=1)


			## Fork position
			sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
			# Pairs of forks
			pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4
			ax4.plot(time / 60., pairsOfForks, linewidth=1, color = color, alpha = alpha)
			ax4.set_yticks(np.arange(0,7))
			ax4.set_ylim([0, 6])


		# y_ticks = ax0.get_yticks()
		# new_y_ticks = y_ticks[0:-1:2]
		# ax0.set_yticks(new_y_ticks)

		ax0.set_xlim([0., time.max() / 60.])
		if varIdx == 0:
			ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax0.yaxis.set_visible(False)
		ax0.set_ylim([400, 6000.])
		ax0.set_title("Glucose minimal\n" + r"$\mu = $" + "44 min", fontsize=FONT_SIZE)
		# ax0.xaxis.set_visible(False)
		#ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

		# nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
		# try:
		# 	T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0] / 60.
		# except IndexError:
		# 	print "No shift detected"
		# 	T_ADD_AA = 1e9

		axes_list = [ax0, ax2, ax4]
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
			ax2.set_ylabel("Origin of\nrepliction\nper cell", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax2.yaxis.set_visible(False)
		ax2.xaxis.set_visible(False)
		ax2.set_ylim([1, 8])

		y_ticks = ax2.get_yticks()
		new_y_ticks = y_ticks[0:-1:2]
		ax2.set_yticks(new_y_ticks)

		if varIdx == 0:
			ax4.set_ylabel("Relative rate\nof dNTP\npolymerization", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax4.yaxis.set_visible(False)
		ax4.set_xlabel("Time (min)", fontsize=FONT_SIZE)
		# ax4.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

		whitePadSparklineAxis(ax0, False)
		whitePadSparklineAxis(ax2, False)
		whitePadSparklineAxis(ax4)

	plt.subplots_adjust(left = 0.22, wspace=0.3)

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
