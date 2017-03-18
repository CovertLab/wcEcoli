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
	fig.set_figheight(5)

	title_list = ["Glucose minimal anaerobic\n" + r"$\mu = $" + "100 min", "Glucose minimal\n" + r"$\mu = $" + "44 min", "Glucose minimal + 20 amino acids\n" + r"$\mu = $" + "22 min"]

	for varIdx in [2]:

		all_cells = ap.get_cells(variant=[varIdx], generation=[2,3,4,5,6,7])

		ax0 = plt.subplot2grid((1,2), (0,1))

		round_init = np.zeros(len(all_cells))

		for idx, simDir in enumerate(all_cells):
			simOutDir = os.path.join(simDir, "simOut")
			
			## Initiation mass
			massPerOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")
			round_init[idx] = np.sum(massPerOric >= 1.)

		#norm_round_init = np.bincount(round_init.astype(np.int))
		#import ipdb; ipdb.set_trace()
		ax0.hist(round_init, bins=[0,1,2,3], normed=True)

		axes_list = [ax0]

		for a in axes_list:
			for tick in a.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE) 
			for tick in a.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE) 

		ax0.set_xlabel("Rounds of chromosome replication\ninitated per cell cycle", fontsize=FONT_SIZE)
		ax0.set_title(title_list[2] + r", $n_{cells}=$" + "{}".format(len(all_cells)), fontsize=FONT_SIZE)

		whitePadSparklineAxis(ax0)

	plt.subplots_adjust(bottom = 0.2, wspace=0.3)

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
