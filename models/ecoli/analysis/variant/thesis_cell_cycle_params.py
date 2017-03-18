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


	bremer_tau = [40, 100, 24]

	bremer_origin_per_cell = [3.36,1.96,  6.54]
	bremer_terminus_per_cell = [3.64,1.23,  9.19]
	bremer_fork_per_cell = [3.64,1.46,  9.19]
	origins_per_cell_at_initiation = [2, 1, 4]


	sim_origin_per_cell = np.zeros(ap.n_variant)
	sim_terminus_per_cell = np.zeros(ap.n_variant)
	sim_fork_per_cell = np.zeros(ap.n_variant)
	sim_origins_per_cell_at_initiation = np.zeros(ap.n_variant)

	sim_origin_per_cell_std = np.zeros(ap.n_variant)
	sim_terminus_per_cell_std = np.zeros(ap.n_variant)
	sim_fork_per_cell_std = np.zeros(ap.n_variant)
	sim_origins_per_cell_at_initiation_std = np.zeros(ap.n_variant)

	for varIdx in range(ap.n_variant):

		all_cells = ap.get_cells(variant=[varIdx], seed=[0], generation=[1,2])
		import cPickle
		sim_data = cPickle.load(open(ap.get_variant_kb(varIdx)))

		mean_num_origin = np.zeros(len(all_cells))
		mean_num_terc = np.zeros(len(all_cells))
		mean_num_forks = np.zeros(len(all_cells))
		num_origin_at_init = np.zeros(len(all_cells))

		for idx, simDir in enumerate(all_cells):
			simOutDir = os.path.join(simDir, "simOut")
			
			numOrigin = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("numberOfOric")
			mean_num_origin[idx] = numOrigin.mean()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeIds = bulkMolecules.readAttribute("objectNames")
			fullChromIndex = moleculeIds.index(sim_data.moleculeGroups.fullChromosome[0])
			numTerC = bulkMolecules.readColumn("counts")[:, fullChromIndex]
			mean_num_terc[idx] = numTerC.mean()


			sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
			pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4
			mean_num_forks[idx] = pairsOfForks.mean()


			massPerOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")
			idxInit = np.where(massPerOric >= 1)[0]
			numOriginAtInit = numOrigin[idxInit]
			num_origin_at_init[idx] = numOriginAtInit


		sim_origin_per_cell[varIdx] = mean_num_origin.mean()
		sim_terminus_per_cell[varIdx] = mean_num_terc.mean()
		sim_fork_per_cell[varIdx] = mean_num_forks.mean()
		sim_origins_per_cell_at_initiation[varIdx] = num_origin_at_init.mean()

		sim_origin_per_cell_std[varIdx] = mean_num_origin.std()
		sim_terminus_per_cell_std[varIdx] = mean_num_terc.std()
		sim_fork_per_cell_std[varIdx] = mean_num_forks.std()
		sim_origins_per_cell_at_initiation_std[varIdx] = num_origin_at_init.std()


	import ipdb; ipdb.set_trace()




	# 	axes_list = [ax0]

	# 	for a in axes_list:
	# 		for tick in a.yaxis.get_major_ticks():
	# 			tick.label.set_fontsize(FONT_SIZE) 
	# 		for tick in a.xaxis.get_major_ticks():
	# 			tick.label.set_fontsize(FONT_SIZE) 

	# 	ax0.set_xlabel("Rounds of chromosome replication\ninitated per cell cycle", fontsize=FONT_SIZE)
	# 	ax0.set_title(title_list[2] + r", $n_{cells}=$" + "{}".format(len(all_cells)), fontsize=FONT_SIZE)

	# 	whitePadSparklineAxis(ax0)

	# plt.subplots_adjust(bottom = 0.2, wspace=0.3)

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
