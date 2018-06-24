#!/usr/bin/env python
"""
Plots environment nutrient concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle
import math

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))
	nAvogadro = sim_data.constants.nAvogadro
	cellDensity = sim_data.constants.cellDensity

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load environment data
	environment = TableReader(os.path.join(simOutDir, "Environment"))
	nutrient_names = environment.readAttribute("objectNames")
	nutrient_concentrations = environment.readColumn("nutrientConcentrations")
	environment.close()


	fig = plt.figure(figsize = (17, 6))
	ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])

	for idx, nutrient_id in enumerate(nutrient_names):
		if (not math.isnan(nutrient_concentrations[0, idx]) and nutrient_concentrations[0, idx]>0):
			ax.plot(time, nutrient_concentrations[:,idx], linewidth=2, label=nutrient_id)

	ax.set_xlabel('Time (sec)')
	ax.set_ylabel('concentration (mmol/L)')

	# # Shrink current axis by 20%
	# box = plt.get_position()
	# plt.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	# Put a legend to the right of the current axis
	ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)



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

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], None)
