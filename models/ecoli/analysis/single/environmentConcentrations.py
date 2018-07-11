#!/usr/bin/env python
"""
Plots environment nutrient concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cPickle
import math
import itertools

from wholecell.analysis.plotting_tools import COLORS_LARGE

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

	# Load mass fractions
	mass = TableReader(os.path.join(simOutDir, "Mass"))

	protein = mass.readColumn("proteinMass")
	tRna = mass.readColumn("tRnaMass")
	rRna = mass.readColumn("rRnaMass")
	mRna = mass.readColumn("mRnaMass")
	dna = mass.readColumn("dnaMass")
	smallMolecules = mass.readColumn("smallMoleculeMass")

	masses = np.vstack([
		protein/protein[0],
		rRna/rRna[0],
		tRna/tRna[0],
		mRna/mRna[0],
		dna/dna[0],
		smallMolecules/smallMolecules[0],
		]).T
	massLabels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA", "Small Mol.s"]

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load environment data
	environment = TableReader(os.path.join(simOutDir, "Environment"))
	nutrient_names = environment.readAttribute("objectNames")
	nutrient_concentrations = environment.readColumn("nutrientConcentrations")
	# volume = environment.readColumn("volume")
	environment.close()


	# Build a mapping from nutrient_name to color
	idToColor = {}
	for nutrient_name, color in itertools.izip(nutrient_names, itertools.cycle(COLORS_LARGE)):
		idToColor[nutrient_name] = color

	fig = plt.figure(figsize = (17, 18))

	ax1 = fig.add_subplot(3, 1, 1)  # The big subplot
	ax2 = fig.add_subplot(3, 1, 2)  # The big subplot
	# ax2_1 = fig.add_subplot(4, 1, 3)
	# ax2_2 = fig.add_subplot(4, 1, 4)

	# mass fractions
	ax1.plot(time, masses, linewidth=2)
	ax1.set_xlabel("Time (min)")
	ax1.set_ylabel("Mass (normalized by t = 0)")
	ax1.set_title("Biomass components")
	ax1.legend(massLabels, loc="best")


	## environment concentrations
	# log below
	for idx, nutrient_name in enumerate(nutrient_names):
		if (not math.isnan(nutrient_concentrations[0, idx]) and np.mean(nutrient_concentrations[:, idx])!=0):
			ax2.plot(time, nutrient_concentrations[:,idx], linewidth=2, label=nutrient_name, color=idToColor[nutrient_name])

	ax2.set_yscale('log')
	ymin, ymax = ax2.get_ylim()
	ax2.set_ylim((ymin, 1))
	ax2.spines['top'].set_visible(False)
	ax2.xaxis.set_ticks_position('bottom')

	divider = make_axes_locatable(ax2)
	axLin = divider.append_axes("top", size="100%", pad=0, sharex=ax2)

	# linear above
	for idx, nutrient_name in enumerate(nutrient_names):
		if (not math.isnan(nutrient_concentrations[0, idx]) and np.mean(nutrient_concentrations[:, idx])!=0):
			axLin.plot(time, nutrient_concentrations[:,idx], linewidth=2, label=nutrient_name, color=idToColor[nutrient_name])

	axLin.set_xscale('linear')
	axLin.set_ylim((1, ymax))

	# Remove bottom axis line
	axLin.spines['bottom'].set_visible(False)
	axLin.xaxis.set_ticks_position('top')
	plt.setp(axLin.get_xticklabels(), visible=False)

	ax2.set_title('environment concentrations -- Linear above, log below')
	ax2.set_xlabel('Time (sec)')
	ax2.set_ylabel('concentration (mmol/L)')

	# place legend
	ax2.legend(bbox_to_anchor=(0.5, -0.25), loc=9, borderaxespad=0., ncol=3, prop={'size': 10})

	plt.subplots_adjust(wspace=0.2, hspace=0.4)






	#
	# for idx, nutrient_name in enumerate(nutrient_names):
	# 	if (not math.isnan(nutrient_concentrations[0, idx]) and np.mean(nutrient_concentrations[:, idx])!=0):
	# 		# Unadjusted
	# 		plt.subplot(4, 1, 2)
	# 		plt.plot(time, nutrient_concentrations[:,idx], linewidth=2, label=nutrient_name, color=idToColor[nutrient_name])
	#
	# 		# Log scale
	# 		plt.subplot(4, 1, 3)
	# 		plt.plot(time, np.log10(nutrient_concentrations[:,idx]), linewidth=2, label=nutrient_name, color=idToColor[nutrient_name])
	#
	#




	#
	#
	#
	# plt.subplot(4, 1, 1)
	# plt.xlabel("Time (min)")
	# plt.ylabel("Mass (normalized by t = 0)")
	# plt.title("Biomass components")
	# plt.legend(massLabels, loc="best")
	#
	# plt.subplot(4, 1, 2)
	# plt.title('environment concentrations')
	# plt.xlabel('Time (sec)')
	# plt.ylabel('concentration (mmol/L)')
	#
	# plt.subplot(4, 1, 3)
	# plt.title('environment log(concentrations)')
	# plt.xlabel('Time (sec)')
	# plt.ylabel('Log10 concentration (mmol/L)')

	# # place legend
	# plt.legend(bbox_to_anchor=(0.5, -0.25), loc=9, borderaxespad=0., ncol=3, prop={'size': 10})
	#
	# plt.subplots_adjust(wspace=0.2, hspace=0.4)

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
