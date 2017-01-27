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

	fig, axesList = plt.subplots(5,1)

	for simDir in ap.get_cells():
		simOutDir = os.path.join(simDir, "simOut")
		
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
		rnaMass = massDataFile.readColumn("rnaMass")
		proteinMass = massDataFile.readColumn("proteinMass")
		smallMoleculeMass = massDataFile.readColumn("smallMoleculeMass")
		cellMass = massDataFile.readColumn("cellMass")
		massDataFile.close()


		rna_protein_ratio = rnaMass / proteinMass

		axesList[0].plot(time / 60., rna_protein_ratio)
		axesList[1].plot(time / 60., rnaMass / cellMass)
		axesList[2].plot(time / 60., proteinMass / cellMass)
		axesList[3].plot(time / 60., smallMoleculeMass / cellMass)
		axesList[4].plot(time / 60., (rnaMass + proteinMass + smallMoleculeMass) / cellMass)

	axesList[4].set_xlabel("Time (min)")

	axesList[0].set_ylabel("RNA/Protein\n(mass/mass)")
	axesList[1].set_ylabel("RNA mass\nfraction")
	axesList[2].set_ylabel("Protein mass\nfraction")
	axesList[3].set_ylabel("Small molec mass\nfraction")
	axesList[4].set_ylabel("Dry mass\nfraction")

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
