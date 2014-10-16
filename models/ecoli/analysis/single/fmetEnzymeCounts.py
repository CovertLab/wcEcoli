#!/usr/bin/env python
"""
Plot counts of enzymes needed for methionine formylation, deformylation, and cleavage

@author: Kalli Kappel
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/15/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants

from wholecell.containers.unique_molecules_data import UniqueMoleculesData

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf"))

	names = h.root.names
	bulkMolecules = h.root.BulkMolecules

	moleculeIds = names.moleculeIDs.read()
	MTFId = "EG11268-MONOMER[c]"
	PDFId = "EG11440-MONOMER[c]"
	MAPId = "EG10570-MONOMER[c]"
	MTFIndex = moleculeIds.index(MTFId)
	PDFIndex = moleculeIds.index(PDFId)
	MAPIndex = moleculeIds.index(MAPId)
	MTFCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, MTFIndex]
	PDFCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, PDFIndex]
	MAPCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, MAPIndex]
	time = bulkMolecules.col("time")

	h.close()

	plt.figure(figsize = (8.5, 11))

	plt.plot(time/60.,MTFCountsBulk,label="MTF Counts")
	plt.plot(time/60.,PDFCountsBulk,label="PDF Counts")
	plt.plot(time/60.,MAPCountsBulk,label="MAP Counts")
	plt.legend(loc="upper left")
	#plt.axis([0,60,0,25])
	plt.xlabel("Time (min)")
	plt.ylabel("Counts of FMET Enzymes")
	plt.title("Counts of FMET Enzymes")

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
