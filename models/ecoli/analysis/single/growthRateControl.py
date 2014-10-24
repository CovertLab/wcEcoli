#!/usr/bin/env python
"""
Plots counts of things assosciated with growth rate control via ppGpp

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/24/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants

FONT = {
		'size'	:	8
		}

def setAxisMaxMin(axis, data):
	ymax = np.max(data)
	ymin = 0
	if ymin == ymax:
		axis.set_yticks([ymin])
	else:
		axis.set_yticks([ymin, ymax])

def sparklineAxis(axis, x, y, tickPos, lineType, color):
	axis.plot(x, y, linestyle = 'steps' + lineType, color = color, linewidth = 2)
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.yaxis.set_ticks_position(tickPos)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')
	axis.tick_params(labelbottom = 'off')
	for tl in axis.get_yticklabels():
		tl.set_color(color)


def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))

	# List of IDs we are interested in
	relA = "RELA-MONOMER[c]"
	spoT = "SPOT-MONOMER[c]"
	ppGpp = "PPGPP[c]"
	bulkIds = [relA, spoT, ppGpp]
	# Other things we are intersted in
	# Ribosome mass fraction
	# Protein mass fraction
	# How much stable vs total RNA transcription is happening

	# Load counts of bulk molecules
	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:
		# Get indexes
		moleculeIds = bulkMoleculesFile.root.names.moleculeIDs.read()
		bulkIndexes = [moleculeIds.index(idx) for idx in bulkIds]

		# Load data
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules
		bulkTime = bulkMolecules.col("time")
		bulkCounts = bulkMolecules.read(0, None, 1, "counts")[:, bulkIndexes]

	with tables.open_file(os.path.join(simOutDir, "Mass.hdf")) as massFile:
		table = massFile.root.Mass
		cellDry = table.col("dryMass")
		protein = table.col("proteinMass")
		rna = table.col("rnaMass")
		tRna = table.col("tRnaMass")
		rRna = table.col("rRnaMass")
		mRna = table.col("mRnaMass")
		#dna = table.col("dnaMass")
		massTime = table.col("time")

	# Calculate mass fractions
	rRnaMassFraction = rRna/cellDry
	proteinMassFraction = protein/cellDry
	tRnaMassFraction = tRna/cellDry
	stableRnaMassFraction = (rRna + tRna)/cellDry

	# Calculate rates of transcription
	# stableRna = tRna + rRna
	# totalRnaTranscriptionRate = np.ediff1d(rna)
	# stableRnaTranscriptionRate = np.ediff1d(stableRna)

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	for idx in xrange(len(bulkIds)):
		protein_axis = plt.subplot(5, 3, idx + 1)

		sparklineAxis(protein_axis, bulkTime / 60., bulkCounts[:, idx], 'left', '-', 'b')
		setAxisMaxMin(protein_axis, bulkCounts[:, idx])

		# Component label
		protein_axis.set_xlabel(bulkIds[idx])

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

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
