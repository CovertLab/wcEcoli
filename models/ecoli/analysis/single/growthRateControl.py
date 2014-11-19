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
from wholecell.utils import units

FONT = {'size'	:	8}
PPGPP_IDX = 2

def setAxisMaxMinY(axis, data):
	ymax = np.nanmax(data)
	ymin = 0
	if ymin == ymax:
		axis.set_yticks([ymin])
	else:
		axis.set_yticks([ymin, ymax])

def setAxisMaxMinX(axis, data):
	xmax = np.nanmax(data)
	xmin = 0
	if xmin == xmax:
		axis.set_xticks([xmin])
	else:
		axis.set_xticks([xmin, xmax])

def sparklineAxis(axis, x, y, tickPos, lineType, color):
	axis.plot(x, y, linestyle = 'steps' + lineType, color = color, linewidth = 2)
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.yaxis.set_ticks_position(tickPos)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')
	#axis.tick_params(labelbottom = 'off')
	for tl in axis.get_yticklabels():
		tl.set_color(color)


def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	## Load Data ##
	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	cellDensity = kb.cellDensity
	nAvogadro = kb.nAvogadro

	# List of IDs we are interested in
	relA = "RELA-MONOMER[c]"
	spoT = "SPOT-MONOMER[c]"
	ppGpp = "PPGPP[c]"
	bulkIds = [relA, spoT, ppGpp]
	# Other things we are intersted in
	# Ribosome mass fraction
	# Protein mass fraction
	# How much stable vs total RNA transcription is happening

	## LOAD DATA ##
	# Load counts of bulk molecules
	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:
		# Get indexes
		moleculeIds = bulkMoleculesFile.root.names.moleculeIDs.read()
		bulkIndexes = [moleculeIds.index(idx) for idx in bulkIds]

		# Load data
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules
		bulkTime = bulkMolecules.col("time")
		bulkCounts = bulkMolecules.read(0, None, 1, "counts")[:, bulkIndexes]

	# Load mass fractions
	with tables.open_file(os.path.join(simOutDir, "Mass.hdf")) as massFile:
		table = massFile.root.Mass
		cellMass = table.col("cellMass")
		growth = table.col("growth")
		cellDry = table.col("dryMass")
		protein = table.col("proteinMass")
		rna = table.col("rnaMass")
		tRna = table.col("tRnaMass")
		rRna = table.col("rRnaMass")
		massTime = table.col("time")
	proteinMassFraction = protein / cellDry
	rRnaMassFraction = rRna / cellDry
	tRnaMassFraction = tRna / cellDry
	rnaMassFraction = rna / cellDry
	rRnaFractionOfRna = rRna / rna
	growthRate = growth / cellDry * 3600
	N = 20
	growthRateRunningMean = np.zeros(growthRate.size)
	for i in range(growthRate.size):
		growthRateRunningMean[i] = np.mean(growthRate[i:i+N])
	dryMass = cellDry
	doubledDryMass = dryMass[0]*2 * np.ones(dryMass.size)

	# Load ratio of stable to total RNA synthesis
	with tables.open_file(os.path.join(simOutDir, "InitiatedTranscripts.hdf")) as massFile:
		table = massFile.root.InitiatedTranscripts
		ratioStableToToalInitalized = table.col("ratioStableToToalInitalized")
		initTime = table.col("time")
	N = 20
	runningMeanRatio = np.zeros(ratioStableToToalInitalized.size)
	for i in range(ratioStableToToalInitalized.size):
		runningMeanRatio[i] = np.mean(ratioStableToToalInitalized[i:i+N])

	# Load other growth rate control data
	with tables.open_file(os.path.join(simOutDir, "GrowthRateControl.hdf")) as massFile:
		table = massFile.root.GrowthRateControl
		totalStalls = table.col("totalStalls")
		synthetaseSaturation = table.col("synthetaseSaturation")
		spoTSaturation = table.col("spoT_saturation")
		gcTime = table.col("time")
	
	synthetaseSaturationMean = synthetaseSaturation.mean(axis = 1)
	synthetaseSaturationMax = synthetaseSaturation.max(axis = 1)
	synthetaseSaturationMin = synthetaseSaturation.min(axis = 1)

	## CALCULATE DATA ##
	# Calculate ppGpp concentration
	cellMass = units.fg * cellMass
	cellVolume = cellMass / cellDensity
	ppGppConc = ((1 / nAvogadro) * (1 / cellVolume) * bulkCounts[:, PPGPP_IDX]).asNumber(units.umol / units.L)

	## START PLOTTING ##
	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)
	NUMBER_ROWS = 9

	# Plot ppGpp concentration
	ppGppConc_axis = plt.subplot(NUMBER_ROWS, 1, 1)

	sparklineAxis(ppGppConc_axis, initTime / 60., ppGppConc, 'left', '-', 'b')
	setAxisMaxMinX(ppGppConc_axis, initTime / 60.)
	ppGppConc_axis.set_ylabel('[ppGpp] uM')

	# Plot ratio of stable to total rna synthesis
	fractionStable_axis = plt.subplot(NUMBER_ROWS, 1, 2)

	sparklineAxis(fractionStable_axis, initTime / 60., ratioStableToToalInitalized, 'left', '-', 'b')
	setAxisMaxMinX(fractionStable_axis, initTime / 60.)
	fractionStable_axis.set_ylabel('$r_s / r_t$')
	fractionStable_axis.plot(initTime / 60, runningMeanRatio, linestyle = '--', color = 'k', linewidth = 2)

	# Plot instantanious growth rate and dry mass
	growthRate_axis = plt.subplot(NUMBER_ROWS, 1, 3)
	sparklineAxis(growthRate_axis, massTime / 60., growthRate, 'left', '-', 'b')
	setAxisMaxMinX(growthRate_axis, massTime / 60.)
	growthRate_axis.set_ylabel('Growth rate\ngDCW/gDCW-hr')
	growthRate_axis.plot(massTime / 60, growthRateRunningMean, linestyle = '--', color = 'k', linewidth = 2)

	# Plot dry mass
	dryMass_axis = plt.subplot(NUMBER_ROWS, 1, 4)
	sparklineAxis(dryMass_axis, massTime / 60., dryMass, 'left', '-', 'b')
	setAxisMaxMinX(dryMass_axis, massTime / 60.)
	dryMass_axis.set_ylabel('Dry mass (g)')
	dryMass_axis.plot(massTime / 60, doubledDryMass, linestyle = '--', color = 'k', linewidth = 2)

	# Plot total stalls
	stall_axis = plt.subplot(NUMBER_ROWS, 1, 5)
	sparklineAxis(stall_axis, gcTime / 60., totalStalls, 'left', '-', 'b')
	setAxisMaxMinY(stall_axis, totalStalls)
	setAxisMaxMinX(stall_axis, gcTime / 60.)
	stall_axis.set_ylabel('Total stalls')

	# Saturation of synthetases
	synthetaseSat_axis = plt.subplot(NUMBER_ROWS, 1, 6)
	sparklineAxis(synthetaseSat_axis, gcTime / 60., synthetaseSaturationMin, 'left', '-', 'b')
	setAxisMaxMinY(synthetaseSat_axis, np.around(synthetaseSaturationMin, decimals=2))
	synthetaseSat_axis.set_ylabel('Synthetase\nsaturation\nminimum')

	# Plot spoT saturation
	spotSat_axis = plt.subplot(NUMBER_ROWS, 1, 7)
	sparklineAxis(spotSat_axis, gcTime / 60., spoTSaturation, 'left', '-', 'b')
	setAxisMaxMinY(spotSat_axis, np.around(spoTSaturation, decimals = 2))
	setAxisMaxMinX(spotSat_axis, gcTime / 60.)
	spotSat_axis.set_ylabel('spoT saturation')

	# Plot proteins if interest
	bulkIds.pop(bulkIds.index("PPGPP[c]"))
	ROW = 7
	for idx in xrange(len(bulkIds)):
		bulkObject_axis = plt.subplot(NUMBER_ROWS, 3, idx + 3*ROW+1)

		sparklineAxis(bulkObject_axis, bulkTime / 60., bulkCounts[:, idx], 'left', '-', 'b')
		setAxisMaxMinY(bulkObject_axis, bulkCounts[:, idx])
		setAxisMaxMinX(bulkObject_axis, bulkTime / 60.)

		# Component label
		bulkObject_axis.set_title(bulkIds[idx])

	# Plot mass fractions of RNA and protein
	ROW = 8
	RNA_axis = plt.subplot(NUMBER_ROWS, 2, ROW*2+1)
	sparklineAxis(RNA_axis, massTime / 60., rRnaFractionOfRna, 'left', '-', 'b')
	setAxisMaxMinX(RNA_axis, massTime / 60.)
	RNA_axis.set_yticks([0., 1.])
	RNA_axis.set_ylabel('rRNA fraction of total RNA')

	protein_axis = plt.subplot(NUMBER_ROWS, 2, ROW*2+2)
	sparklineAxis(protein_axis, massTime / 60., proteinMassFraction, 'left', '-', 'b')
	setAxisMaxMinX(protein_axis, massTime / 60.)
	protein_axis.set_yticks([0., 1.])
	protein_axis.set_ylabel('Protein mass fraction')

	## FORMATTING AND SAVING ##
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
