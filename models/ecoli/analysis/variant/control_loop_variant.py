#!/usr/bin/env python
"""
Plots rRNA control loop parameters

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/17/2017
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

FONT_SIZE = 7

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	allCells = ap.get_cells()

	sim_data = cPickle.load(open(ap.get_variant_kb(0), "rb"))

	fig = plt.figure(figsize = (8.5, 11))

	ax0 = plt.subplot2grid((6,2), (0,0))
	ax1 = plt.subplot2grid((6,2), (1,0))
	ax2 = plt.subplot2grid((6,2), (2,0))
	ax3 = plt.subplot2grid((6,2), (3,0))
	ax4 = plt.subplot2grid((6,2), (4,0))
	ax5 = plt.subplot2grid((6,2), (5,0))

	ax6 = plt.subplot2grid((6,2), (0,1), rowspan=2)
	ax7 = plt.subplot2grid((6,2), (3,1), rowspan=2)
	axesList = [ax0, ax1, ax2, ax3, ax4, ax5]

	ribConcMax = 0.
	ribConcMin = 1000.

	firstPass = True

	for var_idx in range(ap.n_variant):
		sim_data = cPickle.load(open(ap.get_variant_kb(var_idx), "rb"))
		allCells = ap.get_cells(variant=[var_idx])

		for simDir in allCells:
			label1 = None
			label2 = None

			simOutDir = os.path.join(simDir, "simOut")

			# Load time
			initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

			# Load control loop data
			errorInElongationRate = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("errorInElongationRate")
			proportionalTerm = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("proportionalTerm")
			integralTerm = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("integralTerm")
			rRnaSynthRate_expected = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_expected")
			rRnaSynthRate_updated = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_updated")
			bias = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_updated")

			growthRate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
			growthRate = (1 / units.s) * growthRate
			doublingTime = 1 / growthRate * np.log(2)

			effective_elongation_rate = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate")

			expected_doubling_time = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("expectedDoublingTime")
			setpoint_elongation_rate = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("expectedElongationRate")

			nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
			media = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][0][1]
			expected_elongation_rate = sim_data.process.translation.ribosomeElongationRateDict[media]


			cellMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			uniqueMoleculeCounts.close()
			ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass))

			if ribosomeConcentration.asNumber(units.mmol / units.L).max() > ribConcMax:
				ribConcMax = ribosomeConcentration.asNumber(units.mmol / units.L).max()
			if ribosomeConcentration.asNumber(units.mmol / units.L)[10:].min() < ribConcMin:
				ribConcMin = ribosomeConcentration.asNumber(units.mmol / units.L)[10:].min()

			# Plot stuff
			if firstPass:
				label1 = "raw"
				label2 = "bias"

			axesList[0].plot(time / 60., errorInElongationRate, label=label1, color = "blue", alpha = 0.6)
			axesList[0].plot(time / 60., np.zeros(time.size), '--')
			axesList[0].plot(time / 60., errorInElongationRate + bias, label=label2, color = "red", alpha=0.6)
			axesList[0].legend(fontsize=FONT_SIZE, loc=4,frameon=False)
			axesList[0].set_ylabel("Error " + r"$(e_{expected} - e_{actual})$", fontsize=FONT_SIZE)
			axesList[0].set_ylim([0, -14])

			# axesList[1].plot(time / 60., proportionalTerm, label = "proportional", alpha = 0.7)
			# axesList[1].plot(time / 60., integralTerm, label = "integral", alpha = 0.7)
			axesList[1].plot(time / 60., proportionalTerm + integralTerm)
			# axesList[1].legend(fontsize=FONT_SIZE, loc=4,frameon=False)
			axesList[1].set_ylabel("Correction terms", fontsize=FONT_SIZE)

			if firstPass:
				label1 = "expected"
				label2 = "simulated"

			axesList[2].plot(time / 60., rRnaSynthRate_expected, label = label1, color = "blue")
			axesList[2].plot(time / 60., rRnaSynthRate_updated, label = label2, color = "red")
			axesList[2].legend(fontsize=FONT_SIZE, loc=4,frameon=False)
			axesList[2].set_ylabel("rRNA synthesis prob", fontsize=FONT_SIZE)
			#axesList[2].set_ylim([0.16, 0.2])

			axesList[3].plot(time / 60., doublingTime.asNumber(units.min))
			axesList[3].set_ylabel("Inst. doubling\ntime (min)", fontsize=FONT_SIZE)
			axesList[3].plot(time / 60., expected_doubling_time, '--')
			axesList[3].set_ylim([20., 120.])

			axesList[4].plot(time / 60., effective_elongation_rate)
			axesList[4].set_ylabel("Eff. elng.\nrate", fontsize=FONT_SIZE)
			axesList[4].plot(time / 60., setpoint_elongation_rate, '--')
			axesList[4].plot(time / 60., expected_elongation_rate.asNumber(units.aa/units.s) * np.ones(time.size), '--')
			axesList[4].set_ylim([5., 22.])

			axesList[5].plot(time / 60., ribosomeConcentration.asNumber(units.mmol / units.L))
			axesList[5].set_ylabel("[Rib]", fontsize=FONT_SIZE)

			ax6.plot(rRnaSynthRate_expected, rRnaSynthRate_updated, '.')

			doublingTime = doublingTime.asNumber(units.min)
			doublingTime_nan = doublingTime[~np.isnan(doublingTime)]
			counts, bins = np.histogram(doublingTime_nan, bins = 100)
			binwidth=bins[1] - bins[0]
			counts_norm = counts / counts.max()
			ax7.bar(bins[:-1], counts_norm, width=binwidth, linewidth=0)

			firstPass = False

	ax6.plot([0.1,0.2],[0.1,0.2],linestyle='--',color="black")
	ax6.set_xlabel("Expected k_syn")
	ax6.set_ylabel("Simulation k_syn")
	ax6.set_ylim([0.1, 0.2])
	ax6.set_xlim([0.1, 0.2])

	ax7.set_xlim([20., 150.])
	ax7.set_xticks([24., 44., 100.])

	axesList[5].set_ylim([ribConcMin, ribConcMax])

	for a in axesList:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in a.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		# a.spines["right"].set_visible(False)
		a.spines["top"].set_visible(False)
		a.spines["bottom"].set_visible(False)
		a.xaxis.set_visible(False)

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
