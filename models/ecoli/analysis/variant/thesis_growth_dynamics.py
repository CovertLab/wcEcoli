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
WIDTH = 100

# def sparklineAxis(axis):
# 	axis.spines['top'].set_visible(False)
# 	axis.spines['bottom'].set_visible(False)
# 	axis.xaxis.set_ticks_position('none')
# 	axis.tick_params(which = 'both', direction = 'out')


def mm2inch(value):
	return value * 0.0393701

def remove_artifacts(a):
	median = np.median(a)
	std = np.std(a)
	a[a < median - 1.5*std] = np.nan
	a[a > median + 1.5*std] = np.nan

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	
	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	fig = plt.figure()
	fig.set_figwidth(8.5)
	fig.set_figheight(10)

	title_list = ["Glucose minimal anaerobic\n" + r"$\mu = $" + "100 min", "Glucose minimal\n" + r"$\mu = $" + "44 min", "Glucose minimal + 20 amino acids\n" + r"$\mu = $" + "22 min"]

	for varIdxReal in range(ap.n_variant):

		all_cells = ap.get_cells(variant=[varIdxReal], seed=[0], generation=[0,1,2,3,4])

		import cPickle
		simDataFile = ap.get_variant_kb(all_cells[0])
		sim_data = cPickle.load(open(simDataFile, "rb"))

		oriC = sim_data.constants.oriCCenter.asNumber()
		terC = sim_data.constants.terCCenter.asNumber()
		genomeLength = len(sim_data.process.replication.genome_sequence)


		if varIdxReal == 0:
			varIdx = 1
		elif varIdxReal == 1:
			varIdx = 0
		elif varIdxReal == 2:
			varIdx = 2
		else:
			raise Exception()

		ax0 = plt.subplot2grid((7,ap.n_variant), (0,varIdx))
		ax1 = plt.subplot2grid((7,ap.n_variant), (1,varIdx), sharex=ax0)
		ax2 = plt.subplot2grid((7,ap.n_variant), (2,varIdx), sharex=ax0)
		ax3 = plt.subplot2grid((7,ap.n_variant), (3,varIdx), sharex=ax0)
		ax4 = plt.subplot2grid((7,ap.n_variant), (4,varIdx), sharex=ax0)
		ax5 = plt.subplot2grid((7,ap.n_variant), (5,varIdx), sharex=ax0)
		ax6 = plt.subplot2grid((7,ap.n_variant), (6,varIdx), sharex=ax0)


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
			ax0.plot(time / 60., cellMass, color = color, alpha = alpha, linewidth=1)

			# Get active ribosome counts
			uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
			ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
			uniqueMoleculeCounts.close()
			ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass))
			ribosomeConcentration = ribosomeConcentration.asNumber(units.umol / units.L)
			remove_artifacts(ribosomeConcentration)
			ax1.plot(time / 60., ribosomeConcentration, color = color, alpha = alpha, linewidth=1)

			# Instantanious growth rate
			growthRate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
			remove_artifacts(growthRate)
			ax2.plot(time / 60., growthRate * 60., color = color, alpha = alpha, linewidth=1)

			# Get ribosome elongation rate and moving average
			elongationRate = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate")
			elongationRate_smooth = np.convolve(elongationRate, np.ones(WIDTH) / WIDTH, mode = "same")
			remove_artifacts(elongationRate_smooth)
			ax3.plot(time / 60., elongationRate_smooth, color = color, alpha = alpha, linewidth=1)

			# Supplied vs used capacity
			timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
			processTranslationUsage = elongationRate * ribosomeCounts / cellMass / timeStepSec
			remove_artifacts(processTranslationUsage)
			ax4.plot(time / 60., processTranslationUsage, color = color, alpha = alpha, linewidth=1)

			# Synthesis prob
			rRnaSynProb = TableReader(os.path.join(simOutDir, "ControlLoop")).readColumn("rRnaSynthRate_updated")
			ax5.plot(time / 60., rRnaSynProb, color = color, alpha = alpha, linewidth=1)

			# Initiation rates
			rrn16S_produced = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn16S_produced")
			rrn23S_produced = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn23S_produced")
			rrn5S_produced = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("rrn5S_produced")

			rrn16s_init_rate = rrn16S_produced / timeStepSec / cellMass
			rrn23s_init_rate = rrn23S_produced / timeStepSec / cellMass
			rrn5s_init_rate = rrn5S_produced / timeStepSec / cellMass

			net_rrn = rrn16s_init_rate + rrn23s_init_rate + rrn5s_init_rate

			net_rrn_smooth = np.convolve(net_rrn, np.ones(WIDTH) / WIDTH, mode = "same")

			remove_artifacts(net_rrn_smooth)

			ax6.plot(time / 60., net_rrn_smooth, color = color, alpha = alpha, linewidth=1)
 

		ax0.set_xlim([0., time.max() / 60.])
		if varIdx == 0:
			ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax0.yaxis.set_visible(False)
		ax0.set_ylim([400, 6000.])
		ax0.set_title(title_list[varIdx], fontsize=FONT_SIZE)
		# ax0.xaxis.set_visible(False)
		#ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

		# nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
		# try:
		# 	T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0] / 60.
		# except IndexError:
		# 	print "No shift detected"
		# 	T_ADD_AA = 1e9

		axes_list = [ax0, ax1, ax2, ax3, ax4]
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
			ax1.set_ylabel("70S concentration\n(" + r"$\mu$" + "mol/L)", fontsize=FONT_SIZE)
		# if varIdx > 0:
		# 	ax1.yaxis.set_visible(False)
		ax1.set_ylim([16, 30])
		ax1.xaxis.set_visible(False)


		if varIdx == 0:
			ax2.set_ylabel("Instantanious\ngrowth rate\n(fg/fg-min)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax2.yaxis.set_visible(False)
		ax2.xaxis.set_visible(False)
		ax2.set_ylim([0.005, 0.03])

		if varIdx == 0:
			ax3.set_ylabel("Average 70S\nelongation rate", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax3.yaxis.set_visible(False)
		ax3.xaxis.set_visible(False)
		ax3.set_ylim([8,21])

		# y_ticks = ax2.get_yticks()
		# new_y_ticks = y_ticks[0:-1:2]
		# ax2.set_yticks(new_y_ticks)

		if varIdx == 0:
			ax4.set_ylabel("Metabolic supply\n(count/s)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax4.yaxis.set_visible(False)
		ax4.xaxis.set_visible(False)
		ax4.set_ylim([100, 400])


		if varIdx == 0:
			ax5.set_ylabel("Synth prob", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax5.yaxis.set_visible(False)
		ax5.xaxis.set_visible(False)
		ax5.set_ylim([0.12, 0.2])

		if varIdx == 0:
			ax6.set_ylabel("rrn prodction\nrate (1/s-fg)", fontsize=FONT_SIZE)
		if varIdx > 0:
			ax6.yaxis.set_visible(False)
		ax6.set_ylim([0.003, 0.03])


		whitePadSparklineAxis(ax0, False)
		whitePadSparklineAxis(ax1, False)
		whitePadSparklineAxis(ax2, False)
		whitePadSparklineAxis(ax3, False)
		whitePadSparklineAxis(ax4, False)
		whitePadSparklineAxis(ax5, False)
		whitePadSparklineAxis(ax6)

	plt.subplots_adjust(left = 0.18, wspace=0.3, right=0.92)

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
