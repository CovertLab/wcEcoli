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
	std = np.nanstd(a)
	a[a < median - 1.5*std] = np.nan
	a[a > median + 1.5*std] = np.nan

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	
	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	fig = plt.figure()
	fig.set_figwidth(8)
	fig.set_figheight(4)

	ax0 = plt.subplot2grid((1,2), (0,0))
	ax1 = plt.subplot2grid((1,2), (0,1))

	expected_doubling_time = np.zeros(3)
	expected_growth_rate = np.zeros(3)

	all_doubling_time = np.zeros(3)
	all_doubling_time_std = np.zeros(3)
	all_growth_rates = np.zeros(3)
	all_growth_rates_std = np.zeros(3)

	for varIdx in range(ap.n_variant):

		all_cells = ap.get_cells(variant=[varIdx])
		print len(all_cells)

		doubling_times = np.zeros(len(all_cells))
		growth_rates = np.zeros(len(all_cells))

		import cPickle
		simDataFile = ap.get_variant_kb(varIdx)
		sim_data = cPickle.load(open(simDataFile, "rb"))

		for idx, simDir in enumerate(all_cells):
			simOutDir = os.path.join(simDir, "simOut")

			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			doubling_times[idx] = (time[-1] - time[0]) / 60.

			# Instantanious growth rate
			growthRate = TableReader(os.path.join(simOutDir, "Mass")).readColumn("instantaniousGrowthRate")
			growthRate = growthRate * 60.
			growth_rates[idx] = np.nanmean(growthRate)

		all_doubling_time[varIdx] = doubling_times.mean()
		all_doubling_time_std[varIdx] = doubling_times.std()
		expected_doubling_time[varIdx] = sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min)

		all_growth_rates[varIdx] = np.nanmean(growth_rates)
		all_growth_rates_std[varIdx] = np.nanstd(growth_rates)
		expected_growth_rate[varIdx] = np.log(2) / sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min)
		ax0.plot(sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min) * np.ones(doubling_times.size), doubling_times, '.', color = "grey", alpha = 0.5, zorder=1, markersize=10)
		ax1.plot(np.log(2) / sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.min) * np.ones(growth_rates.size), growth_rates, '.', color = "grey", alpha = 0.5, zorder=1, markersize=10)
		ax0.plot([0,200],[0,200], color="black", linewidth=0.5, linestyle='--')
		ax1.plot([0,0.04],[0,0.04], color="black", linewidth=0.5, linestyle='--')



	lines = {'linestyle': 'None'}
	plt.rc('lines', **lines)

	ax0.errorbar(expected_doubling_time, all_doubling_time, yerr=all_doubling_time_std, fmt='', marker='o', markersize=3, linewidth=1, color="blue")
	ax1.errorbar(expected_growth_rate, all_growth_rates, yerr=all_growth_rates_std, fmt='', marker='o', markersize=3, linewidth=1, color="blue")

	print all_growth_rates_std / all_growth_rates



	ax0.set_xlim([0,200])
	ax0.set_ylim([0,200])

	ax1.set_xlim([0,0.04])
	ax1.set_ylim([0,0.04])


	# 	ax0.set_xlim([0., time.max() / 60.])
	# 	if varIdx == 0:
	# 		ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
	# 	if varIdx > 0:
	# 		ax0.yaxis.set_visible(False)
	# 	ax0.set_ylim([400, 6000.])
	# 	ax0.set_title(title_list[varIdx], fontsize=FONT_SIZE)
	# 	# ax0.xaxis.set_visible(False)
	# 	#ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	# 	# nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
	# 	# try:
	# 	# 	T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0] / 60.
	# 	# except IndexError:
	# 	# 	print "No shift detected"
	# 	# 	T_ADD_AA = 1e9

	# 	axes_list = [ax0, ax1, ax2, ax3, ax4, ax5]
	# 	# for a in axes_list:
	# 	# 	shift_time = T_ADD_AA
	# 	# 	width = a.get_xlim()[1] - shift_time
	# 	# 	height = a.get_ylim()[1] - a.get_ylim()[0]
	# 	# 	a.add_patch(
	# 	# 		patches.Rectangle(
	# 	# 			(shift_time, a.get_ylim()[0]),   # (x,y)
	# 	# 			width,          # width
	# 	# 			height,          # height
	# 	# 			alpha = 0.25,
	# 	# 			color = "gray",
	# 	# 			linewidth = 0.
	# 	# 		)
	# 	# 	)

	# 	for a in axes_list:
	# 		for tick in a.yaxis.get_major_ticks():
	# 			tick.label.set_fontsize(FONT_SIZE) 
	# 		for tick in a.xaxis.get_major_ticks():
	# 			tick.label.set_fontsize(FONT_SIZE) 

	# 	if varIdx == 0:
	# 		ax1.set_ylabel("Instantanious\ngrowth rate\n(fg/fg-min)", fontsize=FONT_SIZE)
	# 	if varIdx > 0:
	# 		ax1.yaxis.set_visible(False)
	# 	ax1.xaxis.set_visible(False)
	# 	ax1.set_ylim([0.003, 0.035])

	# 	if varIdx == 0:
	# 		ax2.set_ylabel("70S concentration\n(" + r"$\mu$" + "mol/L)", fontsize=FONT_SIZE)
	# 	if varIdx > 0:
	# 		ax2.yaxis.set_visible(False)
	# 	ax2.set_ylim([16, 30])
	# 	ax2.xaxis.set_visible(False)

	# 	if varIdx == 0:
	# 		ax3.set_ylabel("Average 70S\nelongation rate (aa/s)", fontsize=FONT_SIZE)
	# 	if varIdx > 0:
	# 		ax3.yaxis.set_visible(False)
	# 	ax3.xaxis.set_visible(False)
	# 	ax3.set_ylim([8,22])

	# 	# y_ticks = ax2.get_yticks()
	# 	# new_y_ticks = y_ticks[0:-1:2]
	# 	# ax2.set_yticks(new_y_ticks)

	# 	if varIdx == 0:
	# 		ax4.set_ylabel("Metabolic supply\n(count/s-fg)", fontsize=FONT_SIZE)
	# 	if varIdx > 0:
	# 		ax4.yaxis.set_visible(False)
	# 	ax4.xaxis.set_visible(False)
	# 	ax4.set_ylim([100, 400])


	# 	if varIdx == 0:
	# 		ax5.set_ylabel("rrn prodction\nrate (1/s-fg)", fontsize=FONT_SIZE)
	# 	if varIdx > 0:
	# 		ax5.yaxis.set_visible(False)
	# 	ax5.set_ylim([0.001, 0.01])


	# 	if varIdx == 1:
	# 		ax5.set_xlabel("Time (min)")

	whitePadSparklineAxis(ax0)
	whitePadSparklineAxis(ax1)

	ax0.set_xlabel("Expected", fontsize=FONT_SIZE)
	ax1.set_xlabel("Expected", fontsize=FONT_SIZE)
	ax0.set_ylabel("Simulated", fontsize=FONT_SIZE)

	ax0.set_title("Expected vs simulated doubling time (min)\nn={} cells per condition".format(len(all_cells)), fontsize=FONT_SIZE)
	ax1.set_title("Expected vs simulated growth rate (1/min)\nn={} cells per condition".format(len(all_cells)), fontsize=FONT_SIZE)

	# 	whitePadSparklineAxis(ax2, False)
	# 	whitePadSparklineAxis(ax3, False)
	# 	whitePadSparklineAxis(ax4, False)
	# 	whitePadSparklineAxis(ax5)

	plt.subplots_adjust(bottom =0.2)

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
