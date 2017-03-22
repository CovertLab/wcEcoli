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
	fig.set_figwidth(8.5)
	fig.set_figheight(8.5/3)


	bremer_tau = [40, 100, 24]

	bremer_tau_low_err = [1.36, 6.03, 18.87]
	bremer_tau_high_err = [6.62, 8.63, 30.3]

	bremer_tau_low_err = [6.03, 18.87, 1.36]
	bremer_tau_high_err = [8.63, 30.3, 6.62]

	bremer_rna_mass_per_cell = [77, 20,  211]
	bremer_elng_rate = [18, 12,  21]
	#bremer_synth_ratio = [0.68, 0.41,  0.85]
	bremer_rrn_init_rate = [20*23, 4*12.4, 58*35.9]

	sim_doubling_time = np.zeros(ap.n_variant)
	sim_doubling_time_std = np.zeros(ap.n_variant)

	sim_rna_mass_per_cell = np.zeros(ap.n_variant)
	sim_elng_rate = np.zeros(ap.n_variant)
	sim_rrn_init_rate = np.zeros(ap.n_variant)

	sim_rna_mass_per_cell_std = np.zeros(ap.n_variant)
	sim_elng_rate_std = np.zeros(ap.n_variant)
	sim_rrn_init_rate_std = np.zeros(ap.n_variant)

	for varIdx in range(ap.n_variant):
		print "variant {}".format(varIdx)
		all_cells = ap.get_cells(variant=[varIdx], generation=[0])
		print "Total cells: {}".format(len(all_cells))
		import cPickle
		sim_data = cPickle.load(open(ap.get_variant_kb(varIdx)))

		mean_num_origin = np.zeros(len(all_cells))
		mean_num_terc = np.zeros(len(all_cells))
		mean_num_forks = np.zeros(len(all_cells))
		num_origin_at_init = np.zeros(len(all_cells))
		doubling_time = np.zeros(len(all_cells))

		for idx, simDir in enumerate(all_cells):
			print "cell {} of {}".format(idx, len(all_cells))

			simOutDir = os.path.join(simDir, "simOut")
			
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			doubling_time[idx] = time[-1] - time[0]
			timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

			rnaMass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("rnaMass")

			elngRate = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate")


			transcriptDataFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
			rnaSynth = transcriptDataFile.readColumn("countRnaSynthesized")
			isTRna = sim_data.process.transcription.rnaData["isTRna"]
			isRRna = sim_data.process.transcription.rnaData["isRRna"]
			isMRna = sim_data.process.transcription.rnaData["isMRna"]
			stableRnaSynth = np.sum(rnaSynth[:,isTRna], axis=1) + np.sum(rnaSynth[:,isRRna], axis=1)
			stableRnaSynth = rnaSynth[:,isTRna + isRRna]


			# rnaFraction = (stableRnaSynth.mean(axis=0) / timeStepSec.mean() * sim_data.process.transcription.rnaData["length"][isRRna + isTRna].asNumber()).sum() / 1e5
			rrn_init_rate = (rnaSynth[:, isRRna].sum(axis=1) / timeStepSec).mean() * 60. / 3
			# totalRnaSynth = np.sum(rnaSynth, axis=1).astype(float)
			# rnaFraction = stableRnaSynth / totalRnaSynth
			# mRnaSynth = np.sum(rnaSynth[:,isMRna], axis=1)
			# mRnaFraction = mRnaSynth / totalRnaSynth
			# import ipdb; ipdb.set_trace()

		sim_rna_mass_per_cell[varIdx] = rnaMass.mean()
		sim_elng_rate[varIdx] = elngRate.mean()
		sim_rrn_init_rate[varIdx] = np.nanmean(rrn_init_rate)
		sim_doubling_time[varIdx] = doubling_time.mean() / 60.

		sim_rna_mass_per_cell_std[varIdx] = mean_num_origin.std()
		sim_elng_rate_std[varIdx] = mean_num_terc.std()
		sim_rrn_init_rate_std[varIdx] = np.nanstd(rrn_init_rate)
		sim_doubling_time_std[varIdx] = doubling_time.std() / 60.
	
	ax0 = plt.subplot2grid((1,3), (0,0))
	ax1 = plt.subplot2grid((1,3), (0,1), sharex=ax0)
	ax2 = plt.subplot2grid((1,3), (0,2), sharex=ax0)

	lines = {'linestyle': 'None'}
	plt.rc('lines', **lines)

	ax0.errorbar(sim_doubling_time, sim_rna_mass_per_cell, yerr=sim_rna_mass_per_cell_std, xerr=sim_doubling_time_std, label="Simulation", color="black", fmt='', marker='o', markersize=3, linewidth=0.5)
	ax0.errorbar(bremer_tau, bremer_rna_mass_per_cell, yerr=np.array(bremer_rna_mass_per_cell) * 0.06, xerr = [bremer_tau_low_err, bremer_tau_high_err], label="Bremer & Dennis 1996", color="blue", marker='o', markersize=3, linewidth=0.5)
	ax0.set_title("RNA mass per cell (fg)", fontsize=FONT_SIZE)
	ax0.set_xlim([20, 105])

	ax1.errorbar(sim_doubling_time, sim_elng_rate, yerr=sim_elng_rate_std, xerr=sim_doubling_time_std, label="Simulation", color="black", fmt='', marker='o', markersize=3, linewidth=0.5)
	ax1.errorbar(bremer_tau, bremer_elng_rate, yerr=np.array(bremer_elng_rate) * 0.06,xerr = [bremer_tau_low_err, bremer_tau_high_err], label="Bremer & Dennis 1996", color="blue", marker='o', markersize=3, linewidth=0.5)
	ax1.set_title("Ribosome elongation\nrate (aa/s/ribosome)", fontsize=FONT_SIZE)
	ax1.set_xlabel("Doubling time (min)", fontsize=FONT_SIZE)

	ax2.errorbar(sim_doubling_time, sim_rrn_init_rate, yerr=sim_rrn_init_rate_std, xerr=sim_doubling_time_std, label="Simulation", color="black", fmt='', marker='o', markersize=3, linewidth=0.5)
	ax2.errorbar(bremer_tau, bremer_rrn_init_rate, yerr=np.array(bremer_rrn_init_rate) * 0.06, xerr = [bremer_tau_low_err, bremer_tau_high_err], label="Bremer & Dennis 1996", color="blue", marker='o', markersize=3, linewidth=0.5)
	ax2.set_title("Rate of rrn initiation (1/min)", fontsize=FONT_SIZE)

	ax2.legend(loc=1, frameon=True, fontsize=7)

	
	axes_list = [ax0, ax1, ax2]#, ax3]

	whitePadSparklineAxis(ax0)
	whitePadSparklineAxis(ax1)
	whitePadSparklineAxis(ax2)
	# whitePadSparklineAxis(ax3)

	for a in axes_list:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in a.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 

	plt.subplots_adjust(bottom = 0.2, wspace=0.3, top=0.88)

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
