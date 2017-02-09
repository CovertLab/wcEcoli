#!/usr/bin/env python
"""
Plot RNAse counts

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/14/2015
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import scipy.fftpack
import cPickle

import wholecell.utils.constants

from wholecell.io.tablereader import TableReader

from wholecell.utils.sparkline import whitePadSparklineAxis


def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")

	sim_data = cPickle.load(open(simDataFile, "rb"))

	endoRnaseIds = sim_data.process.rna_decay.endoRnaseIds
	exoRnaseIds = sim_data.moleculeGroups.exoRnaseIds
	RNase_IDS = np.concatenate((endoRnaseIds, exoRnaseIds))

	endoRnase_RnaIDs = sim_data.moleculeGroups.endoRnase_RnaIDs
	exoRnase_RnaIDs = sim_data.moleculeGroups.exoRnase_RnaIDs
	RNase_RnaIDS = np.concatenate((endoRnase_RnaIDs, exoRnase_RnaIDs))
	RNase_IDS = np.concatenate((RNase_IDS, RNase_RnaIDS))

	
	# RNase P and PNP
	RNase_IDS = ['EG10862_RNA[c]', 'EG10743_RNA[c]', 'EG10862-MONOMER[c]', 'EG10743-MONOMER[c]']


	rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNase_IDS], np.int)
	rnapRnaCounts = bulkMolecules.readColumn("counts")[:, rnapRnaIndexes]
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	bulkMolecules.close()

	# plt.figure(figsize = (8.5, 11))
	# plt.rc('xtick', labelsize=7) 
	# plt.rc('ytick', labelsize=5)

	plt.figure(figsize = (8.5, 5))
	plt.rc('xtick', labelsize=12) 
	plt.rc('ytick', labelsize=12)

	count = 0
	count_bis = len(RNase_IDS) / 2

	freq_observed = []

	# import ipdb; ipdb.set_trace()

	for subplotIdx in xrange(0, len(RNase_IDS)):
		if not subplotIdx % 2:
			rnapRnaCountsIdx = count
			count += 1
		if subplotIdx % 2:
			rnapRnaCountsIdx = count_bis
			count_bis += 1

		ax = plt.subplot(len(RNase_IDS) / 2, 2, 1 + subplotIdx)

		whitePadSparklineAxis(ax)

		plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])

		if not subplotIdx >= len(RNase_IDS) - 2:
			frame = plt.gca()
			for xlabel_i in frame.axes.get_xticklines():
				xlabel_i.set_visible(True)
			for xlabel_i in frame.axes.get_xticklabels():
				xlabel_i.set_visible(False)

		# if subplotIdx >= len(RNase_IDS) - 2:
		# 	plt.xlabel("Time (min)", fontsize = 7)

		# if not subplotIdx % 2:	
		# 	plt.ylabel("RNA counts", fontsize = 10)
		# if subplotIdx % 2:
		# 	plt.ylabel("Protein counts", fontsize = 10)

		# plt.title(RNase_IDS[rnapRnaCountsIdx], fontsize = 10) 

		max_yticks = 4
		yloc = plt.MaxNLocator(max_yticks)
		ax.yaxis.set_major_locator(yloc)

		signal = rnapRnaCounts[:, subplotIdx]
		if subplotIdx == 17:
			np.savetxt(os.path.join(plotOutDir, 'PNPase-MONOMER[c].txt'), signal)
		if subplotIdx == 35:
			np.savetxt(os.path.join(plotOutDir, 'PNPase-RNA[c].txt'), signal)

		# identifying periodicity on RNA and protein copy numbers
		fourier = np.fft.fft(signal)

		# copmuting mean and std values of the power spectral density
		M = np.mean(abs(fourier))
		S = np.std(abs(fourier))

		# computing frequencies 
		n = signal.size
		timestep = 1 # second
		freq = np.fft.fftfreq(n, d = timestep)
		fft_freq = sorted(zip(abs(fourier),freq))[n-6:n]

		# identifing peaks (frequency and period) with maximum PSD
		for i in xrange(0,len(fft_freq)):
			if fft_freq[i][1] > 0.: # only positive frequencies
				if 1. / fft_freq[i][1] < 3600.: # only periods lower than the doubling time
					if abs(fft_freq[i][0] - M) / S > 3: # strong and significant fft
						if subplotIdx == 17: # PNPase monomer
							# print RNase_IDS[subplotIdx], 1. / fft_freq[i][1] / 60. # period (min)
							freq_observed.append([1, 1. / fft_freq[i][1] / 60.])

						if subplotIdx == 35: # PNPase rna
							# print RNase_IDS[subplotIdx], 1. / fft_freq[i][1] / 60. # period (min)
							freq_observed.append([0, 1. / fft_freq[i][1] / 60.])


	index_generation = plotOutDir[plotOutDir.index('generation_') + 11 : plotOutDir.index('generation_') + 17]
	fileFrequencies = 'PNPase_frequencies_' + index_generation + '.txt'
	np.savetxt(os.path.join(plotOutDir, fileFrequencies), freq_observed)


	plt.subplots_adjust(hspace = 0.75, top = 0.95, bottom = 0.05)
	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
