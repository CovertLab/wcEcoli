#!/usr/bin/env python
"""
If metabolite changes over the course of generations, plot trend
Adapted from limitedMetabolites.py

@date: Created 5/17/2017
@author: Nicole Ferraro
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import math
import os
from re import findall

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.plotting_tools import COLORS_LARGE, COLORS_SMALL
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
import cPickle

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
	# Only get first daughter
	#allDir = ap.get_cells(daughter=[0])
	allDir = ap.get_cells()

	sim_data = cPickle.load(open(simDataFile, "rb"))

	try:
		gene_name = sim_data.genetic_perturbations.keys()[0]
	except:
		gene_name = "control"

	metaboliteNames = np.array(sorted(sim_data.process.metabolism.concDict.keys()))

	metaboliteFoldChange = np.array([])
	cellCycleLen = []
	change_names = []
	gens = []
	# change to range for length of allDir so not hard coded
	for gen,simDir in enumerate(allDir):

		matches = findall('generation_\d{6}', simDir)
		gens.append(int(matches[0][-6:]))

		simOutDir = os.path.join(simDir, "simOut")

		enzymeKineticsData = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		metaboliteCounts = enzymeKineticsData.readColumn("metaboliteCountsFinal")
		normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		cellCycleLen.append(time[-1] - time[0])
		total_changes = normalizedCounts[-1,:] - normalizedCounts[0,:]
		total_changes[np.isnan(total_changes)] = 0
		enzymeKineticsData.close()

		if gen == 0:
			change_inds = np.where(total_changes < 1.9)[0]
			change_names = metaboliteNames[change_inds]
			metaboliteFoldChange = total_changes[change_inds]
		
		elif len(change_names) > 0:
			change_inds = []
			for name in change_names:
				change_inds.append(np.where(metaboliteNames == name)[0][0])
			metaboliteFoldChange = np.vstack((metaboliteFoldChange, total_changes[change_inds]))

	# Only plot if any metabolites were impacted (fold change < 1.9 over a cell cycle for any generation)

	if metaboliteFoldChange.shape[0] > 0:
		colors = COLORS_SMALL 
		fig, axesList = plt.subplots(2)
		fig.set_size_inches(11, 11)
		ax0 = axesList[0]
		#ax = plt.subplot(1, 1, 1)
		ax0.set_color_cycle(colors)
		max_y = math.ceil(np.max(metaboliteFoldChange))
		ax0.set_ylim([0,max_y])
		num_mets = metaboliteFoldChange.shape[1]
		gen_array = np.array([gens]*num_mets).transpose()

		for i in range(0,metaboliteFoldChange.shape[1]):
			ax0.scatter(gen_array[:,i], metaboliteFoldChange[:,i], color=COLORS_SMALL[i])
		
		ax0.set_xlabel("Generations")
		ax0.set_ylabel("Metabolite Fold Change over Cell Cycle")
		
		ax0.legend(change_names, fontsize=12)
		ax0.set_title(gene_name)

		ax1 = axesList[1]
		ax1.scatter(gens, cellCycleLen)
		ax1.set_color_cycle(colors)
		ax1.set_xlabel('Generations')
		ax1.set_ylabel('Cell Cycle Len (s)')
		ax1.set_title(gene_name + ' cell cycle length')
		
		from wholecell.analysis.analysis_tools import exportFigure
		exportFigure(fig, plotOutDir, plotOutFileName,metadata)
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
