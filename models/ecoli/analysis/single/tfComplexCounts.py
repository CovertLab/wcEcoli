#!/usr/bin/env python
"""
Plots TF complexes to see if they exist

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/7/15
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

IGNORE_FIRST_PERCENTAGE = 0.1

from math import log10, floor
def round_to_1(x):
	if x < 0:
		x = x*-1
	return -1*round(x, -int(floor(log10(x))))

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	complexIds = [x["id"] for x in kb.moleculeGroups.tfComplexCounts]
	complexNames = [x["name"] for x in kb.moleculeGroups.tfComplexCounts]

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Calculate concentration data
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")

	complexIdxs = np.array([bulkMoleculeIds.index(x) for x in complexIds])
	complexCounts = bulkMolecules.readColumn("counts")[:, complexIdxs]

	bulkMolecules.close()
	
	fig = plt.figure(figsize = (11, 11))
	rows = 13
	cols = 12
	idx = 0
	for idx in range(len(complexIds)):
		ax = plt.subplot(rows, cols, idx+1)

		# deviation = 1-poolConc[:,idx]/concSetpoint[:,idx]
		counts = complexCounts[1:, idx]
		ax.plot(time[1:] / 60., counts, linewidth=1, label="counts", color='k')

		# Highlights low counts
		bbox = None
		if np.mean(counts) < 10:
			bbox = {'facecolor':'red', 'alpha':0.5, 'pad':1}
		ax.set_title('{} - {}\navg: {}'.format(complexIds[idx][:-3], complexNames[idx], int(np.mean(counts))), fontsize=6, bbox=bbox)

		# Sets ticks so that they look pretty
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=6)
		ax.set_xticks([])
		ymin = counts[counts.shape[0] * IGNORE_FIRST_PERCENTAGE:].min()
		ymax = counts[counts.shape[0] * IGNORE_FIRST_PERCENTAGE:].max()
		ax.set_ylim([ymin, ymax])
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels([str(counts.min()), str(counts.max())])

	# Create legend
	ax = plt.subplot(rows, cols,len(complexIds) + 2)
	ax.plot(0, 0, linewidth=2, label="count", color='k')
	ax.legend(loc = 10,prop={'size':10})
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_title("Highlights low counts", fontsize=12, bbox={'facecolor':'red', 'alpha':0.5, 'pad':1})

	# Save
	plt.subplots_adjust(hspace = 1, wspace = 1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

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
