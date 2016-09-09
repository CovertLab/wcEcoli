#!/usr/bin/env python
"""
Plots fraction of mRNAs transcribed (out of all genes to be transcribed) for all seeds.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/8/2016
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get IDs of mRNAs
	sim_data = cPickle.load(open(simDataFile, "rb"))
	rnaIds = sim_data.process.transcription.rnaData["id"]
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	mRnaIds = np.where(isMRna)[0]
	mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

	synthProb = sim_data.process.transcription.rnaSynthProb["basal"]
	mRnaSynthProb = np.array([synthProb[x] for x in mRnaIds])

	expression = sim_data.process.transcription.rnaExpression["basal"]
	mRnaExpression = np.array([expression[x] for x in mRnaIds])

	degRate = sim_data.process.transcription.rnaData["degRate"]
	mRnaDegRate = np.array([degRate[x].asNumber(1 / units.s) for x in mRnaIds])


	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells()
	numMRnas = mRnaNames.shape[0]
	numCells = all_cells.shape[0]
	numerical_zero = 1e-03

	# Get number of mRNA transcripts produced
	transcribedCounts = []
	for n, simDir in enumerate(all_cells):
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mRnaIndexes = np.array([moleculeIds.index(x) for x in mRnaNames])
		mRnaCounts = bulkMolecules.readColumn("counts")[:, mRnaIndexes]
		bulkMolecules.close()

		rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
	 	countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
	 	countMRnaDegraded = countRnaDegraded[:, mRnaIds]

		
		# Calculate number of transcripts produced
		mRnaDegraded_skip0 = countMRnaDegraded[1:, :] # remove first timestep
		mRnaCounts_skip0 = mRnaCounts[1:, :]
		mRnaCounts_skiplast = mRnaCounts[:-1, :]

		mRnaProducedCounts = mRnaCounts_skip0 - (mRnaCounts_skiplast - mRnaDegraded_skip0)
		mRnaProducedCountsSumOverTime = mRnaProducedCounts.sum(axis = 0)

		transcribedCounts.append(mRnaProducedCountsSumOverTime)

	transcribedCounts = np.array(transcribedCounts)
	transcribedCountsSumOverCells = transcribedCounts.sum(axis = 0)
	transcribedFreq = transcribedCountsSumOverCells / float(numCells)


	# Plot
	fig = plt.figure(figsize = (14, 10))
	ax = plt.subplot(1, 1, 1)

	zerosIndex = np.where(transcribedFreq == 0)
	log10_transcribedFreq = np.log10(transcribedFreq)
	log10_transcribedFreq[zerosIndex] = np.log10(numerical_zero)

	ax.scatter(np.log10(mRnaSynthProb), log10_transcribedFreq, facecolors = "none", edgecolors = "b")
	ax.set_title("Correlation of synthesis probability and number of transcripts produced\nn = %s cells" % numCells, fontsize = 12)
	ax.set_xlabel("log_10(synthesis probability)")
	ax.set_ylabel("log_10(Average number of transcripts produced per generation)")
	ax.tick_params(which = "both", direction = "out", top = "off")
	ax.spines["top"].set_visible(False)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


	# from bokeh.plotting import figure, output_file, show, ColumnDataSource, save
	# from bokeh.models import HoverTool

	# import ipdb; ipdb.set_trace()
	# # Output to html file
	# output_file(str(plotOutFileName) + ".html")

	# # Hover
	# source = ColumnDataSource(data = dict(
	# 	x = np.log10(mRnaSynthProb),
	# 	y = np.log10(transcribedFreq),
	# 	mRnaId = mRnaNames,
	# 	synthProb = mRnaSynthProb,
	# 	degRate = mRnaDegRate,
	# 	expression = mRnaExpression,
	# 	)
	# )

	# hover = HoverTool(
	# 	tooltips = [
	# 		("mRna ID", "@mRnaId"),
	# 		("synthesis probability", "@synthProb"),
	# 		("degradation rate", "@degRate"),
	# 		("expression", "@expression"),
	# 		]
	# 	)

	# # Create a plot with title and axis labels
	# plot = figure(
	# 	title = "Correlation of synthesis probability and number of transcripts produced\nn = %s cells" % numCells,
	# 	title_text_font_size = ["8pt"], 
	# 	x_axis_label = "log_10(synthesis probability)", 
	# 	y_axis_label = "log_10(Average number of transcripts produced per generation)",
	# 	width = 800,
	# 	height = 500,
	# 	tools = [hover, "reset", "wheel_zoom", "box_zoom", "pan", "resize", "lasso_select", "tap", "save"],
	# 	)

	# plot.scatter("x", "y", source = source,
	# 	size = 5,	 
	# 	fill_color = "navy", 
	# 	fill_alpha = 0.5, 
	# 	line_color = None, )


	# plot.xaxis.axis_label_text_font_size = ["6pt"]
	# plot.yaxis.axis_label_text_font_size = ["6pt"]

	# save(plot)



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
