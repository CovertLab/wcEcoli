"""
Compare tRNA supply vs AA demand for each timestep

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/4/18
"""

from __future__ import absolute_import


import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle
import bokeh.io
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import (HoverTool, BoxZoomTool, LassoSelectTool, PanTool,
	WheelZoomTool, ResizeTool, UndoTool, RedoTool)

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot


UPSHIFT_VARIANT_INDEX = 2
UPSHIFT_LABEL = "000002_add_aa"
PRESHIFT_GENS = [0, 1, 2, 3]
POSTSHIFT_GENS = [5, 6, 7, 8]

def getAminoAcidDemand(simDirs):
	"""
	Computes average amino acid demand over the timesteps.
	Note: removes entry for selenocysteine.
	"""
	out = []
	for simDir in simDirs:
		simOutDir = os.path.join(simDir, "simOut")

		# Get amino acid usage
		growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))
		aasUsed = growthLimits.readColumn("aasUsed")
		growthLimits.close()

		# Divide each row by it's sum
		aasUsed_normalized = aasUsed / aasUsed.sum(axis = 1, keepdims = True)

		# Save data
		if not len(out):
			out = aasUsed_normalized
		else:
			out = np.vstack((out, aasUsed_normalized))

	# Remove selenocysteine's entry
	mask = np.ones(out.shape[1], dtype = bool)
	mask[-2] = False
	return out[:, mask]

def plotHist(raw_data, ax, color):
	mask = ~np.isnan(raw_data)
	data = raw_data[mask]
	bins = np.arange(data.min(), data.max(), 0.001)
	ax.hist(data, color = color, alpha = 0.5, bins = bins, log = True)
	return

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Check that the sim is a nutrient upshift
		ap = AnalysisPaths(inputDir, variant_plot = True)
		if UPSHIFT_VARIANT_INDEX not in ap.get_variants():
			print "This plot only runs for the nutrient timeseries variant at index: %s" % UPSHIFT_VARIANT_INDEX
			return

		sim_data = cPickle.load(open(ap.get_variant_kb(UPSHIFT_VARIANT_INDEX), "rb"))
		if not hasattr(sim_data, "nutrientsTimeSeriesLabel"):
			print "This plot only runs for a nutrient timeseries variant."
			return

		if sim_data.nutrientsTimeSeriesLabel != UPSHIFT_LABEL:
			print "This plot only runs for the nutrient timeseries labelled: %s" % UPSHIFT_LABEL
			return

		# Get cells
		preshift_cells = ap.get_cells(variant = [UPSHIFT_VARIANT_INDEX], generation = PRESHIFT_GENS)
		postshift_cells = ap.get_cells(variant = [UPSHIFT_VARIANT_INDEX], generation = POSTSHIFT_GENS)

		# Get amino acid usage
		preshift_demand = getAminoAcidDemand(preshift_cells)
		postshift_demand = getAminoAcidDemand(postshift_cells)

		# Plot
		fig, axesList = plt.subplots(5, 4, figsize = (8.5, 11))
		axesList = axesList.flatten()
		aaIds = sim_data.moleculeGroups.aaIDs
		mask = np.ones(len(aaIds), dtype = bool)
		mask[-2] = False
		aaIds = np.array(aaIds)[mask]

		for i, ax in enumerate(axesList):
			plotHist(preshift_demand[:, i], ax, "b")
			plotHist(postshift_demand[:, i], ax, "r")
			ax.set_xticks(ax.get_xlim())
			ax.set_title(aaIds[i], fontsize = 10)

			# Report mean average
			ax.text(0.5, 1.2, "%0.3f" % np.nanmean(preshift_demand[:, i]), color = "b", transform = ax.transAxes, ha = "center", fontsize = 10)
			ax.text(0.5, 1.35, "%0.3f" % np.nanmean(postshift_demand[:, i]), color = "r", transform = ax.transAxes, ha = "center", fontsize = 10)

		plt.xlabel("blue = preshift\nred = postshift")
		plt.subplots_adjust(wspace = 0.5, hspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
