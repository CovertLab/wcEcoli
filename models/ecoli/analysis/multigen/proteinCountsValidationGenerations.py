#!/usr/bin/env python
"""
Compare protein counts to Wisniewski 2014 data set

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/3/2015
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
import cPickle
from scipy.stats import pearsonr

import mpld3
from mpld3 import plugins, utils

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.containers.bulk_objects_container import BulkObjectsContainer

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.utils import units

from wholecell.utils.sparkline import whitePadSparklineAxis

# TODO: account for complexation
def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):


	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)


	sim_data = cPickle.load(open(simDataFile, "rb"))

	# import ipdb; ipdb.set_trace();
	validation_data = cPickle.load(open(validationDataFile, "rb"))

	ids_complexation = sim_data.process.complexation.moleculeNames
	ids_complexation_complexes = [ids_complexation[i] for i in np.where((sim_data.process.complexation.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_equilibrium_complexes = [ids_equilibrium[i] for i in np.where((sim_data.process.equilibrium.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_translation = sim_data.process.translation.monomerData["id"].tolist()
	ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))
	bulkContainer = BulkObjectsContainer(ids_protein, dtype = np.float64)
	view_complexation = bulkContainer.countsView(ids_complexation)
	view_complexation_complexes = bulkContainer.countsView(ids_complexation_complexes)
	view_equilibrium = bulkContainer.countsView(ids_equilibrium)
	view_equilibrium_complexes = bulkContainer.countsView(ids_equilibrium_complexes)
	view_translation = bulkContainer.countsView(ids_translation)
	view_validation = bulkContainer.countsView(validation_data.protein.wisniewski2014Data["monomerId"].tolist())
	view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())


	# Get all cells
	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)


	allDir = ap.get_cells()

	View_Validation_Schmidt = []


	# Steady-state analysis
	# data =  np.loadtxt(plotOutDir + "dataSupFig2-3Translation model-gene-depedentRIB-RNAaffinity.txt")
	# productionRate = data[:,0]
	# decayRate = data[:,1]

	# fig = plt.figure(figsize = (8, 8))

	# axis = plt.subplot(1,1,1)

	# axis.scatter(
	# 	productionRate,
	# 	decayRate,
	# 	alpha=.2,
	# 	s=60)
	# print pearsonr(productionRate, decayRate)[0]

	# maxLine = np.ceil( max(productionRate.max(), decayRate.max())) * 1.03

	# minLine = np.ceil( min(productionRate.min(), decayRate.min())) * 1.02

	# plt.plot([minLine, maxLine], [minLine, maxLine], '-k')

	# plt.xlim(xmin=minLine, xmax=maxLine)
	# plt.ylim(ymin=minLine, ymax=maxLine)


	# axis.spines["right"].set_visible(False)
	# axis.spines["top"].set_visible(False)
	# axis.spines["left"].set_position(("outward", 10))
	# axis.spines["bottom"].set_position(("outward", 10))
	# # axis.set_yticks(axis.get_ylim())
	# # axis.set_xticks(axis.get_xlim())
	# axis.tick_params(right = "off")
	# axis.tick_params(top = "off")
	# axis.tick_params(which = "both", direction = "out")

	# axis.set_xlim([minLine, maxLine])
	# axis.set_ylim([minLine, maxLine])

	# from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	# exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	# exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)
	# plt.close("all")


	# Simulated vs predicted protein coutns
	fig = plt.figure(figsize = (8, 8))

	for simDir in allDir:
		print simDir

		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_protein], np.int)
		proteinCountsBulk = bulkMolecules.readColumn("counts")[:, proteinIndexes]
		bulkMolecules.close()

		# Account for monomers
		bulkContainer.countsIs(proteinCountsBulk.mean(axis = 0))
		# bulkContainer.countsIs(proteinCountsBulk[-1,:]) # final time point

		# Account for unique molecules
		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		rnaPolyIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
		nActiveRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
		nActiveRnaPoly = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnaPolyIndex]
		uniqueMoleculeCounts.close()
		bulkContainer.countsInc(nActiveRibosome.mean(), sim_data.moleculeGroups.s30_fullComplex + sim_data.moleculeGroups.s50_fullComplex)
		bulkContainer.countsInc(nActiveRnaPoly.mean(), sim_data.moleculeGroups.rnapFull)

		# Account for small-molecule bound complexes
		view_equilibrium.countsInc(
			np.dot(sim_data.process.equilibrium.stoichMatrixMonomers(), view_equilibrium_complexes.counts() * -1)
			)

		# Account for monomers in complexed form
		view_complexation.countsInc(
			np.dot(sim_data.process.complexation.stoichMatrixMonomers(), view_complexation_complexes.counts() * -1)
			)

		wisniewskiCounts = validation_data.protein.wisniewski2014Data["avgCounts"]
		proteinIds = validation_data.protein.wisniewski2014Data["monomerId"].tolist()

		# Wisniewski Counts
		# points = ax[0].scatter(np.log10(wisniewskiCounts + 1), np.log10(view_validation.counts() + 1), c='w', edgecolor = 'k', alpha=.7)
		# ax[0].set_xlabel("Log 10(Wisniewski 2014 Counts)")
		# ax[0].set_title("Pearson r: %0.2f" % pearsonr(np.log10(view_validation.counts() + 1), np.log10(wisniewskiCounts + 1))[0])

		labels = list(proteinIds)
		# tooltip = plugins.PointLabelTooltip(points, labels)
		# plugins.connect(fig, tooltip)

		view_validation_schmidt = bulkContainer.countsView(validation_data.protein.schmidt2015Data["monomerId"].tolist())
		View_Validation_Schmidt.append(view_validation_schmidt.counts())
		

	View_Validation_Schmidt = (np.array(View_Validation_Schmidt)).mean(axis = 0)

	# Schmidt Counts
	schmidtLabels = validation_data.protein.schmidt2015Data["monomerId"]
	schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]


	axis = plt.subplot(1,1,1)

	axis.scatter(
		np.log10(schmidtCounts + 1),
		np.log10(View_Validation_Schmidt + 1),
		alpha=.2,
		s=60)
	print pearsonr( np.log10(View_Validation_Schmidt + 1), np.log10(schmidtCounts + 1) )[0]
	
	maxLine = np.ceil( 
					max((np.log10(schmidtCounts + 1)).max(), 
					(np.log10(View_Validation_Schmidt + 1)).max())
				)
	plt.plot([0, maxLine], [0, maxLine], '-k')

	plt.xlim(xmin=0, xmax=maxLine)
	plt.ylim(ymin=0, ymax=maxLine)


	axis.spines["right"].set_visible(False)
	axis.spines["top"].set_visible(False)
	axis.spines["left"].set_position(("outward", 10))
	axis.spines["bottom"].set_position(("outward", 10))
	# axis.set_yticks(axis.get_ylim())
	# axis.set_xticks(axis.get_xlim())
	axis.tick_params(right = "off")
	axis.tick_params(top = "off")
	axis.tick_params(which = "both", direction = "out")

	axis.set_xlim([-0.07, maxLine])
	axis.set_ylim([-0.07, maxLine])

	# NOTE: This Pearson correlation goes up (at the time of writing) about 0.05 if you only
	# include proteins that you have translational efficiencies for


	from wholecell.analysis.analysis_tools import exportFigure, exportHtmlFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	exportHtmlFigure(fig, plt, plotOutDir, plotOutFileName, metadata)
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
