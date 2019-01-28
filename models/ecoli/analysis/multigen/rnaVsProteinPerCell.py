"""
Plots average RNA counts per cell vs average protein counts per cell.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/1/2017
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from models.ecoli.analysis import multigenAnalysisPlot


PLOT_ZEROS_ON_LINE = 2.5e-6
SKIP_FIRST_N_TIME_STEPS = 5 # complexes do not form yet during these time steps
COLOR_FREQ_0 = "y"
COLOR_FREQ_1 = "r"
COLOR_FREQ_SUBGEN = "b"

def replace_zeros(array):
	mask = array == 0
	array[mask] = PLOT_ZEROS_ON_LINE
	return array

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Load sim data
		sim_data = cPickle.load(open(simDataFile, "rb"))

		# Load rna IDs, ordered by monomer IDs
		rnaIds = sim_data.process.transcription.rnaData["id"][sim_data.relation.rnaIndexToMonomerMapping]

		# Load complex IDs
		complexIds = sim_data.process.complexation.ids_complexes

		# Get data from all generations to:
		# 1) identify sub-generationally transcribed genes
		# 2) get average RNA per cell
		# 3) get average monomers (sum from de-complexing) per cell
		# 4) get average complexes per cell
		rnaSumIsNonZero = []
		avgRnaCounts = []
		avgMonomerCountsFromDeComplexing = []
		avgComplexCounts = []
		for i, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			# Read from bulk molecules
			rnaCounts, complexCounts = read_bulk_molecule_counts(simOutDir, (rnaIds, complexIds))

			# Sum RNA counts over timesteps
			rnaCountsSum = rnaCounts.sum(axis=0)

			# Flag where the sum is nonzero (True if nonzero, False if zero)
			rnaSumIsNonZero.append(rnaCountsSum != 0)

			# Record average RNA per cell
			avgRnaCounts.append(rnaCounts.mean(axis=0))

			# Record average monomers (sum from de-complexing) per cell
			monomerCountsReader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			monomerCounts = monomerCountsReader.readColumn("monomerCounts")
			if i == 0:
				avgMonomerCountsFromDeComplexing.append(monomerCounts[SKIP_FIRST_N_TIME_STEPS:, :].mean(axis=0))
			else:
				avgMonomerCountsFromDeComplexing.append(monomerCounts.mean(axis=0))

			# Record average complexes per cell
			avgComplexCounts.append(complexCounts.mean(axis=0))

		avgRnaCounts = np.average(avgRnaCounts, axis=0)
		avgMonomerCountsFromDeComplexing = np.average(avgMonomerCountsFromDeComplexing, axis=0)
		avgComplexCounts = np.average(avgComplexCounts, axis=0)

		# Compute transcription frequency (over generations)
		transcriptionFreq = np.mean(rnaSumIsNonZero, axis=0)

		# Generate transcription frequency color-code for use in plotting
		colors = [COLOR_FREQ_0 if freq == 0 else COLOR_FREQ_1 if freq == 1 else COLOR_FREQ_SUBGEN for freq in transcriptionFreq]

		# Plot
		fig, axesList = plt.subplots(1, 2, figsize = (11, 8.5))
		ax0, ax1 = axesList

		# Subplot 0: average RNA vs average monomers (sum from de-complexing) per cell
		avgRnaCounts = replace_zeros(avgRnaCounts)
		avgMonomerCountsFromDeComplexing = replace_zeros(avgMonomerCountsFromDeComplexing)

		ax0.scatter(np.log10(avgRnaCounts), np.log10(avgMonomerCountsFromDeComplexing), c = colors, s = 8, alpha = 0.4)
		ax0.set_xlabel("log10 <RNA> per cell")
		ax0.set_ylabel("log10 <Monomers> per cell")
		ax0.set_title("Monomers (sum from de-complexing)")

		# Subplot 1: average RNA vs average functional units
		avgComplexCounts = replace_zeros(avgComplexCounts)
		translation = sim_data.process.translation
		for complexId, counts in zip(complexIds, avgComplexCounts):
			subunitIds = sim_data.process.complexation.getMonomers(complexId)["subunitIds"]

			for subunitId in subunitIds:

				# Skip RNA subunits in complexes
				if subunitId not in translation.monomerData["id"]:
					continue

				index = np.where(translation.monomerData["id"] == subunitId)[0][0]
				freq = transcriptionFreq[index]
				color = COLOR_FREQ_0 if freq == 0 else COLOR_FREQ_1 if freq == 1 else COLOR_FREQ_SUBGEN
				avgRnaCount = avgRnaCounts[index]
				ax1.scatter(np.log10(avgRnaCount), np.log10(counts), c=color, s=8, alpha=0.4)

		ax1.set_xlabel("log10 <RNA> per cell")
		ax1.set_ylabel("log10 <Functional Unit> per cell")
		ax1.set_title("Functional Units"
					  "\n(x points per gene,"
					  "\nx = num complexes the monomer participates in)"
					  "\n(avg across {0} generations)".format(len(allDir)))

		# Format plots
		for ax in axesList:
			ax.tick_params(which = "both", direction = "out")
			ax.axvline(np.log10(PLOT_ZEROS_ON_LINE), linewidth = 0.1, color = "k")
			ax.axhline(np.log10(PLOT_ZEROS_ON_LINE), linewidth = 0.1, color = "k")
		plt.subplots_adjust(hspace = 0.5, wspace = 0.5, left = 0.1, bottom = 0.3, top = 0.7, right = 0.95)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
