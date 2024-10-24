"""
Plots frequency of observing at least 1 protein.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		allDir = self.ap.get_cells()

		proteinPresence = []
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			# Get boolean protein presence
			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			proteinCounts = monomerCounts.readColumn("monomerCounts")
			meanProteinCounts = proteinCounts.mean(axis=0)

			proteinPresence.append(meanProteinCounts != 0)

		proteinPresence = np.array(proteinPresence)

		# Plot
		fig = plt.figure(figsize = (12, 12))
		ax = plt.subplot(1, 1, 1)
		nGens = len(allDir)
		ax.hist(np.mean(proteinPresence, axis = 0), nGens)
		ax.set_xlabel("Frequency of observing at least 1 protein copy in 1 generation", fontsize = 14)
		ax.set_ylabel("Number of proteins", fontsize = 14)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
