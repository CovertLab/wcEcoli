import os

from matplotlib import pyplot as plt
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

START_GEN = 0
END_GEN = 5

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		# Get all cells
		allDir = self.ap.get_cells(generation = np.arange(START_GEN, END_GEN))

		massNames = [
					"cellMass",
					"waterMass",
					"dryMass",
					"smallMoleculeMass",
					"proteinMass",
					"tRnaMass",
					"rRnaMass",
					'mRnaMass',
					"dnaMass",
					"membrane_mass",
					]

		cleanNames = [
					"Cell\nmass",
					"Water\nmass",
					"Dry\nmass",
					"Small\nmolecule\nmass",
					"Protein\nmass",
					"tRNA\nmass",
					"rRNA\nmass",
					"mRNA\nmass",
					"DNA\nmass",
					"Membrane\nmass",
					]

		fig, axesList = plt.subplots(len(massNames), sharex = True)
		fig.set_size_inches(14, 3 * len(massNames))

		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			for idx, massType in enumerate(massNames):
				massToPlot = mass.readColumn(massNames[idx])
				axesList[idx].plot(time / 60. , massToPlot, linewidth = 2)

				axesList[idx].set_ylabel(cleanNames[idx] + " (fg)")

		axesList[0].set_title("Mass by Category")
		axesList[len(massNames) - 1].set_xlabel("Time (min)")

		xticks = list(range(0, int(time[-1] / 60) + 1, 25))
		for ax in axesList:
			ax.set_xticks(xticks)
			# plt.setp(ax.get_xticklabels(), visible=True)
			ax.tick_params(axis='x', labelbottom=True, labelsize=6)
			ax.grid(True)

	
		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)
		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
