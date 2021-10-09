"""
Show the trajectory of growth rates and cell composition measurements.  Useful
for variants with a shift.

TODO
	- select which to show and make a cohort plot?
	- add axis labels
	- moving average window
"""

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


def plot(x, y):
	plt.plot(x[0], y[0], 'og')
	plt.plot(x[-1], y[-1], 'or')
	plt.plot(x, y)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True).squeeze()
		protein = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
		rna = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
		growth_means = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True, fun=np.mean).squeeze()
		protein_means = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True, fun=np.mean).squeeze()
		rna_means = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True, fun=np.mean).squeeze()
		mass_means = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True, fun=np.mean).squeeze()

		# Process data
		ratio = rna / protein
		growth_ma = np.convolve(growth, np.ones(200) / 200, mode='valid')
		ratio_ma = np.convolve(ratio, np.ones(200) / 200, mode='valid')

		# Create plot
		plt.figure(figsize=(4, 12))

		plt.subplot(3, 1, 1)
		plot(mass_means, growth_means)

		plt.subplot(3, 1, 2)
		plot(rna_means / protein_means, growth_means)

		plt.subplot(3, 1, 3)
		plot(ratio, growth)
		plot(ratio_ma, growth_ma)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
