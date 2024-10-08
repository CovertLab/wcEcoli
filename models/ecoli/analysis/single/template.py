"""
Template for single analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = self.read_pickle_file(simDataFile)
		validation_data = self.read_pickle_file(validationDataFile)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		names = ['ATP[c]']  # Replace with desired list of names
		(counts,) = read_bulk_molecule_counts(simOutDir, (names,))

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
