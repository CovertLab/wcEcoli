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
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		growth_limits_reader = TableReader(os.path.join(simOutDir, 'GrowthLimits'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		ppgpp_concentration = growth_limits_reader.readColumn('ppgpp_conc')

		plt.figure()
		plt.plot(time / 60., ppgpp_concentration, label='ppGpp concentration')
		plt.xlabel('Time (min)')
		plt.ylabel('ppGpp concentration (uM)')

		### Create Plot ###
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == '__main__':
	Plot().cli()
