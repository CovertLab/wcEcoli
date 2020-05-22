"""
Template for parca analysis plots

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/22/20
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import filepath


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception('inputDir does not currently exist as a directory')

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
