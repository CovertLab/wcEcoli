"""
Plot count of stalled RNA polymerases from stalled transcription elongation
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import range


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		stallRNAP_counts_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
		stallRNAP_counts = stallRNAP_counts_reader.readColumn('didStall')

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		plt.plot(time / 60., stallRNAP_counts)
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("Stalled RNA Polymerase")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
