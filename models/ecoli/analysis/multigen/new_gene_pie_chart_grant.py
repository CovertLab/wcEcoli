"""
Plot mRNA and protein counts for new genes across multiple generations, as well
as plots to analyze the impact of new gene expression, including growth rate,
RNAP and ribosome counts, and ppGpp concentration.
"""
# TODO: update file header comment

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
# noinspection PyUnresolvedReferences
import numpy as np
from numpy import inf

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

SLICE_COLOR = (66/255, 170/255, 154/255)

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		variant = 7
		slice_sizes = [0.023 * 100, 100 - 0.023 * 100]
		plot_suffix = "_" + str(variant) + "_rnap_fraction"
		labels = ["gfp", "other"]

		# Plotting
		plt.figure(figsize = (4, 4))

		# plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
		fig, ax = plt.subplots()
		ax.pie(slice_sizes, labels=labels, autopct='%1.1f%%')

		exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
