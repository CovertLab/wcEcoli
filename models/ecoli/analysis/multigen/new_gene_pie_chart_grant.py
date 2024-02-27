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
		labels = ["Exogenous Gene", "Native Genes"]
		colors = [SLICE_COLOR, "#dfdfdf"]

		# Plotting

		# Pie Chart
		plt.figure(figsize = (4, 4))
		fig, ax = plt.subplots()
		ax.pie(slice_sizes, labels=labels, colors=colors) # autopct='%1.1f%%'
		exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
		plt.close("all")

		# Stacked Bar Chart
		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False
		mpl.rcParams['axes.spines.left'] = False
		mpl.rcParams['axes.spines.bottom'] = False

		fig, ax = plt.subplots(figsize=(8.5/2, 11/5))
		ax.invert_yaxis()
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
		ax.set_xlim(0, 100)

		bar_widths = np.array(slice_sizes)
		bar_cum = bar_widths.cumsum()
		bar_starts = bar_cum - bar_widths
		for i in range(len(slice_sizes)):
			rects = ax.barh(
				0, bar_widths[i], left = bar_starts[i],
				height = 0.5, label=labels[i], color=colors[i])
			if i == (len(slice_sizes) - 1):
				text_color = colors[0]
				ax.bar_label(rects, labels = [str(slice_sizes[0]) + "%"], label_type='center', color=text_color)
		ax.legend(ncols=len(labels), bbox_to_anchor=(0, 1),
				  loc='lower left', fontsize='small')

		exportFigure(plt, plotOutDir, plotOutFileName + "_bar" + plot_suffix, metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
