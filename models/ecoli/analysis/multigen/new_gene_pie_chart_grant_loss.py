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

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		data_desc = "_loss_rnaps" # TODO: CHANGE
		plot_suffix = data_desc
		data = { # TODO: CHANGE
			"No GFP Production": 3895,
			"Lower GFP Production": 3488,
			"Higher GFP Production": 1928
			}
		ylab = "RNAP Counts" # TODO: CHANGE
		colors = [
			(136/255, 205/255, 240/255),
			(188/255, 140/255, 191/255),
			(221/255, 203/255, 119/255)
			]

		# Plotting
		# # Pie Chart
		# plt.figure(figsize = (4, 4))
		# fig, ax = plt.subplots()
		# ax.pie(slice_sizes, labels=labels, colors=colors) # autopct='%1.1f%%'
		# exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
		# plt.close("all")

		# Stacked Bar Chart
		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False
		# mpl.rcParams['axes.spines.left'] = False
		# mpl.rcParams['axes.spines.bottom'] = False

		fig = plt.figure(figsize=(7.25, 6))
		plt.bar(list(data.keys()), list(data.values()), color=colors)
		plt.ylabel(ylab)
			# if i == 0:
			# 	text_color = colors[0]
			# 	ax.bar_label(rects, labels = [" -" + str(percent_loss) + "%"], label_type='edge', color=text_color, fontsize=14)
		# ax.legend(ncols=len(labels), bbox_to_anchor=(0, 1),
		# 		  loc='lower left', fontsize='small')

		exportFigure(plt, plotOutDir, plotOutFileName + "_bar" + plot_suffix, metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
