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
		data_desc = "_loss_ribo_detail" # TODO: CHANGE
		plot_suffix = data_desc
		data = { # TODO: CHANGE
			"New Gene": np.array([ 0, 3367, 3778]),
			"RNAP Subunits": np.array([204, 45, 79]),
			"Ribosomal Proteins": np.array([3114, 720, 1162]),
			"rRNA": np.array([0, 0, 0]),
			"Other": np.array([13955, 3582, 5751])
			}
		ylab = "Ribosome Counts" # TODO: CHANGE
		species = (
			"Wildtype (no GFP)",
			"Normal GFP Degradation",
			"Rapid GFP Degradation"
			)
		colors = [
			(66 / 255, 170 / 255, 154 / 255),
			(136/255, 205/255, 240/255),
			(188/255, 140/255, 191/255),
			(221/255, 203/255, 119/255),
			"#dfdfdf"
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

		fig, ax = plt.subplots()
		width = 0.5
		counter = 0
		bottom = np.zeros(3)
		# plt.bar(species, list(data.values()), color=colors)
		for category, allocation in data.items():
			p = ax.bar(species, allocation, width, label=category, bottom=bottom, color = colors[counter])
			bottom += allocation
			counter += 1
		plt.ylabel(ylab)
		ax.legend(loc="upper center")
			# if i == 0:
			# 	text_color = colors[0]
			# 	ax.bar_label(rects, labels = [" -" + str(percent_loss) + "%"], label_type='edge', color=text_color, fontsize=14)
		# ax.legend(ncols=len(labels), bbox_to_anchor=(0, 1),
		# 		  loc='lower left', fontsize='small')

		exportFigure(plt, plotOutDir, plotOutFileName + "_bar" + plot_suffix, metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
