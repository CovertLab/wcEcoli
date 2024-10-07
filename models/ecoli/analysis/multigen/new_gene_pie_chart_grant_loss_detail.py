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

		# RIBOSOMES
		data_desc = "_loss_ribo_detail" # TODO: CHANGE
		plot_suffix = data_desc
		data = { # TODO: CHANGE
			"New Gene": np.array([ 0, 2838, 3211]),
			"RNAP Subunits": np.array([193, 100, 34]),
			"Ribosomal Proteins": np.array([2997, 1594, 594]),
			"rRNA": np.array([0, 0, 0]),
			"Other": np.array([13655, 7503, 3000])
			}
		ylab = "Active Ribosome Counts" # TODO: CHANGE
		species = (
			"Wildtype (no GFP)",
			"Lower GFP Production", # Variant 6
			"Higher GFP Production" # Variant 16
			)

		# # RNAPs
		# data_desc = "_loss_rnap_detail" # TODO: CHANGE
		# plot_suffix = data_desc
		# data = { # TODO: CHANGE
		# 	"New Gene": np.array([ 0, 12, 19]),
		# 	"RNAP Subunits": np.array([7, 5, 2]),
		# 	"Ribosomal Proteins": np.array([65, 41, 21]),
		# 	"rRNA": np.array([301, 185, 89]),
		# 	"Other": np.array([393, 249, 133])
		# 	}
		# ylab = "Active RNA Polymerase Counts" # TODO: CHANGE
		# species = (
		# 	"Wildtype (no GFP)",
		# 	"Lower GFP Production", # Variant 6
		# 	"Higher GFP Production" # Variant 16
		# 	)

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
