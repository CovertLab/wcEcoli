"""
Compare cell properties at varying levels of ppGpp concentration.  Useful with
with the ppgpp_conc variant and in comparison with data presented in Zhu et al.
2019. https://academic.oup.com/nar/article/47/9/4684/5420536.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


MEAN = 'mean'
STD = 'std'
MARKERS = ['o', 's']  # Expand for more overlays
PANEL_WIDTH = 3.8
PANEL_HEIGHT = 1.7

START_GEN = 0
GFP_ON = 8
GFP_measured = 16
END_GEN = 24


VARIANT_1 = 3
VARIANT_2 = 6

VARIANT_1_SEED = 1
VARIANT_2_SEED = 1

LINE_COLOR = (66/255, 170/255, 154/255)
LINE_COLOR2 = (152/255, 78/255, 163/255)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

			# start here
		cell_paths1 = self.ap.get_cells(
			variant=np.array([VARIANT_1]), generation=np.arange(START_GEN, END_GEN),
			seed=np.array([VARIANT_1_SEED]))
		cell_paths1_1 = self.ap.get_cells(
			variant=np.array([VARIANT_1]), generation=np.arange(START_GEN, GFP_ON),
			seed=np.array([VARIANT_1_SEED]))
		cell_paths1_2 = self.ap.get_cells(
			variant=np.array([VARIANT_1]), generation=np.arange(GFP_ON, GFP_measured),
			seed=np.array([VARIANT_1_SEED]))
		cell_paths1_3 = self.ap.get_cells(
			variant=np.array([VARIANT_1]), generation=np.arange(GFP_measured, END_GEN),
			seed=np.array([VARIANT_1_SEED]))

		cell_paths2 = self.ap.get_cells(
			variant=np.array([VARIANT_2]), generation=np.arange(START_GEN, END_GEN),
			seed=np.array([VARIANT_2_SEED]))
		cell_paths2_1 = self.ap.get_cells(
			variant=np.array([VARIANT_2]), generation=np.arange(START_GEN, GFP_ON),
			seed=np.array([VARIANT_2_SEED]))
		cell_paths2_2 = self.ap.get_cells(
			variant=np.array([VARIANT_2]), generation=np.arange(GFP_ON, GFP_measured),
			seed=np.array([VARIANT_2_SEED]))
		cell_paths2_3 = self.ap.get_cells(
			variant=np.array([VARIANT_2]), generation=np.arange(GFP_measured, END_GEN),
			seed=np.array([VARIANT_2_SEED]))

		# sim_dir = cell_paths1[0]
		# simOutDir = os.path.join(sim_dir, 'simOut')

		# Load data
		# time1 = read_stacked_columns(
		# 	cell_paths1, 'Main', 'time', ignore_exception=True)
		time_no_first1 = read_stacked_columns(
			cell_paths1, 'Main', 'time', remove_first=True, ignore_exception=True)
		# time1_1 = read_stacked_columns(
		# 	cell_paths1_1, 'Main', 'time', ignore_exception=True)
		time_no_first1_1 = read_stacked_columns(
			cell_paths1_1, 'Main', 'time', remove_first=True, ignore_exception=True)
		# time1_2 = read_stacked_columns(
		# 	cell_paths1_2, 'Main', 'time', ignore_exception=True)
		time_no_first1_2 = read_stacked_columns(
			cell_paths1_2, 'Main', 'time', remove_first=True, ignore_exception=True)
		# time1_3 = read_stacked_columns(
		# 	cell_paths1_3, 'Main', 'time', ignore_exception=True)
		time_no_first1_3 = read_stacked_columns(
			cell_paths1_3, 'Main', 'time', remove_first=True, ignore_exception=True)

		# time2 = read_stacked_columns(
		# 	cell_paths2, 'Main', 'time', ignore_exception=True)
		time_no_first2 = read_stacked_columns(
			cell_paths2, 'Main', 'time', remove_first=True, ignore_exception=True)
		# time2_1 = read_stacked_columns(
		# 	cell_paths2_1, 'Main', 'time', ignore_exception=True)
		time_no_first2_1 = read_stacked_columns(
			cell_paths2_1, 'Main', 'time', remove_first=True, ignore_exception=True)
		# time2_2 = read_stacked_columns(
		# 	cell_paths2_2, 'Main', 'time', ignore_exception=True)
		time_no_first2_2 = read_stacked_columns(
			cell_paths2_2, 'Main', 'time', remove_first=True, ignore_exception=True)
		# time2_3 = read_stacked_columns(
		# 	cell_paths2_3, 'Main', 'time', ignore_exception=True)
		time_no_first2_3 = read_stacked_columns(
			cell_paths2_3, 'Main', 'time', remove_first=True, ignore_exception=True)

		# time1 = time1 - time1[0]
		time_no_first1 = time_no_first1 - time_no_first1[0]
		# time1_1 = time1_1 - time1_1[0]
		time_no_first1_1 = time_no_first1_1 - time_no_first1_1[0]
		# time1_2 = time1_2 - time1_2[0]
		time_no_first1_2 = time_no_first1_2 - time_no_first1_2[0]
		# time1_3 = time1_3 - time1_3[0]
		time_no_first1_3 = time_no_first1_3 - time_no_first1_3[0]

		# time2 = time2 - time2[0]
		time_no_first2 = time_no_first2 - time_no_first2[0]
		# time2_1 = time2_1 - time2_1[0]
		time_no_first2_1 = time_no_first2_1 - time_no_first2_1[0]
		# time2_2 = time2_2 - time2_2[0]
		time_no_first2_2 = time_no_first2_2 - time_no_first2_2[0]
		# time2_3 = time2_3 - time2_3[0]
		time_no_first2_3 = time_no_first2_3 - time_no_first2_3[0]

		# ppGpp

		ppGpp_concentration = read_stacked_columns(
			cell_paths1, "GrowthLimits", "ppgpp_conc", remove_first=True,
			ignore_exception=True)

		ppGpp_concentration2 = read_stacked_columns(
			cell_paths2, "GrowthLimits", "ppgpp_conc", remove_first=True,
			ignore_exception=True)

		avg_ppgpp_concentration1_1 = read_stacked_columns(
			cell_paths1_1, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True, fun=lambda x: np.mean(x)).squeeze()
		avg_ppgpp_concentration1_2 = read_stacked_columns(
			cell_paths1_2, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True, fun=lambda x: np.mean(x)).squeeze()
		avg_ppgpp_concentration1_3 = read_stacked_columns(
			cell_paths1_3, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True, fun=lambda x: np.mean(x)).squeeze()

		avg_ppgpp_concentration2_1 = read_stacked_columns(
			cell_paths2_1, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True, fun=lambda x: np.mean(x)).squeeze()
		avg_ppgpp_concentration2_2 = read_stacked_columns(
			cell_paths2_2, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True, fun=lambda x: np.mean(x)).squeeze()
		avg_ppgpp_concentration2_3 = read_stacked_columns(
			cell_paths2_3, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True, fun=lambda x: np.mean(x)).squeeze()

		plt.plot(time_no_first1 / 60., ppGpp_concentration, color=LINE_COLOR)
		plt.plot(time_no_first2 / 60., ppGpp_concentration2, color=LINE_COLOR2)

		plt.plot(time_no_first1_1 / 60., avg_ppgpp_concentration1_1, color=LINE_COLOR)
		plt.plot(time_no_first1_2 / 60., avg_ppgpp_concentration1_2, color=LINE_COLOR)
		plt.plot(time_no_first1_3 / 60., avg_ppgpp_concentration1_3, color=LINE_COLOR)
		plt.plot(time_no_first2_1 / 60., avg_ppgpp_concentration2_1, color=LINE_COLOR)
		plt.plot(time_no_first2_2 / 60., avg_ppgpp_concentration2_2, color=LINE_COLOR)
		plt.plot(time_no_first2_3 / 60., avg_ppgpp_concentration2_3, color=LINE_COLOR)

		# end here


		exportFigure(
			plt, plotOutDir, plotOutFileName +  "_"
							 + "VAR_" + str(VARIANT_1) + "_SEED_" + str(VARIANT_1_SEED) + "_"
							 + "VAR_" + str(VARIANT_2) + "_SEED_" + str(VARIANT_2_SEED) + "_GENS_"
							 + str(START_GEN) + "_" + str(END_GEN), metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
