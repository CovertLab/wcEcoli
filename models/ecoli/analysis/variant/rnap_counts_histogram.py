"""
Plot number of ribosomes for all generations, early generations, and/or late
(i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt
import os.path

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure,\
	stacked_cell_identification,\
	read_stacked_bulk_molecules, stacked_cell_threshold_mask
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as \
	COLORS, labeled_indexable_hist

# 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_timeout_cells = 1

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index
MAX_CELL_INDEX = 8 # do not include any generation >= this index

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

MAX_CELL_LENGTH = 180
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# print("Running analysis script with exclude_timeout_cells=",
		# 	  exclude_timeout_cells,
		# 	  " and exclude_early_gens=", exclude_early_gens)

		# Data extraction
		rnap_counts = {}
		generations = {}
		variants = self.ap.get_variants()
		for variant in variants:
			if variant >= MAX_VARIANT:
				continue

			all_cells = self.ap.get_cells(variant=[variant],
										  only_successful=True)
			if len(all_cells) == 0:
				continue

			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', MAX_CELL_LENGTH,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			all_cells_gens = np.array([int(os.path.basename(os.path.dirname(
				cell_path))[-6:]) for cell_path in all_cells])[
				exclude_timeout_cell_mask]
			generations[variant] = all_cells_gens

			# RNA polymerase counts
			rnapId = ["APORNAP-CPLX[c]"]
			(rnapCountsBulk,) = read_stacked_bulk_molecules(all_cells, (rnapId,))

			cell_id_vector = stacked_cell_identification(all_cells, 'Main', 'time')
			cell_ids, idx, cell_total_timesteps = np.unique(
				cell_id_vector, return_inverse=True, return_counts=True)
			sum_rnap_counts = np.bincount(idx, weights=rnapCountsBulk)
			avg_rnap_counts = sum_rnap_counts / cell_total_timesteps

			rnap_counts[variant] = avg_rnap_counts[exclude_timeout_cell_mask]

		# Plotting
		std_bin_width = 250
		std_sf = 0
		std_xlim = [1000,7000]
		std_xlab = 'RNA Polymerase Counts'

		data_start = [0]
		data_end = [MAX_CELL_INDEX]
		plot_label = ['_all_gens']
		if exclude_early_gens:
			data_start += [0,MIN_LATE_CELL_INDEX]
			data_end += [MIN_LATE_CELL_INDEX,MAX_CELL_INDEX]
			plot_label += ['_early_gens', '_late_gens']

		for j in range(len(data_start)):
			_, axes = plt.subplots(1, 1, figsize=(10, 5))
			labeled_indexable_hist(self, axes, rnap_counts, generations,
								   data_start[j], data_end[j], COLORS,
								   std_xlab, bin_width=std_bin_width, sf=std_sf)
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + plot_label[j],
						 metadata)

			axes.set_xlim(std_xlim)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_label[j] +
						 '_trimmed', metadata)

if __name__ == "__main__":
	Plot().cli()
