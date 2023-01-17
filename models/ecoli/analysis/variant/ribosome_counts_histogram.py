"""
Plot number of ribosomes for all generations, early generations, and/or late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, index_of_first, labeled_indexable_hist
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS

exclude_timeout_cells = 1 # 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_early_gens = 1 # 1 to plot early (before MIN_LATE_CELL_INDEX), and late generationss in addition to all generations

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index
MAX_CELL_INDEX = 8 # do not include any generation >= this index
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized
MAX_CELL_LENGTH = 180
if exclude_timeout_cells:
	MAX_CELL_LENGTH += 1000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		print("Running analysis script with exclude_timeout_cells=", exclude_timeout_cells,
			  " and exclude_early_gens=", exclude_early_gens)

		# Data extraction
		print("---Data Extraction---")
		ribosome_counts = {}

		variants = self.ap.get_variants()
		min_variant = min(variants)
		for variant in variants:
			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			exclude_timeout_index = index_of_first(dt,MAX_CELL_LENGTH)

			# Ribosome counts
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')
				uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')

			avg_ribosome_counts = read_stacked_columns(all_cells, 'UniqueMoleculeCounts','uniqueMoleculeCounts',
													   fun=lambda x: np.mean(x[:,ribosomeIndex],axis=0))
			ribosome_counts[variant] = avg_ribosome_counts[:exclude_timeout_index]

		# Plotting
		print("---Plotting---")
		std_bin_width = 250
		std_sf = 0
		std_xlim = [10000,28000]
		std_xlab = 'Active Ribosome Counts'

		data_start = [0]
		data_end = [MAX_CELL_INDEX]
		plot_label = ['_all_gens']
		if exclude_early_gens:
			data_start += [0,MIN_LATE_CELL_INDEX]
			data_end += [MIN_LATE_CELL_INDEX,MAX_CELL_INDEX]
			plot_label += ['_early_gens', '_late_gens']

		for j in range(len(data_start)):
			_, axes = plt.subplots(1, 1, figsize=(10, 5))
			labeled_indexable_hist(self, axes, ribosome_counts, data_start[j], data_end[j], COLORS,
								   std_xlab, bin_width=std_bin_width, sf=std_sf)
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + plot_label[j], metadata)

			axes.set_xlim(std_xlim)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_label[j] + '_trimmed', metadata)

if __name__ == "__main__":
	Plot().cli()
