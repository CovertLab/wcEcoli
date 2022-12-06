"""
Plot number of ribosomes for all generations, early generations, and late (i.e. not early) generations
"""

import numpy as np
from matplotlib import pyplot as plt

import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, stacked_cell_identification
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS as COLORS


FONT_SIZE=9
MIN_LATE_CELL_INDEX = 4 # generations before this may not be representative of dynamics due to how they are initialized

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def hist(self, ax, data, xlabel, bin_width=1., xlim=None, sf=1):
		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			bins = max(1, int(np.ceil((variant_data.max() - variant_data.min()) / bin_width)))
			mean = variant_data.mean()
			std = variant_data.std()
			ax.hist(variant_data, bins, color=color, alpha=0.5,
				label=f'Var {variant}: {mean:.{sf}f} +/- {std:.{sf+1}f}')
			ax.axvline(mean, color=color, linestyle='--', linewidth=1)

		if xlim:
			ax.set_xlim(xlim)
		self.remove_border(ax)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ribosome_counts = {}
		ribosome_counts_early_gens = {}
		ribosome_counts_late_gens = {}

		print("---Data Extraction---")
		for variant in self.ap.get_variants():

			if variant >= 10:
				continue

			print("Variant: ", variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			if variant == 0:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')
				uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')

			unique_molecule_stacked_counts = read_stacked_columns(all_cells,'UniqueMoleculeCounts','uniqueMoleculeCounts')
			ribosome_stacked_counts = unique_molecule_stacked_counts[:,ribosomeIndex]

			cell_id_vector = stacked_cell_identification(all_cells, 'Main', 'time')
			cell_ids, idx, cell_total_timesteps = np.unique(cell_id_vector, return_inverse=True, return_counts=True)
			sum_ribosome_counts = np.bincount(idx, weights=ribosome_stacked_counts)
			avg_ribosome_counts = sum_ribosome_counts / cell_total_timesteps

			ribosome_counts[variant] = avg_ribosome_counts

			all_cells_gens = [int(c.split("/")[-2][-6:]) for c in all_cells]
			early_cell_index = [i for i,v in enumerate(all_cells_gens) if v < MIN_LATE_CELL_INDEX]
			late_cell_index = [i for i,v in enumerate(all_cells_gens) if v >= MIN_LATE_CELL_INDEX]

			ribosome_counts_early_gens[variant] = ribosome_counts[variant][early_cell_index]
			ribosome_counts_late_gens[variant] = ribosome_counts[variant][late_cell_index]

		print("---Plotting---")
		# LATE GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], ribosome_counts_late_gens, 'Active Ribosome Counts',bin_width=250,sf=0)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_lategens', metadata)

		axes[0].set_xlim([10000, 28000])
		exportFigure(plt, plotOutDir, plotOutFileName + '_lategens_trimmed', metadata)

		# EARLY GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], ribosome_counts_early_gens, 'Active Ribosome Counts',bin_width=250,sf=0)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_earlygens', metadata)

		axes[0].set_xlim([10000, 28000])
		exportFigure(plt, plotOutDir, plotOutFileName + '_earlygens_trimmed', metadata)

		# ALL GENS
		_, axes = plt.subplots(2, 1, figsize=(10, 10))
		self.hist(axes[0], ribosome_counts, 'Active Ribosome Counts',bin_width=250,sf=0)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_allgens', metadata)

		axes[0].set_xlim([10000, 28000])
		exportFigure(plt, plotOutDir, plotOutFileName + '_allgens_trimmed', metadata)


if __name__ == "__main__":
	Plot().cli()
