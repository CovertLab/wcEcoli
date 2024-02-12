"""
Plot histograms of ppGpp concentrations and compares the distributions across
different variants.
"""

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.analysis.plotting_tools import (
	DEFAULT_MATPLOTLIB_COLORS as COLORS)

FONT_SIZE = 9
CONCENTRATION_BOUNDS = [0, 200]  # (uM)
N_BINS = 40

# Set True to exclude cells that hit time limit
EXCLUDE_TIMEOUT_CELLS = True

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 4

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		max_cell_length = 180  # mins
		if not EXCLUDE_TIMEOUT_CELLS:
			max_cell_length += 1000000  # Arbitrary large number

		# Data extraction
		ppgpp_concentrations = {}
		n_total_gens = self.ap.n_generation
		variant_indexes = self.ap.get_variants()

		# Loop through all variant indexes
		for variant_index in variant_indexes:
			# Get all cells (within the generation range) of this variant index
			all_cells = self.ap.get_cells(
				variant=[variant_index],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Get mask to exclude cells that timed out
			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', max_cell_length,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()

			# Get average ppGpp concentrations from each cell
			avg_ppgpp_concentrations = read_stacked_columns(
				all_cells, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True, fun=lambda x: np.mean(x)).squeeze()
			ppgpp_concentrations[variant_index] = avg_ppgpp_concentrations[
				exclude_timeout_cell_mask]

		n_variants = len(ppgpp_concentrations)

		fig = plt.figure(figsize=(8, 2 * (n_variants + 1)))
		bins = np.linspace(
			CONCENTRATION_BOUNDS[0],
			CONCENTRATION_BOUNDS[1],
			N_BINS + 1
			)

		# First subplot overlays the histograms of all variants
		ax0 = fig.add_subplot(n_variants + 1, 1, 1)
		subplots = []

		for i, (variant_index, rc) in enumerate(ppgpp_concentrations.items()):
			color = COLORS[i % len(COLORS)]
			ax0.hist(
				rc, bins=bins, color=color, alpha=0.5,
				label=f'Var {variant_index} (n={len(rc)}, {np.mean(rc):.1f} $\pm$ {np.std(rc):.1f})')

			# Add vertical line at mean ppGpp concentrations
			ax0.axvline(np.mean(rc), color=color, ls='--', lw=3, alpha=0.5)

			# Later subplots show the histograms for each variant separately
			ax = fig.add_subplot(n_variants + 1, 1, i + 2, sharex=ax0)
			ax.hist(
				rc, bins=bins, color=color, alpha=0.5,
				)
			ax.tick_params(labelsize=FONT_SIZE)
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			subplots.append(ax)

		# Set ylims of all subplots to be equal to the first subplot
		ax0_ylim = ax0.get_ylim()
		for ax in subplots:
			ax.set_ylim(ax0_ylim)

		ax0.legend()

		ax0.set_xlim(*CONCENTRATION_BOUNDS)
		ax0.set_xticks(
			np.linspace(CONCENTRATION_BOUNDS[0], CONCENTRATION_BOUNDS[1], 11)
			)
		ax0.tick_params(labelsize=FONT_SIZE)

		ax0.set_xlabel('ppGpp concentrations (uM)', fontsize=FONT_SIZE)
		ax0.spines["top"].set_visible(False)
		ax0.spines["right"].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
