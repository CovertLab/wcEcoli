"""
Template for variant analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


def plot(ax, x, y, xlabel, ylabel):
	ax.plot(x[0], y[0], 'og')
	ax.plot(x[-1], y[-1], 'or')
	ax.plot(x, y)

	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel(ylabel, fontsize=8)
	ax.tick_params(labelsize=6)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		# Create plot
		_, axes = plt.subplots(3, 1, figsize=(4, 12))

		for variant in variants:
			for seed in ap.get_seeds(variant):
				cell_paths = ap.get_cells(variant=[variant], seed=[seed])

				# Load data
				growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True).squeeze()
				protein = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
				rna = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
				growth_means = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True, fun=np.mean).squeeze()
				protein_means = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True, fun=np.mean).squeeze()
				rna_means = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True, fun=np.mean).squeeze()
				mass_means = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True, fun=np.mean).squeeze()

				# Process data
				ratio = rna / protein
				growth_ma = np.convolve(growth, np.ones(200) / 200, mode='valid')
				ratio_ma = np.convolve(ratio, np.ones(200) / 200, mode='valid')

				plot(axes[0], mass_means, growth_means, 'Cell cycle mass', 'Cell cycle growth')
				plot(axes[1], rna_means / protein_means, growth_means, 'Cell cycle RNA/protein', 'Cell cycle growth')
				# plot(axes[2], ratio, growth, 'RNA/protein', 'Growth')
				plot(axes[2], ratio_ma, growth_ma, 'RNA/protein', 'Growth')

		for ax in axes:
			self.remove_border(ax)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
