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


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)

		# To handle the --operons=both case, a variant analysis class needs to
		# * Call ap = AnalysisPaths(all_variant_plot=True).
		#   *OR* call ap = AnalysisPaths(variant_plot=True) and handle only the
		#   primary case, which in an --operons=both run is operons OFF.
		# * Don't assume the combined variant index values are contiguous, e.g.
		#   don't use range(ap.n_variant).
		# * Call filepath.split_variant_index() to split the combined variant
		#   index (each element of ap.get_variants()) into the operon part and
		#   the base variant index.
		# * Use the primary (kb/) or secondary (kb-poly/) sim_data files and
		#   validation_data files as appropriate per variant.
		# * See AnalysisPaths utilities .variants, .base_variants, .operons,
		#   .get_cell_seed(), .get_cell_generation(), etc.
		sim_data1, sim_data2 = ap.read_sim_data_files(simDataFile)
		validation_data1, validation_data2 = ap.read_validation_data_files(
			validationDataFile)
		variants = ap.get_variants()

		for variant in variants:
			# Load modified variant sim_data
			with open(ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)

			cell_paths = ap.get_cells(variant=[variant])

			# Load data
			## Simple stacking functions for data from all cells
			names = ['ATP[c]']  # Replace with desired list of names
			time = read_stacked_columns(cell_paths, 'Main', 'time')
			(counts,) = read_stacked_bulk_molecules(cell_paths, (names,))

			## Or iterate on each cell if additional processing is needed
			for sim_dir in cell_paths:
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Listeners used
				main_reader = TableReader(os.path.join(simOutDir, 'Main'))

				# Load data
				time = main_reader.readColumn('time')

				(counts,) = read_bulk_molecule_counts(simOutDir, (names,))

			_ = time, counts  # Avoid "unused" warnings in this template.

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Avoid "unused" warnings in this template.
		_ = np, sim_data1, sim_data2, validation_data1, validation_data2
		_ = variant_sim_data


if __name__ == "__main__":
	Plot().cli()
