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
	read_stacked_bulk_molecules, read_stacked_columns, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

LINE_COLOR = (66/255, 170/255, 154/255)

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)

		plot_suffixes = [""]
		standard_xlim = (0,2000)
		total_plots = 100 # TODO Modularize and get rid of this magic number

		for i in range(len(plot_suffixes)):

			# Get time marker where GFP induced
			dt = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			num_time_steps = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: len(x)).squeeze()
			dt_to_plot = np.repeat(dt, num_time_steps)
			doubling_time_index = np.where(dt_to_plot == dt[8])[0][0]
			analysis_gen_index = np.where(dt_to_plot == dt[16])[0][0]

			# tRNA counts
			transcription = sim_data.process.transcription
			uncharged_trna_ids = transcription.uncharged_trna_names
			charged_trna_ids = transcription.charged_trna_names

			# TODO: delete after testing
			uncharged_trna_ids = uncharged_trna_ids[0:5]
			charged_trna_ids = charged_trna_ids[0:5]

			def plot_trna_counts(trna_ids, title_prefix, suffix, time, cell_paths, output_dir,
								 output_name, line_color, total_plots, standard_xlim,
								 doubling_time_index,
								 analysis_gen_index, ax_title, metadata):
				mpl.rcParams['axes.spines.right'] = False
				mpl.rcParams['axes.spines.top'] = False
				plt.figure(figsize=(12, total_plots * 3))
				ax1 = plt.subplot(total_plots, 1, 1)
				plot_num = 1

				for trna_id in trna_ids:
					plt.subplot(total_plots, 1, plot_num, sharex=ax1)
					plt.axvline(time[doubling_time_index] / 60., color='gray', linestyle='--',
							   lw=0.5)
					plt.axvline(time[analysis_gen_index] / 60., color='gray', linestyle='--', lw=0.5)

					counts = read_stacked_bulk_molecules(cell_paths, ([trna_id],),
														 ignore_exception=True)
					plt.plot(time / 60., counts[0].T, color=line_color)

					if suffix in ("_standard_axes_both", "_standard_axes_y"):
						plt.ylim(-1, 4.5)
					if suffix in ("_standard_axes_both", "_standard_axes_x"):
						plt.xlim(standard_xlim)

					plt.xlabel("Time (min)")
					plt.ylabel("Counts: " + trna_id, fontsize="x-small")
					plt.title(f"{title_prefix} tRNA Counts: {trna_id}")
					plot_num += 1

				print(f"Total number of plots made: {len(trna_ids)}")
				plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
				exportFigure(plt, output_dir, output_name + suffix, metadata)
				plt.close("all")

			for plot_suffix in plot_suffixes:
				plot_trna_counts(
					uncharged_trna_ids, "Uncharged", plot_suffix, time, cell_paths,
					plotOutDir, plotOutFileName + "_uncharged_tRNA", LINE_COLOR,
					total_plots, standard_xlim, doubling_time_index, analysis_gen_index,
					"Uncharged", metadata
					)

				plot_trna_counts(
					charged_trna_ids, "Charged", plot_suffix, time, cell_paths,
					plotOutDir, plotOutFileName + "_charged_tRNA", LINE_COLOR,
					total_plots, standard_xlim, doubling_time_index, analysis_gen_index,
					"Charged", metadata
					)


if __name__ == '__main__':
	Plot().cli()
