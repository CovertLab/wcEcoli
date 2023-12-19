"""
Plots the number of collisions between RNAPs and replisomes that occur on each
gene. Only the top N genes with the most collisions are plotted.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


PLOT_TOP_N_GENES = 25

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells(only_successful=True)
		n_cells = len(cell_paths)

		# Load data
		n_headon = read_stacked_columns(
			cell_paths, 'RnapData', 'n_headon_collisions')
		n_codirectional = read_stacked_columns(
			cell_paths, 'RnapData', 'n_codirectional_collisions')

		# Get total number of collisions per cell
		n_headon_per_cell = n_headon.sum() / n_cells
		n_codirectional_per_cell = n_codirectional.sum() / n_cells

		# Initialize list of collision coordinates
		headon_coordinates = []
		codirectional_coordinates = []

		for sim_dir in cell_paths:
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Read collision coordinates for this cell
			rnap_data_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
			headon_this_cell = rnap_data_reader.readColumn(
				'headon_collision_coordinates')
			codirectional_this_cell = rnap_data_reader.readColumn(
				'codirectional_collision_coordinates')

			# Flatten collision coordinates
			headon_this_cell = headon_this_cell[
				np.logical_not(np.isnan(headon_this_cell))].flatten()
			codirectional_this_cell = codirectional_this_cell[
				np.logical_not(np.isnan(codirectional_this_cell))].flatten()

			# Add to list
			headon_coordinates.extend(headon_this_cell.tolist())
			codirectional_coordinates.extend(codirectional_this_cell.tolist())

		headon_coordinates = np.array(headon_coordinates)
		codirectional_coordinates = np.array(codirectional_coordinates)

		# Get gene data from sim_data
		gene_ids = sim_data.process.transcription.cistron_data['gene_id']
		gene_coordinates = sim_data.process.transcription.cistron_data['replication_coordinate']
		gene_directions = sim_data.process.transcription.cistron_data['is_forward']
		gene_lengths = sim_data.process.transcription.cistron_data["length"].asNumber(units.nt)
		is_rRNA = sim_data.process.transcription.cistron_data['is_rRNA']

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Compute boundaries for each gene
		gene_boundaries = []
		for coord, is_forward, length in zip(
				gene_coordinates, gene_directions, gene_lengths):
			if is_forward:
				# Check for genes that loop around terC
				assert coord + length <= replichore_lengths[0]
				gene_boundaries.append((coord, coord + length))
			else:
				assert coord - length >= -replichore_lengths[1]
				gene_boundaries.append((coord - length, coord))

		# Initialize arrays
		n_codirectional_per_gene = np.zeros_like(gene_coordinates)
		n_headon_per_gene = np.zeros_like(gene_coordinates)

		# Identify which gene each collision occurred in
		for coord in codirectional_coordinates:
			for i, boundaries in enumerate(gene_boundaries):
				if boundaries[0] < coord < boundaries[1]:
					n_codirectional_per_gene[i] += 1
					continue

		for coord in headon_coordinates:
			for i, boundaries in enumerate(gene_boundaries):
				if boundaries[0] < coord < boundaries[1]:
					n_headon_per_gene[i] += 1
					continue

		# Divide by total number of cells
		n_codirectional_per_gene = n_codirectional_per_gene / n_cells
		n_headon_per_gene = n_headon_per_gene / n_cells

		# Sort by number of collisions
		codirectional_rank = np.argsort(n_codirectional_per_gene)[::-1][:PLOT_TOP_N_GENES]
		headon_rank = np.argsort(n_headon_per_gene)[::-1][:PLOT_TOP_N_GENES]

		# Get common names of top N genes
		codirectional_top_genes = [sim_data.common_names.get_common_name(gene)
			for gene in gene_ids[codirectional_rank]]
		headon_top_genes = [sim_data.common_names.get_common_name(gene)
			for gene in gene_ids[headon_rank]]
		codirectional_rRNA_indexes = np.where(is_rRNA[codirectional_rank])[0]
		headon_rRNA_indexes = np.where(is_rRNA[headon_rank])[0]

		# Plot
		plt.figure(figsize=(13, 3.5))

		ax = plt.subplot(1, 2, 1)
		ax.bar(list(range(PLOT_TOP_N_GENES)), n_codirectional_per_gene[codirectional_rank],
			   color="darkblue")
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(codirectional_top_genes, rotation=90)
		for i in codirectional_rRNA_indexes:
			ax.get_xticklabels()[i].set_color('red')
		ax.set_title("Co-directional (Average per cell = %.1f)" % (n_codirectional_per_cell,))
		ax.set_ylabel("Number of collisions")
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		ax = plt.subplot(1, 2, 2)
		ax.bar(range(PLOT_TOP_N_GENES), n_headon_per_gene[headon_rank],
			   color="crimson")
		ax.set_xticks(range(PLOT_TOP_N_GENES))
		ax.set_xticklabels(headon_top_genes, rotation=90)
		for i in headon_rRNA_indexes:
			ax.get_xticklabels()[i].set_color('red')
		ax.set_title("Head-on (Average per cell = %.1f)" % (n_headon_per_cell,))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
