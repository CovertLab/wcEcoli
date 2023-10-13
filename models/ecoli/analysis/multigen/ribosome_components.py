"""
Plots the timetrace of counts for each of the components of the ribosomal
subunits (rRNAs and ribosomal proteins).
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		# Load IDs of ribosome components from sim_data
		s30_protein_ids = sim_data.molecule_groups.s30_proteins
		s30_16s_rRNA_ids = sim_data.molecule_groups.s30_16s_rRNA
		s30_full_complex_id = [sim_data.molecule_ids.s30_full_complex]
		s50_protein_ids = sim_data.molecule_groups.s50_proteins
		s50_23s_rRNA_ids = sim_data.molecule_groups.s50_23s_rRNA
		s50_5s_rRNA_ids = sim_data.molecule_groups.s50_5s_rRNA
		s50_full_complex_id = [sim_data.molecule_ids.s50_full_complex]

		# Read free counts of all components
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(s30_protein_counts, s30_16s_rRNA_counts, s50_protein_counts,
			s50_23s_rRNA_counts, s50_5s_rRNA_counts
			) = read_stacked_bulk_molecules(
			cell_paths,
			(s30_protein_ids, s30_16s_rRNA_ids, s50_protein_ids,
			s50_23s_rRNA_ids, s50_5s_rRNA_ids))

		# Read counts of free 30S and 50S subunits
		(s30_full_complex_counts, s50_full_complex_counts
			) = read_stacked_bulk_molecules(
			cell_paths,(s30_full_complex_id, s50_full_complex_id))

		# Read counts of active ribosomes
		simOutDir = os.path.join(cell_paths[0], 'simOut')
		unique_molecule_counts_reader = TableReader(
			os.path.join(simOutDir, 'UniqueMoleculeCounts'))
		active_ribosome_index = unique_molecule_counts_reader.readAttribute(
			'uniqueMoleculeIds').index('active_ribosome')
		active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				)[:, active_ribosome_index]

		# Calculate total counts of all components
		s30_limiting_protein_counts = (
			s30_protein_counts.min(axis=1)
			+ s30_full_complex_counts + active_ribosome_counts)
		s30_16s_rRNA_total_counts = (
			s30_16s_rRNA_counts.sum(axis=1)
			+ s30_full_complex_counts + active_ribosome_counts)
		s50_limiting_protein_counts = (
			s50_protein_counts.min(axis=1)
			+ s50_full_complex_counts + active_ribosome_counts)
		s50_23s_rRNA_total_counts = (
			s50_23s_rRNA_counts.sum(axis=1)
			+ s50_full_complex_counts + active_ribosome_counts)
		s50_5s_rRNA_total_counts = (
			s50_5s_rRNA_counts.sum(axis=1)
			+ s50_full_complex_counts + active_ribosome_counts)

		# Plot timetraces of all component counts
		plt.figure(figsize=(10, 3))

		# 30S components
		ax1 = plt.subplot(2, 1, 1)
		ax1.plot(time / 60, s30_limiting_protein_counts, clip_on=False,
			label='limiting r-protein', c='#cccccc', ls='--', lw=2.5)
		ax1.plot(time / 60, s30_16s_rRNA_total_counts, clip_on=False,
			label='16S rRNA', c='C0')
		ax1.set_ylabel('30S component\ncounts')
		ax1.spines["top"].set_visible(False)
		ax1.spines["right"].set_visible(False)
		ax1.spines["bottom"].set_position(("outward", 10))
		ax1.spines["left"].set_position(("outward", 10))
		ax1.spines["bottom"].set_visible(False)
		ax1.get_xaxis().set_visible(False)
		ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
		ax1.set_xlim([0, time[-1] / 60])
		ax1.set_ylim([0, 60000])

		# 50S components
		ax2 = plt.subplot(2, 1, 2)
		ax2.plot(time / 60, s50_limiting_protein_counts, clip_on=False,
			label='limiting r-protein', c='#cccccc', ls='--', lw=2.5)
		ax2.plot(time / 60, s50_23s_rRNA_total_counts, clip_on=False,
			label='23S rRNA', c='C1')
		ax2.plot(time / 60, s50_5s_rRNA_total_counts, clip_on=False,
			label='5S rRNA', c='C2')
		ax2.set_xlabel('Time (min)')
		ax2.set_ylabel('50S component\ncounts')
		ax2.spines["top"].set_visible(False)
		ax2.spines["right"].set_visible(False)
		ax2.spines["bottom"].set_position(("outward", 10))
		ax2.spines["left"].set_position(("outward", 10))
		ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
		ax2.set_xlim([0, time[-1] / 60])
		ax2.set_ylim([0, 60000])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
