"""
Plot mRNA and protein counts for new genes across multiple generations
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		### TODO flag new gene mRNAs and proteins more efficiently
		# Extract mRNA ids for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
		new_gene_mRNA_ids = [k for k, v in mRNA_idx.items() if k.startswith('NG')]
		new_gene_mRNA_indexes = [v for k, v in mRNA_idx.items() if k.startswith('NG')]
		assert len(new_gene_mRNA_ids) != 0, 'no new gene mRNAs found'

		# Extract protein ids for each new gene
		monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx = {monomer: i for i, monomer in enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
		new_gene_monomer_ids = [k for k, v in monomer_idx.items() if k.startswith('NG')]
		assert len(new_gene_monomer_ids) != 0, 'no new gene proteins found'

		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(new_gene_monomer_counts,) = read_stacked_bulk_molecules(cell_paths, new_gene_monomer_ids)
		all_mRNA_stacked_counts = read_stacked_columns(cell_paths, 'mRNACounts', 'mRNA_counts')
		new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]

		# Plotting
		plt.figure(figsize = (8.5, 11))

		# Protein Counts
		plt.subplot(2, 1, 1)
		if len(new_gene_monomer_ids):
			plt.plot(time / 60., new_gene_monomer_counts, label=new_gene_monomer_ids[0])
		else:
			for m in range(len(new_gene_monomer_ids)):
				plt.plot(time / 60., new_gene_monomer_counts[:,m], label = new_gene_monomer_ids[m])
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("New Gene Protein Counts")
		plt.legend()

		# mRNA Counts
		plt.subplot(2, 1, 2)
		if len(new_gene_mRNA_ids) == 1:
			plt.plot(time / 60., new_gene_mRNA_counts, label=new_gene_mRNA_ids[0])
		else:
			for r in range(len(new_gene_mRNA_ids)):
				plt.plot(time / 60., new_gene_mRNA_counts[:,r], label = new_gene_mRNA_ids[r])
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA Counts")
		plt.title("New Gene mRNA Counts")
		plt.legend()

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")



if __name__ == '__main__':
	Plot().cli()
