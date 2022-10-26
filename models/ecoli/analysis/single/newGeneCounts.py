"""
Plot mRNA and protein counts for new genes
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import range


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		### TODO: CHECK THAT THE NEW GENES VARIANT WAS ON FOR THIS SIM

		# Extract mRNA counts for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
		### TODO: FIND THESE NEW GENE IDS USING NAMING CONVENTION NGXXX
		new_gene_mRNA_IDS = ['NG001_RNA[c]']
		new_gene_mrna_indexes = np.array([mRNA_idx[newGeneId] for newGeneId in new_gene_mRNA_IDS], int)
		new_gene_mrna_counts = mRNA_counts[:, new_gene_mrna_indexes]

		# Extract protein counts for each new gene
		monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		monomer_counts = monomer_counts_reader.readColumn('monomerCounts')
		monomer_idx = {monomer: i for i, monomer in enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
		### TODO: FIND THESE NEW GENE IDS USING NAMING CONVENTION NGXXX
		new_gene_monomer_IDS = ['GFP-MONOMER[m]']
		new_gene_monomer_indexes = np.array([monomer_idx[newGeneId] for newGeneId in new_gene_monomer_IDS], int)
		new_gene_monomer_counts = monomer_counts[:, new_gene_monomer_indexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Plotting
		plt.figure(figsize = (8.5, 11))
		### TODO: MAKE THESE PLOTS LOOP THROUGH ALL NEW GENES

		# Protein Counts
		plt.subplot(2, 1, 1)
		plt.plot(time / 60., new_gene_monomer_counts)
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("New Gene Protein Counts")

		# mRNA Counts
		plt.subplot(2, 1, 2)
		plt.plot(time / 60., new_gene_mrna_counts)
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA Counts")
		plt.title("New Gene mRNA Counts")

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()

