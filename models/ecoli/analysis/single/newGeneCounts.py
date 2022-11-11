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
from six.moves import cPickle,range


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Extract mRNA counts for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}

		new_gene_mRNA_IDS  = [k for k,v in mRNA_idx.items() if k.startswith('NG')]
		new_gene_mrna_indexes = [v for k, v in mRNA_idx.items() if k.startswith('NG')]

		assert len(new_gene_mRNA_IDS) != 0, 'no new gene mRNAs found'
		new_gene_mrna_counts = mRNA_counts[:, new_gene_mrna_indexes]

		# Extract protein counts for each new gene
		monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		monomer_counts = monomer_counts_reader.readColumn('monomerCounts')
		monomer_idx = {monomer: i for i, monomer in enumerate(monomer_counts_reader.readAttribute('monomerIds'))}

		new_gene_monomer_IDS = [k for k, v in monomer_idx.items() if k.startswith('NG')]
		new_gene_monomer_indexes = [v for k, v in monomer_idx.items() if k.startswith('NG')]

		assert len(new_gene_monomer_IDS) != 0, 'no new gene proteins found'
		new_gene_monomer_counts = monomer_counts[:, new_gene_monomer_indexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Plotting
		plt.figure(figsize = (8.5, 11))

		# Protein Counts
		plt.subplot(2, 1, 1)
		for m in range(len(new_gene_monomer_IDS)):
			plt.plot(time / 60., new_gene_monomer_counts[:,m], label = new_gene_monomer_IDS[m])
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.title("New Gene Protein Counts")
		plt.legend()

		# mRNA Counts
		plt.subplot(2, 1, 2)
		for r in range(len(new_gene_mRNA_IDS)):
			plt.plot(time / 60., new_gene_mrna_counts[:,r], label = new_gene_mRNA_IDS[r])
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA Counts")
		plt.title("New Gene mRNA Counts")
		plt.legend()

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()

