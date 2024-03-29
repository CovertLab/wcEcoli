"""
Plots counts of 30S rRNA, associated proteins, and complexes
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import sparklineAxis, setAxisMaxMinY
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot

FONT = {
	'size':	8
	}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Load data from KB
		sim_data = self.read_pickle_file(simDataFile)
		proteinIds = sim_data.molecule_groups.s30_proteins
		cistron_ids = [sim_data.process.translation.monomer_data['cistron_id'][np.where(sim_data.process.translation.monomer_data['id'] == pid)[0][0]] for pid in proteinIds]
		rRnaIds = sim_data.molecule_groups.s30_16s_rRNA
		complexIds = [sim_data.molecule_ids.s30_full_complex]

		# Load count data for mRNAs
		RNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_cistron_counts = RNA_counts_reader.readColumn('mRNA_cistron_counts')
		all_mRNA_cistron_indexes = {rna: i for i, rna in enumerate(RNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		rna_indexes = np.array([all_mRNA_cistron_indexes[rna] for rna in cistron_ids], int)
		rna_cistron_counts = mRNA_cistron_counts[:, rna_indexes]
		(freeProteinCounts, freeRRnaCounts, complexCounts) = read_bulk_molecule_counts(
			simOutDir, (proteinIds, rRnaIds, complexIds))
		complexCounts = complexCounts.reshape(-1, 1)

		# Load data
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
		activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

		plt.figure(figsize = (8.5, 14))
		plt.rc('font', **FONT)

		for idx in range(len(proteinIds)):
			rna_axis = plt.subplot(10, 3, idx + 1)

			sparklineAxis(rna_axis, time / 60., rna_cistron_counts[:, idx], 'left', '-', 'b')
			setAxisMaxMinY(rna_axis, rna_cistron_counts[:, idx])

			protein_axis = rna_axis.twinx()
			sparklineAxis(protein_axis, time / 60., freeProteinCounts[:, idx], 'right', '-', 'r')
			setAxisMaxMinY(protein_axis, freeProteinCounts[:, idx])

			# Component label
			rna_axis.set_title(proteinIds[idx][:-3], fontsize=8)

		for idx in range(len(rRnaIds)):
			rna_axis = plt.subplot(10, 3, idx + len(proteinIds) + 1)

			sparklineAxis(rna_axis, time / 60., freeRRnaCounts[:, idx], 'left', '-', 'b')
			setAxisMaxMinY(rna_axis, freeRRnaCounts[:, idx])

			# Component label
			rna_axis.set_title(rRnaIds[idx][:-3], fontsize=8)

		for idx in range(len(complexIds)):
			complex_axis = plt.subplot(10, 3, idx + len(proteinIds) + len(rRnaIds) + 1)

			sparklineAxis(complex_axis, time / 60., complexCounts[:, idx], 'left', '-', 'r')
			setAxisMaxMinY(complex_axis, complexCounts[:, idx])

			# Component label
			complex_axis.set_title(complexIds[idx][:-3], fontsize=8)

		# Plot number of ribosomes
		ribosome_axis = plt.subplot(10, 3, len(proteinIds) + len(rRnaIds) + len(complexIds) + 1)
		sparklineAxis(ribosome_axis, time / 60., activeRibosome, 'left', '-', 'r')
		setAxisMaxMinY(ribosome_axis, activeRibosome)
		ribosome_axis.set_title('Active ribosome', fontsize=8)

		# Save
		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
