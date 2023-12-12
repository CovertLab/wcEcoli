"""
Plot tRNA counts
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get the names of rnas from the KB
		sim_data = self.read_pickle_file(simDataFile)
		transcription = sim_data.process.transcription
		uncharged_trna_ids = transcription.uncharged_trna_names
		charged_trna_ids = transcription.charged_trna_names
		cistron_is_tRNA = transcription.cistron_data['is_tRNA']

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mol_indices = {mol: i for i, mol in enumerate(moleculeIds)}

		uncharged_indices = np.array([mol_indices[moleculeId] for moleculeId in uncharged_trna_ids], int)
		charged_indices = np.array([mol_indices[moleculeId] for moleculeId in charged_trna_ids], int)

		bulk_counts = bulkMolecules.readColumn("counts")
		rna_counts = bulk_counts[:, uncharged_indices] + bulk_counts[:, charged_indices]

		plt.figure(figsize = (8.5, 11))

		counts = rna_counts[-1, :]
		expectedCountsArbitrary = transcription.cistron_tu_mapping_matrix.dot(
			transcription.rna_expression[sim_data.condition])[cistron_is_tRNA]
		expectedCounts = expectedCountsArbitrary/expectedCountsArbitrary.sum() * counts.sum()

		maxLine = 1.1 * max(expectedCounts.max(), counts.max())
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("Expected tRNA count (scaled to total)")
		plt.ylabel("Actual tRNA count (at final time step)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
