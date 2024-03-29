"""
Plot protein monomer counts
"""

import os
import pickle

import numpy as np
from scipy.stats import pearsonr
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get the names of proteins from the KB

		sim_data = self.read_pickle_file(simDataFile)

		monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		avgCounts = monomerCounts.readColumn("monomerCounts").mean(axis=0)

		relativeCounts = avgCounts / avgCounts.sum()

		protein_loss_rates = (
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
			+ sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)
			)
		expectedCountsArbitrary = normalize(
			sim_data.process.transcription.cistron_tu_mapping_matrix.dot(
				sim_data.process.transcription.rna_expression[sim_data.condition])[
				sim_data.relation.cistron_to_monomer_mapping]
			* sim_data.process.translation.translation_efficiencies_by_monomer
			/ protein_loss_rates
			)

		expectedCountsRelative = expectedCountsArbitrary / expectedCountsArbitrary.sum()

		plt.figure(figsize = (8.5, 11))

		maxLine = 1.1 * max(np.log10(expectedCountsRelative.max() + 1), np.log10(relativeCounts.max() + 1))
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1), 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("log10(Expected protein distribution (from fitting))")
		plt.ylabel("log10(Actual protein distribution (average over life cycle))")
		plt.title("PCC (of log values): %0.2f" % pearsonr(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1))[0])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
