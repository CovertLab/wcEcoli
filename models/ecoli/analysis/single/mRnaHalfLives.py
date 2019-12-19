"""
Plot mRNA half lives (observed vs. actual)

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/30/2015
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get the expected degradation rates from KB
		sim_data = cPickle.load(open(simDataFile))
		mRNA_ids = sim_data.process.transcription.rnaData['id']
		isMRna = sim_data.process.transcription.rnaData["isMRna"]
		expected_degradation_rate = sim_data.process.transcription.rnaData['degRate'][isMRna].asNumber()

		# Get length of simulation
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		sim_length = main_reader.readColumn('time')[-1]

		# Read counts of mRNAs
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		mRNA_counts_mean = mRNA_counts.mean(axis=0)

		# Read number of degradation events
		rna_degradation_reader = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		n_RNA_degraded = rna_degradation_reader.readColumn('countRnaDegraded')
		n_mRNA_degraded_total = n_RNA_degraded.sum(axis = 0)[isMRna]

		observed_rates = []
		expected_rates = []
		ids = []

		# Loop through each gene
		for i in range(len(mRNA_counts_mean)):
			# Only account for mRNAs with mean counts above 3
			if mRNA_counts_mean[i] > 3:
				observed_rates.append(
					(n_mRNA_degraded_total[i]/sim_length) / mRNA_counts_mean[i]
					)
				expected_rates.append(expected_degradation_rate[i])
				ids.append(mRNA_ids[i])

		if len(expected_rates) == 0:
			print("Skipping analysis - RNA counts not sufficient for analysis")
			return

		observed_rates = np.array(observed_rates)
		expected_rates = np.array(expected_rates)

		plt.figure(figsize = (8, 8))
		maxLine = 1.1*np.max(np.concatenate((observed_rates, expected_rates)))

		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expected_rates, observed_rates,
			'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("Expected RNA decay")
		plt.ylabel("Observed RNA decay")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
