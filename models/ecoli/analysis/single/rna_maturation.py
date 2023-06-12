"""
Plots dynamics of variables related to the RNA maturation process.
"""

import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rna_maturation_reader =	TableReader(os.path.join(simOutDir, 'RnaMaturationListener'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		time_min = time / 60
		total_maturation_events = rna_maturation_reader.readColumn('total_maturation_events')
		total_degraded_ntps = rna_maturation_reader.readColumn('total_degraded_ntps')
		mature_rnas_generated = rna_maturation_reader.readColumn(
			'mature_rnas_generated')
		maturation_enzyme_counts = rna_maturation_reader.readColumn(
			'maturation_enzyme_counts')
		enzyme_ids = rna_maturation_reader.readAttribute('enzyme_ids')

		if mature_rnas_generated.ndim == 2:
			total_mature_rnas_generated = mature_rnas_generated.sum(axis=1)
		else:
			total_mature_rnas_generated = mature_rnas_generated

		plt.figure(figsize=(6, 8))
		ax0 = plt.subplot(4, 1, 1)
		ax0.plot(time_min, total_maturation_events)
		ax0.set_ylabel('Number of\nmaturation events')

		ax1 = plt.subplot(4, 1, 2)
		ax1.plot(time_min, total_mature_rnas_generated)
		ax1.set_ylabel('Number of mature RNAs\ngenerated')

		ax2 = plt.subplot(4, 1, 3)
		ax2.plot(time_min, total_degraded_ntps)
		ax2.set_ylabel('Number of NTPs\ndegraded')

		ax3 = plt.subplot(4, 1, 4)
		ax3.plot(time_min, maturation_enzyme_counts)
		ax3.legend(enzyme_ids)
		ax3.set_ylabel('RNA maturation\nenzyme counts')
		ax3.set_xlabel('Time (min)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
