"""
Plots dynamics of variables related to the RNA maturation process.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

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

		if mature_rnas_generated.ndim == 2:
			total_mature_rnas_generated = mature_rnas_generated.sum(axis=1)
		else:
			total_mature_rnas_generated = mature_rnas_generated

		plt.figure(figsize=(6, 6))
		ax = plt.subplot(3, 1, 1)
		ax.plot(time_min, total_maturation_events)
		ax.set_ylabel('Number of\nmaturation events')

		ax = plt.subplot(3, 1, 2)
		ax.plot(time_min, total_mature_rnas_generated)
		ax.set_ylabel('Number of mature RNAs\ngenerated')

		ax = plt.subplot(3, 1, 3)
		ax.plot(time_min, total_degraded_ntps)
		ax.set_ylabel('Number of NTPs\ndegraded')
		ax.set_xlabel('Time (min)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
