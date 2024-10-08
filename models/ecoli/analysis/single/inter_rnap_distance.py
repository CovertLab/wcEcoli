"""
Computes the average distance between transcribing RNA polymerases for each
gene (transcription unit), and compares the distance to the known size of the
molecular footprint of RNAP.
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


SAMPLE_SIZE = 25  # Number of genes to plot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Get RNAP footprint size
		RNAP_footprint_size = sim_data.process.transcription.active_rnap_footprint_size.asNumber(units.nt)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rnap_data_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
		rna_synth_prob_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = (units.s)*(main_reader.readColumn('time') - initial_time)
		n_rna_init_events = rnap_data_reader.readColumn('rnaInitEvent')
		promoter_copy_numbers = rna_synth_prob_reader.readColumn('promoter_copy_number')
		TU_ids = rna_synth_prob_reader.readAttribute('rnaIds')

		# Get RNAP elongation rate for sim condition
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		rnap_elong_rate = sim_data.process.transcription.rnaPolymeraseElongationRateDict[nutrients]

		# Calculate the total number of initiation events that happen to each
		# gene per gene copy throughout the cell cycle
		sum_n_rna_init_events = n_rna_init_events[1:, :].astype(np.float64).sum(axis=0)
		avg_promoter_copy_number = np.mean(promoter_copy_numbers[1:, :], axis=0)
		n_total_rna_init_events_per_copy = sum_n_rna_init_events / avg_promoter_copy_number

		# Divide by length of cell cycle to get average initiation rate
		avg_init_rate = (1./time[-1])*n_total_rna_init_events_per_copy

		# Divide elongation rate with initiation rate to get average distance
		# between RNAPs in nucleotides
		avg_inter_rnap_distance = (rnap_elong_rate/avg_init_rate).asNumber(units.nt)

		# Sort from shortest to longest
		sorted_index = avg_inter_rnap_distance.argsort()
		sorted_TU_ids = [TU_ids[i][:-3] for i in sorted_index]  # [c] stripped
		avg_inter_rnap_distance.sort()

		# Mark genes with RNAPs that are too close to each other
		n_too_close = (avg_inter_rnap_distance[:SAMPLE_SIZE] < RNAP_footprint_size).sum()
		bar_colors = ["r"]*n_too_close + ["b"]*(SAMPLE_SIZE - n_too_close)

		# Plot the first n genes with shortest distances
		plt.figure(figsize=(4, 6))
		plt.barh(np.arange(SAMPLE_SIZE), avg_inter_rnap_distance[:SAMPLE_SIZE],
			tick_label=sorted_TU_ids[:SAMPLE_SIZE],
			color=bar_colors)
		plt.xlabel("Average distance between RNAPs (nt)")
		plt.axvline(RNAP_footprint_size, linestyle='--', color='k')

		# Add values to each bar
		for i, v in enumerate(avg_inter_rnap_distance[:SAMPLE_SIZE]):
			if np.isfinite(v):
				plt.text(v - 1, i, "{0:.1f}".format(v), color='white', fontsize=5,
					horizontalalignment='right', verticalalignment='center')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
