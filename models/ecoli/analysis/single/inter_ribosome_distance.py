"""
Computes the average distance between translating ribosomes for each
mRNA, and compares the distance to the known size of the
molecular footprint of ribosomes.
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

		# Get ribosome footprint size
		ribosome_footprint_size = sim_data.process.translation.active_ribosome_footprint_size.asNumber(units.nt)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		ribosome_data_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		monomer_counts_reader = TableReader(os.path.join(simOutDir, 'MonomerCounts'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = (units.s)*(main_reader.readColumn('time') - initial_time)
		ribosome_init_event_per_monomer = ribosome_data_reader.readColumn('ribosome_init_event_per_monomer')
		cistron_counts = mRNA_counts_reader.readColumn('mRNA_cistron_counts')
		monomer_ids = monomer_counts_reader.readAttribute('monomerIds')

		# Map monomer ids to cistron indices
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		monomer_to_cistron_id_dict = dict(
			zip(monomer_sim_data['id'], monomer_sim_data['cistron_id']))
		cistron_idx_dict = {rna: i for i, rna in
			enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		cistron_indices_for_monomers = [
			cistron_idx_dict[monomer_to_cistron_id_dict[monomer_id]]
			for monomer_id in monomer_ids]
		cistron_counts_for_monomers = cistron_counts[:, cistron_indices_for_monomers]

		# Get ribosome elongation rate for this sim condition
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		ribosome_elong_rate = sim_data.process.translation.ribosomeElongationRateDict[nutrients]

		# Calculate the total number of initiation events that happen to each
		# transcript per mRNA throughout the cell cycle
		n_total_ribosome_init_events_per_mRNA = np.divide(
			ribosome_init_event_per_monomer[1:, :].astype(np.float64),
			cistron_counts_for_monomers[1:, :]).sum(axis=0)

		# Divide by length of cell cycle to get average initiation rate
		avg_init_rate = (1. / time[-1]) * n_total_ribosome_init_events_per_mRNA

		# Divide elongation rate with initiation rate to get average distance
		# between RNAPs in nucleotides
		avg_inter_ribosome_distance = (ribosome_elong_rate / avg_init_rate).asNumber(units.nt)

		# Sort from shortest to longest
		sorted_index = avg_inter_ribosome_distance.argsort()
		sorted_monomer_ids = [monomer_ids[i] for i in sorted_index]
		avg_inter_ribosome_distance.sort()

		# Mark genes with RNAPs that are too close to each other
		n_too_close = (avg_inter_ribosome_distance[:SAMPLE_SIZE] < ribosome_footprint_size).sum()
		bar_colors = ["r"]*n_too_close + ["b"]*(SAMPLE_SIZE - n_too_close)

		# Plot the first n genes with shortest distances
		plt.figure(figsize=(8, 6))
		plt.barh(np.arange(SAMPLE_SIZE), avg_inter_ribosome_distance[:SAMPLE_SIZE],
			tick_label=sorted_monomer_ids[:SAMPLE_SIZE],
			color=bar_colors)
		plt.xlabel("Average distance between ribosomes (nt)")
		plt.axvline(ribosome_footprint_size, linestyle='--', color='k')

		# Add values to each bar
		for i, v in enumerate(avg_inter_ribosome_distance[:SAMPLE_SIZE]):
			if np.isfinite(v):
				plt.text(v - 1, i, "{0:.1f}".format(v), color='white', fontsize=5,
					horizontalalignment='right', verticalalignment='center')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == '__main__':
	Plot().cli()
