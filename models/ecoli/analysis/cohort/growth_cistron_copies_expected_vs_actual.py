"""
Generates a scatter plot of growth-associated mRNA cistron copy numbers (genes
 that encode for ribosomal proteins and RNAPs) that are expected from the
expression levels calculated in the ParCa vs actual copies in the simulations.
Both values are normalized such that the sum of copy numbers of all plotted mRNA
species sum to one.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


# First generation (counting from zero) from which to gather RNA counts.
FIRST_GENERATION = 0

FIGSIZE = (4, 4)
BOUNDS = [1e-3, 1e-1]
NUMERICAL_ZERO = 1e-10

EXPECTED_COUNT_CONDITION = 'basal'


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		n_gens = self.ap.n_generation
		cell_paths = self.ap.get_cells(
			generation=list(range(FIRST_GENERATION, n_gens))
			)

		simOutDir = os.path.join(cell_paths[0], 'simOut')
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')

		# Calculate expected copies from sim_data
		is_mRNA = sim_data.process.transcription.cistron_data['is_mRNA']
		all_cistron_ids = sim_data.process.transcription.cistron_data['id']
		all_mRNA_cistron_ids = all_cistron_ids[is_mRNA]
		mRNA_cistron_id_to_index = {
			cistron_id: i for i, cistron_id in enumerate(all_mRNA_cistron_ids)
			}
		rprotein_cistron_ids = all_cistron_ids[
			sim_data.process.transcription.cistron_data['is_ribosomal_protein']]
		rnap_cistron_ids = all_cistron_ids[
			sim_data.process.transcription.cistron_data['is_RNAP']]
		rprotein_indexes = np.array(
			[mRNA_cistron_id_to_index[x] for x in rprotein_cistron_ids])
		rnap_indexes = np.array(
			[mRNA_cistron_id_to_index[x] for x in rnap_cistron_ids])

		expected_mRNA_counts = sim_data.process.transcription.fit_cistron_expression[EXPECTED_COUNT_CONDITION][is_mRNA]
		expected_rprotein_counts = expected_mRNA_counts[rprotein_indexes]
		expected_rnap_counts = expected_mRNA_counts[rnap_indexes]

		# Normalize counts to total sum
		t = expected_rprotein_counts.sum() + expected_rnap_counts.sum()
		expected_rprotein_counts /= t
		expected_rnap_counts /= t

		assert np.all(
			mRNA_cistron_ids == sim_data.process.transcription.cistron_data['id'][is_mRNA])

		# Initialize array for actual counts
		all_actual_rprotein_counts = np.zeros((len(rprotein_indexes), len(cell_paths)))
		all_actual_rnap_counts = np.zeros((len(rnap_indexes), len(cell_paths)))

		# Load data
		for i, sim_dir in enumerate(cell_paths):
			simOutDir = os.path.join(sim_dir, 'simOut')

			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
			mRNA_cistron_counts = mRNA_counts_reader.readColumn('mRNA_cistron_counts')

			# Get average count across all timesteps
			all_actual_rprotein_counts[:, i] = mRNA_cistron_counts[:, rprotein_indexes].mean(axis=0)
			all_actual_rnap_counts[:, i] = mRNA_cistron_counts[:, rnap_indexes].mean(axis=0)

		# Get average count across all sims
		actual_rprotein_counts = all_actual_rprotein_counts.mean(axis=1)
		actual_rnap_counts = all_actual_rnap_counts.mean(axis=1)

		# Normalize counts to total sum
		t = actual_rprotein_counts.sum() + actual_rnap_counts.sum()
		actual_rprotein_counts /= t
		actual_rnap_counts /= t

		plt.figure(figsize=FIGSIZE)

		plt.scatter(
			expected_rprotein_counts + NUMERICAL_ZERO,
			actual_rprotein_counts + NUMERICAL_ZERO,
			s=2, c='r',
			label='ribosomal proteins')
		plt.scatter(
			expected_rnap_counts + NUMERICAL_ZERO,
			actual_rnap_counts + NUMERICAL_ZERO,
			s=2, c='b',
			label='RNAP subunits')

		plt.title('Expected vs actual RNA copies')
		plt.xlabel('Expected normalized copies')
		plt.ylabel('Actual normalized copies')
		plt.xlim(BOUNDS)
		plt.ylim(BOUNDS)
		plt.xscale('log')
		plt.yscale('log')
		plt.legend(loc=1, prop={'size': 8})

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
