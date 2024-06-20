"""
Plots and compares the various components of the model that contribute to the
regulation of rRNA operons.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns)
from wholecell.io.tablereader import TableReader

TU_ID_TO_RRNA_OPERON_ID = {
	'TU0-1181[c]': 'rrnA',
	'TU0-1182[c]': 'rrnB',
	'TU0-1183[c]': 'rrnC',
	'TU0-1191[c]': 'rrnD',
	'TU0-1186[c]': 'rrnE',
	'TU0-1187[c]': 'rrnG',
	'TU0-1189[c]': 'rrnH',
	}


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		rna_data = sim_data.process.transcription.rna_data
		rna_ids_orig = rna_data['id']
		rna_coordinates = rna_data['replication_coordinate']

		rRNA_ids = rna_data['id'][rna_data['is_rRNA']]

		# If the list of rRNA IDs does not match the rRNA TU IDs, skip analysis
		# (likely using rRNA operon structure variants)
		if set(rRNA_ids) != set(TU_ID_TO_RRNA_OPERON_ID.keys()):
			print(
				'Skipping analysis - this analysis should be run for simulations with WT rRNA operon structures.')
			return

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Convert coordiantes to theta values (assuming oriC is at np.pi/2)
		rna_thetas = -np.pi*(
			rna_coordinates
			/ replichore_lengths[(rna_coordinates < 0).astype(int)]
			)

		# Reorder RNA IDs in order of theta
		theta_argsort = rna_thetas.argsort()
		rna_thetas_sorted = np.sort(rna_thetas)
		rna_ids = rna_ids_orig[theta_argsort]

		variants = self.ap.get_variants()
		variant_index_to_p_per_copy = {}
		variant_index_to_n_copy_numbers = {}
		variant_index_to_p = {}

		for variant in variants:
			cell_paths = self.ap.get_cells(
				variant=[variant], only_successful=True)

			# Confirm RnaSynthProb uses same RNA ID array
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			rna_synth_prob_rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')

			assert np.all(rna_synth_prob_rna_ids == rna_ids_orig)

			# Read values
			all_rna_synth_probs = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob')
			all_copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'promoter_copy_number')
			rna_synth_prob = all_rna_synth_probs.mean(axis=0)
			copy_numbers = all_copy_numbers.mean(axis=0)
			rna_synth_prob_per_copy = (
				all_rna_synth_probs / all_copy_numbers).mean(axis=0)

			# Normalize and reorder values
			rna_synth_prob /= rna_synth_prob.sum()
			rna_synth_prob_per_copy /= rna_synth_prob_per_copy.sum()

			rna_synth_prob = rna_synth_prob[theta_argsort]
			rna_synth_prob_per_copy = rna_synth_prob_per_copy[theta_argsort]
			copy_numbers = copy_numbers[theta_argsort]

			# Write to dictionary
			variant_index_to_p_per_copy[variant] = rna_synth_prob_per_copy
			variant_index_to_n_copy_numbers[variant] = copy_numbers
			variant_index_to_p[variant] = rna_synth_prob

		# Plot
		n_variants = len(variants)
		fig, axes = plt.subplots(
			n_variants, 3, subplot_kw={'projection': 'polar'},
			figsize=(15, 5*n_variants))

		def setup_ax(ax):
			ax.set_theta_offset(np.pi / 2)  # Put 0 deg at the top instead of the right side
			ax.grid(False)
			ax.set_yticklabels([])
			ax.axis("off")

			# Plot chromosome outline
			ax.plot(np.linspace(0, 2 * np.pi, 1000), np.repeat(1, 1000), color='k')
			ax.vlines(0, 0.98, 1.02, color='k')
			ax.vlines(np.pi, 0.98, 1.02, color='k')

			# Add labels for oriC and terC
			ax.text(0, 0.87, "oriC", ha="center", va="center")
			ax.text(np.pi, 0.87, "terC", ha="center", va="center")

		for variant_index in variant_index_to_p_per_copy.keys():
			setup_ax(axes[variant_index][0])

			# Plot locations of rRNA transcripts
			axes[variant_index][0].plot(
				rna_thetas_sorted,
				1 + 10 * variant_index_to_p_per_copy[variant_index], color='C0')
			plt.tight_layout()

			setup_ax(axes[variant_index][1])

			# Plot locations of rRNA transcripts
			axes[variant_index][1].plot(
				rna_thetas_sorted,
				1 + 0.05 * variant_index_to_n_copy_numbers[variant_index], color='C1')
			plt.tight_layout()

			setup_ax(axes[variant_index][2])

			# Plot locations of rRNA transcripts
			axes[variant_index][2].plot(
				rna_thetas_sorted,
				1 + 10 * variant_index_to_p[variant_index], color='k')
			plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
