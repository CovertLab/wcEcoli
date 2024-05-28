"""
Compares the averaged densities of RNAPs on genes encoding for rRNAs across
variants.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


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

		# Load IDs of rRNA transcription units
		rna_data = sim_data.process.transcription.rna_data
		rRNA_ids = rna_data['id'][rna_data['is_rRNA']]

		# If the list of rRNA IDs does not match the rRNA TU IDs, skip analysis
		# (likely using rRNA operon structure variants)
		if set(rRNA_ids) != set(TU_ID_TO_RRNA_OPERON_ID.keys()):
			print(
				'Skipping analysis - this analysis should be run for simulations with WT rRNA operon structures.')
			return

		rnap_footprint_size = sim_data.process.transcription.active_rnap_footprint_size.asNumber(
			units.nt)

		variants = self.ap.get_variants()
		variant_index_to_distance_each_rrna = {}
		variant_index_to_avg_distance = {}

		for variant in variants:
			with open(self.ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)

			# Get RNAP elongation rate used in simulations
			nutrients = variant_sim_data.conditions[variant_sim_data.condition]["nutrients"]
			rnap_elong_rate = variant_sim_data.process.transcription.rnaPolymeraseElongationRateDict[
				nutrients].asNumber(units.nt / units.s)

			# Get cell paths
			cell_paths = self.ap.get_cells(
				variant=[variant], only_successful=True)

			# Get rna ID attributes from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rnap_data_reader = TableReader(
				os.path.join(sim_out_dir, 'RnapData'))
			rnap_data_rna_ids = rnap_data_reader.readAttribute('rnaIds')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			rna_synth_prob_rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')

			# Get indexes of rRNA operons in both listeners
			rnap_data_rrna_indexes = np.array([
				rnap_data_rna_ids.index(key)
				for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])
			rna_synth_prob_rrna_indexes = np.array([
				rna_synth_prob_rna_ids.index(key)
				for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])

			# Get number of rRNA initiation events and copy numbers
			initiation_events = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True, remove_first=True,
				fun=lambda x: x[:, rnap_data_rrna_indexes]
				)
			copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'promoter_copy_number',
				ignore_exception=True, remove_first=True,
				fun=lambda x: x[:, rna_synth_prob_rrna_indexes])

			# Calculate total number of initiation events per copy
			initiation_events_per_copy_each_rrna = (
				initiation_events / copy_numbers
				).sum(axis=0)

			# Get mask for rRNAs with zero initiation events throughout all
			# sims - these rRNAs are assumed to be nonfunctional (most likely
			# knocked out)
			nonfunctional_rrna_mask = (
				initiation_events_per_copy_each_rrna == 0)

			# Calculate initiation events per copy amongst all functional rRNAs
			initiation_events_per_copy_all_functional_rrnas = (
				initiation_events[:, ~nonfunctional_rrna_mask].sum(axis=1)
				/ copy_numbers[:, ~nonfunctional_rrna_mask].sum(axis=1)
				).sum()

			# Get total simulation time
			generation_lengths = read_stacked_columns(
				cell_paths, 'Main', 'time',
				ignore_exception=True,
				fun=lambda x: x[-1] - x[0])
			total_simulation_time = generation_lengths.sum()

			# Calculate average distance between adjacent RNAPs for each operon
			# and for all operons
			with np.errstate(divide='ignore'):
				avg_rrna_distance_each_rrna = (
					rnap_elong_rate * total_simulation_time
					/ initiation_events_per_copy_each_rrna)
			avg_rrna_distance_all_rrnas = (
				rnap_elong_rate * total_simulation_time
				/ initiation_events_per_copy_all_functional_rrnas)

			variant_index_to_distance_each_rrna[variant] = avg_rrna_distance_each_rrna
			variant_index_to_avg_distance[variant] = avg_rrna_distance_all_rrnas

		fig = plt.figure(figsize=(6, 4))

		# Plot distances for each rRNA operon and the average
		ax = fig.add_subplot(1, 1, 1)
		for variant_index, distances in variant_index_to_distance_each_rrna.items():
			ax.scatter(
				np.arange(8),
				list(distances) + [variant_index_to_avg_distance[variant_index]],
				marker='x', label=f'Variant {variant_index}')

		# Draw dotted line at RNAP's molecular footprint
		ax.axhline(rnap_footprint_size, ls='--', lw=2, color='k')

		ax.set_xticks(np.arange(8))
		ax.set_xticklabels([v for v in TU_ID_TO_RRNA_OPERON_ID.values()] + ['avg'])
		ax.set_xlim([-0.8, 7.8])
		ax.set_ylim([0, 600])
		ax.set_ylabel('Avg distance between RNAPs (nt)')
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
