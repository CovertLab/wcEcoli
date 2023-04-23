"""
Produces a txt file for all, early, and late generations containing a list of
transcription units and the index of variants where those genes were
overcrowded by RNA polymerases. Here, overcrowded by RNA polymerases is
defined as the actual probability of RNA synthesis being less than the target
probability of RNA synthesis on average for at least one generation in at
least one seed for that variant index.

TODO: Filter sims that timed out
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns, \
	stacked_cell_threshold_mask

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE=9
MAX_VARIANT = 43 # do not include any variant >= this index
MAX_CELL_INDEX = 16 # do not include any generation >= this index

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
MIN_LATE_CELL_INDEX = 4

"""
Count number of sims that reach this generation (remember index 7 
corresponds to generation 8)
"""
COUNT_INDEX = 15

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Data extraction
		print("---Data Extraction---")
		overcrowded_tu_index_dicts = [{}]
		if exclude_early_gens == 1:
			overcrowded_tu_index_dicts.append({}) # Early gens
			overcrowded_tu_index_dicts.append({}) # Late gens

		variants = self.ap.get_variants()
		min_variant = min(variants)

		all_generations = {}

		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			all_cells_gens = \
				np.array([int(os.path.basename(os.path.dirname(
					cell_path))[-6:]) for cell_path in all_cells])
			all_generations[variant] = all_cells_gens

			# Get tu ids
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				rnap_reader = TableReader(
					os.path.join(simOutDir, 'RnaSynthProb'))
				rna_ids = rnap_reader.readAttribute('rnaIds')

			# RNA polymerase overcrowding
			avg_actual_rna_synth_prob = read_stacked_columns(all_cells,
				'RnaSynthProb','actual_rna_synth_prob',
				fun=lambda x: np.mean(x, axis=0))
			avg_target_rna_synth_prob = read_stacked_columns(all_cells,
				'RnaSynthProb','target_rna_synth_prob',
				fun=lambda x: np.mean(x, axis=0))

			avg_overcrowded_tu_indexes = np.where(
				sum(avg_actual_rna_synth_prob <
					avg_target_rna_synth_prob) > 0)[0]

			if exclude_early_gens == 1:
				early_cell_mask = all_generations[variant] < MIN_LATE_CELL_INDEX
				if len(early_cell_mask) == 1:
					early_cell_mask = early_cell_mask[0]
				avg_overcrowded_tu_indexes_early = np.where(
					sum((avg_actual_rna_synth_prob <
						 avg_target_rna_synth_prob)[early_cell_mask]) > 0)[0]

				late_cell_mask = all_generations[variant] >= MIN_LATE_CELL_INDEX
				if len(late_cell_mask) == 1:
					late_cell_mask = late_cell_mask[0]
				avg_overcrowded_tu_indexes_late = np.where(
					sum((avg_actual_rna_synth_prob <
						 avg_target_rna_synth_prob)[late_cell_mask]) > 0)[0]

			# Record that these genes were overcrowded this variant
			for tu_index in avg_overcrowded_tu_indexes:
				if tu_index not in overcrowded_tu_index_dicts[0]:
					overcrowded_tu_index_dicts[0][tu_index] = []
				overcrowded_tu_index_dicts[0][tu_index].append(
					variant)

			if exclude_early_gens == 1:
				for tu_index in avg_overcrowded_tu_indexes_early:
					if tu_index not in overcrowded_tu_index_dicts[1]:
						overcrowded_tu_index_dicts[1][tu_index] = []
					overcrowded_tu_index_dicts[1][tu_index].append(
						variant)

				for tu_index in avg_overcrowded_tu_indexes_late:
					if tu_index not in overcrowded_tu_index_dicts[2]:
						overcrowded_tu_index_dicts[2][tu_index] = []
					overcrowded_tu_index_dicts[2][tu_index].append(
						variant)

		# Determine the ids corresponding to these tus
		overcrowded_tu_indexes = overcrowded_tu_index_dicts[0].keys()
		overcrowded_gene_indexes_to_ids_dict = {tu_index:
			rna_ids[tu_index] for tu_index in overcrowded_tu_indexes}

		# Output
		print("---Writing Output---")
		output_descr = ["_all_gens"]
		if exclude_early_gens == 1:
			output_descr += ["_early_gens", "_late_gens"]

		for j in range(len(overcrowded_tu_index_dicts)):
			with open(os.path.join(plotOutDir, plotOutFileName +
				output_descr[j] + ".txt"), 'w') as f:
				for tu_id in overcrowded_tu_index_dicts[j]:
					overcrowded_variant_string = ", ".join(map(str,
						overcrowded_tu_index_dicts[j][tu_id]))
					f.write("" + overcrowded_gene_indexes_to_ids_dict[
						tu_id] + ": " + overcrowded_variant_string + "\n")

if __name__ == "__main__":
	Plot().cli()
