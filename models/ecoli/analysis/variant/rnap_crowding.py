"""
Produces a txt file for all, early, and late generations containing a list of
transcription units and the index of variants where those genes were
overcrowded by RNA polymerases. Here, overcrowded by RNA polymerases is
defined as the actual probability of RNA synthesis being less than the target
probability of RNA synthesis on average for at least one generation in at
least one seed for that variant index.
"""

import numpy as np
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import read_stacked_columns

FONT_SIZE=9
MAX_VARIANT = 43 # do not include any variant >= this index

"""
Plot data from generations [MIN_CELL_INDEX, MAX_CELL_INDEX)
Note that early generations may not be representative of dynamics 
due to how they are initialized
"""
MIN_CELL_INDEX = 4
MAX_CELL_INDEX = 16

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# Data extraction
		print("---Data Extraction---")
		overcrowded_tu_index_dict = {}

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

			all_cells_gens = np.array([
				int(os.path.basename(os.path.dirname(
				cell_path))[-6:]) for cell_path in all_cells])
			all_generations[variant] = all_cells_gens
			cell_mask = np.logical_and(
				(all_generations[variant] >=MIN_CELL_INDEX),
				(all_generations[variant] < MAX_CELL_INDEX))
			if len(cell_mask) == 1:
				cell_mask = cell_mask.reshape(1)
			if sum(cell_mask) < 1:
				continue

			# Get tu ids
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				rnap_reader = TableReader(
					os.path.join(simOutDir, 'RnaSynthProb'))
				rna_ids = rnap_reader.readAttribute('rnaIds')

			# RNA polymerase overcrowding
			avg_actual_rna_synth_prob = read_stacked_columns(
				all_cells, 'RnaSynthProb', 'actual_rna_synth_prob',
				fun=lambda x: np.mean(x, axis=0))
			avg_target_rna_synth_prob = read_stacked_columns(
				all_cells, 'RnaSynthProb', 'target_rna_synth_prob',
				fun=lambda x: np.mean(x, axis=0))

			avg_overcrowded_tu_indexes = np.where(
				sum((avg_actual_rna_synth_prob <
					avg_target_rna_synth_prob)[cell_mask]) > 0)[0]

			# Record that these genes were overcrowded this variant
			for tu_index in avg_overcrowded_tu_indexes:
				if tu_index not in overcrowded_tu_index_dict:
					overcrowded_tu_index_dict[tu_index] = []
				overcrowded_tu_index_dict[tu_index].append(variant)

		# Determine the ids corresponding to these tus
		overcrowded_tu_indexes = overcrowded_tu_index_dict.keys()
		overcrowded_gene_indexes_to_ids_dict = {
			tu_index: rna_ids[tu_index] for tu_index in overcrowded_tu_indexes}

		# Output
		print("---Writing Output---")
		plot_suffix = "_gens_" + str(MIN_CELL_INDEX) + "_through_" + str(MAX_CELL_INDEX)

		with open(os.path.join(plotOutDir, plotOutFileName +
			plot_suffix + ".txt"), 'w') as f:
			for tu_id in overcrowded_tu_index_dict:
				overcrowded_variant_string = ", ".join(map(str,
					overcrowded_tu_index_dict[tu_id]))
				f.write(f"{overcrowded_gene_indexes_to_ids_dict[tu_id]}:"
					f" {overcrowded_variant_string}\n")

if __name__ == "__main__":
	Plot().cli()
