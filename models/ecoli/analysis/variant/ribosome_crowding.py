"""
Produces a txt file for selected generations containing a list of
genes and the index of variants where those genes were overcrowded by
ribosomes. Here, overcrowded by ribosomes is defined as the actual
probability of translation per transcript being less than the target
probability of translation per transcript on average for at least one
generation in at least one seed for that variant index.
"""

import numpy as np
import pickle
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
		overcrowded_monomer_index_dict = {}

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

			# Get monomer ids
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				ribosome_reader = TableReader(
					os.path.join(simOutDir, 'RibosomeData'))
				monomer_ids = ribosome_reader.readAttribute('monomerIds')

			avg_actual_prob_translation_per_transcript = read_stacked_columns(
				all_cells, 'RibosomeData', 'actual_prob_translation_per_transcript',
				fun=lambda x: np.mean(x, axis = 0))
			avg_target_prob_translation_per_transcript = read_stacked_columns(
				all_cells, 'RibosomeData', 'target_prob_translation_per_transcript',
				fun=lambda x: np.mean(x, axis = 0))

			avg_overcrowded_monomer_indexes = np.where(
				sum((avg_actual_prob_translation_per_transcript <
					avg_target_prob_translation_per_transcript)[cell_mask]) > 0)[0]

			# Record that these genes were overcrowded this variant
			for monomer_index in avg_overcrowded_monomer_indexes:
				if monomer_index not in overcrowded_monomer_index_dict:
					overcrowded_monomer_index_dict[monomer_index] = []
				overcrowded_monomer_index_dict[monomer_index].append(variant)

		# Determine the gene ids corresponding to these proteins
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		monomer_to_mRNA_id_dict = dict(
			zip(monomer_sim_data['id'], monomer_sim_data['cistron_id']))
		mRNA_to_gene_id_dict = dict(zip(mRNA_sim_data['id'], mRNA_sim_data['gene_id']))

		overcrowded_monomer_indexes = overcrowded_monomer_index_dict.keys()
		overcrowded_gene_indexes_to_ids_dict = {
			monomer_index: mRNA_to_gene_id_dict.get(monomer_to_mRNA_id_dict.get(
			monomer_ids[monomer_index])) for monomer_index in
			overcrowded_monomer_indexes}

		# Output
		print("---Writing Output---")
		plot_suffix = "_gens_" + str(MIN_CELL_INDEX) + "_through_" + str(MAX_CELL_INDEX)

		with open(os.path.join(plotOutDir, plotOutFileName +
			plot_suffix + ".txt"), 'w') as f:
			for monomer_id in overcrowded_monomer_index_dict:
				overcrowded_variant_string = ", ".join(
					map(str,overcrowded_monomer_index_dict[monomer_id]))
				f.write(f"{overcrowded_gene_indexes_to_ids_dict[monomer_id]}:"
					f" {overcrowded_variant_string}\n")

if __name__ == "__main__":
	Plot().cli()
