"""
Produces a txt file for all, early, and late generations containing a list of
genes and the index of variants where those genes were overcrowded by
ribosomes. Here, overcrowded by ribosomes is defined as the actual
probability of translation per transcript being less than the target
probability of translation per transcript on average for at least one
generation in at least one seed for that variant index.

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
		overcrowded_monomer_index_dicts = [{}]
		if exclude_early_gens == 1:
			overcrowded_monomer_index_dicts.append({}) # Early gens
			overcrowded_monomer_index_dicts.append({}) # Late gens

		variants = self.ap.get_variants()
		min_variant = min(variants)

		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant],
										  only_successful=True)
			if len(all_cells) == 0:
				continue

			# Get monomer ids
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				ribosome_reader = TableReader(
					os.path.join(simOutDir, 'RibosomeData'))
				monomer_ids = ribosome_reader.readAttribute('monomerIds')

			avg_actual_prob_translation_per_transcript = read_stacked_columns(all_cells,
				'RibosomeData', 'actual_prob_translation_per_transcript',
				fun=lambda x: np.mean(x, axis = 0))
			avg_target_prob_translation_per_transcript = read_stacked_columns(
				all_cells,
				'RibosomeData', 'target_prob_translation_per_transcript',
				fun=lambda x: np.mean(x, axis = 0))

			# Get indexes of proteins corresponding to mRNAs that on
			# average were overcrowded in any generation for any seed
			overcrowded_gens, avg_overcrowded_monomer_indexes = np.where(
				avg_actual_prob_translation_per_transcript <
					avg_target_prob_translation_per_transcript)

			# Record that these genes were overcrowded this variant
			for i in range(len(avg_overcrowded_monomer_indexes)):
				gen = overcrowded_gens[i]
				monomer_index = avg_overcrowded_monomer_indexes[i]
				if monomer_index not in overcrowded_monomer_index_dicts[0]:
					for j in range(len(overcrowded_monomer_index_dicts)):
						overcrowded_monomer_index_dicts[j][monomer_index] = []

				for j in range(len(overcrowded_monomer_index_dicts)):
					if j == 1 and gen >= MIN_LATE_CELL_INDEX:
						continue
					if j == 2 and gen < MIN_LATE_CELL_INDEX:
						continue
					if variant not in overcrowded_monomer_index_dicts[j][
						monomer_index]:
						overcrowded_monomer_index_dicts[j][monomer_index].append(
							variant)

		# Determine the gene ids corresponding to these proteins
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		monomer_to_mRNA_id_dict = dict(zip(monomer_sim_data['id'],
										   monomer_sim_data[
											   'cistron_id']))
		mRNA_to_gene_id_dict = dict(zip(mRNA_sim_data['id'],
										mRNA_sim_data['gene_id']))

		overcrowded_monomer_indexes = overcrowded_monomer_index_dicts[0].keys()
		overcrowded_gene_indexes_to_ids_dict = {monomer_index:
			mRNA_to_gene_id_dict.get(monomer_to_mRNA_id_dict.get(
				monomer_ids[monomer_index])) for monomer_index in
			overcrowded_monomer_indexes}

		# Output
		print("---Writing Output---")
		output_descr = ["_all_gens"]
		if exclude_early_gens == 1:
			output_descr += ["_early_gens", "_late_gens"]

		for j in range(len(overcrowded_monomer_index_dicts)):
			with open(os.path.join(plotOutDir, plotOutFileName +
				output_descr[j] + ".txt"), 'w') as f:
				for monomer_id in overcrowded_monomer_index_dicts[j]:
					overcrowded_variant_string = ", ".join(map(str,
						overcrowded_monomer_index_dicts[j][monomer_id]))
					f.write("" + overcrowded_gene_indexes_to_ids_dict[
						monomer_id] + ": " + overcrowded_variant_string + "\n")

if __name__ == "__main__":
	Plot().cli()
