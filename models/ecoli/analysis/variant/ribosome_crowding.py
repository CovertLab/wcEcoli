"""
Produces a txt file containing a list of genes and the index of variants
where those genes were overcrowded by ribosomes. Here, overcrowded by
ribosomes is defined as the actual probability of translation per
transcript being less than the target probability of translation per
transcript on average for at least one generation in at least one seed for
that variant index.

TODO: Filter early vs late gens, sims that times out
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure,\
	read_stacked_columns, stacked_cell_threshold_mask
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS\
	as COLORS, labeled_indexable_hist, labeled_indexable_scatter

FONT_SIZE=9
MAX_VARIANT = 43 # do not include any variant >= this index

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Data extraction
		print("---Data Extraction---")
		overcrowded_monomer_index_dict = {}

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
			avg_overcrowded_monomer_indexes = np.where(
				sum(avg_actual_prob_translation_per_transcript <
				avg_target_prob_translation_per_transcript) > 0)[0]

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
		monomer_to_mRNA_id_dict = dict(zip(monomer_sim_data['id'],
										   monomer_sim_data[
											   'cistron_id']))
		mRNA_to_gene_id_dict = dict(zip(mRNA_sim_data['id'],
										mRNA_sim_data['gene_id']))

		overcrowded_monomer_indexes = overcrowded_monomer_index_dict.keys()
		overcrowded_gene_indexes_to_ids_dict = {monomer_index:
			mRNA_to_gene_id_dict.get(monomer_to_mRNA_id_dict.get(
				monomer_ids[monomer_index])) for monomer_index in
			overcrowded_monomer_indexes}

		# Output
		print("---Writing Output---")

		print(os.path.join(plotOutDir, plotOutFileName + ".txt"))

		with open(os.path.join(plotOutDir, plotOutFileName + ".txt"),'w') as f:
			for monomer_id in overcrowded_monomer_index_dict:
				overcrowded_variant_string = ", ".join(map(str,
					overcrowded_monomer_index_dict[
					monomer_id]))
				f.write("" + overcrowded_gene_indexes_to_ids_dict[
					monomer_id] + ": " + overcrowded_variant_string + "\n")

if __name__ == "__main__":
	Plot().cli()
