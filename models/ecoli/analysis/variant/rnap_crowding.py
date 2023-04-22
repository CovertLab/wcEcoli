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

		for variant in variants:

			if variant >= MAX_VARIANT:
				continue

			print("Variant: ",variant)
			all_cells = self.ap.get_cells(variant=[variant],
										  only_successful=True)
			if len(all_cells) == 0:
				continue

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

			# Get indexes of proteins corresponding to mRNAs that on
			# average were overcrowded in any generation for any seed
			overcrowded_gens, avg_overcrowded_tu_indexes = np.where(
				avg_actual_rna_synth_prob <
					avg_target_rna_synth_prob)

			# Record that these genes were overcrowded this variant
			for i in range(len(avg_overcrowded_tu_indexes)):
				gen = overcrowded_gens[i]
				tu_index = avg_overcrowded_tu_indexes[i]
				if tu_index not in overcrowded_tu_index_dicts[0]:
					for j in range(len(overcrowded_tu_index_dicts)):
						overcrowded_tu_index_dicts[j][tu_index] = []

				for j in range(len(overcrowded_tu_index_dicts)):
					if j == 1 and gen >= MIN_LATE_CELL_INDEX:
						continue
					if j == 2 and gen < MIN_LATE_CELL_INDEX:
						continue
					if variant not in overcrowded_tu_index_dicts[j][
						tu_index]:
						overcrowded_tu_index_dicts[j][tu_index].append(
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
