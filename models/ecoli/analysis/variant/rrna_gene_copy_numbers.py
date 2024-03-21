"""
Compares the total copy numbers of all genes encoding rRNAs between variants.
"""

import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader


GENE_ID_TO_RRNA_OPERON_ID = {
	'EG30084': 'rrnA',
	'EG30085': 'rrnB',
	'EG30086': 'rrnC',
	'EG30087': 'rrnD',
	'EG30088': 'rrnE',
	'EG30089': 'rrnG',
    'EG30090': 'rrnH',
	}

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if self.ap.n_generation <= 4:
			print('Not enough generations to run analysis.')
			return

		variants = self.ap.get_variants()
		n_variants = len(variants)
		variant_index_to_copy_number = {}

		for variant in variants:
			cell_paths = self.ap.get_cells(variant=[variant])

			# Get gene_ids attribute from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			gene_ids = rna_synth_prob_reader.readAttribute('gene_ids')

			# Get indexes of 16S genes (first gene in each operon)
			rrna_gene_indexes = np.array([
				gene_ids.index(key) for key in GENE_ID_TO_RRNA_OPERON_ID.keys()
				])

			# Get copy numbers of 16S genes
			rrna_gene_copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'gene_copy_number',
				ignore_exception=True, fun=lambda x: x[:, rrna_gene_indexes])
			avg_total_rrna_gene_copy_number = rrna_gene_copy_numbers.sum(axis=1).mean()

			variant_index_to_copy_number[variant] = avg_total_rrna_gene_copy_number

		fig = plt.figure(figsize=(0.8*n_variants, 3))

		# Plot one bar for each variant
		ax0 = fig.add_subplot(1, 1, 1)
		for i, avg_copy_numbers in enumerate(variant_index_to_copy_number.values()):
			ax0.bar(
				i, avg_copy_numbers, width=0.7, alpha=0.5, color='C0')

		ax0.set_xticks(np.arange(n_variants))
		ax0.set_xticklabels(variants)
		ax0.set_xlim([-0.8, n_variants])
		ax0.set_xlabel('Variant indexes')
		ax0.set_ylabel('Total rRNA gene copy numbers')
		ax0.spines['top'].set_visible(False)
		ax0.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
