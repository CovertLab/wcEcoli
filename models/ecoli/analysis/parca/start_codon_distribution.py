"""
Plots the distribution of start codons and N-terminal amino acids for all
protein-encoding annotated genes used by the model.
"""

from collections import Counter
import pickle

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		if sim_data.operons_on:
			print('The operons option must be set to off for this plot -- skipping analysis.')
			return

		# Get distributions of start codons and N-terminal amino acids
		monomer_data = sim_data.process.translation.monomer_data
		start_codons = [
			str(seq[:3]) for seq
			in sim_data.getter.get_sequences(monomer_data['cistron_id'])
			]
		n_terminal_amino_acids = [
			seq[0] for seq in sim_data.getter.get_sequences(
				[protein_id[:-3] for protein_id in monomer_data['id']]
				)
			]
		start_codon_counter = Counter(start_codons)
		n_terminal_amino_acid_counter = Counter(n_terminal_amino_acids)

		# Plot distributions
		plt.figure(figsize=(6, 6))
		ax0 = plt.subplot(1, 2, 1)
		ax0.bar(
			range(len(start_codon_counter)),
			start_codon_counter.values())
		ax0.set_xticks(range(len(start_codon_counter)))
		ax0.set_xticklabels(start_codon_counter.keys())
		ax0.set_xlabel("Start codon")
		ax0.set_ylabel("Number of genes")

		ax1 = plt.subplot(1, 2, 2, sharey=ax0)
		ax1.bar(
			range(len(n_terminal_amino_acid_counter)),
			n_terminal_amino_acid_counter.values())
		ax1.set_xticks(range(len(n_terminal_amino_acid_counter)))
		ax1.set_xticklabels(n_terminal_amino_acid_counter.keys())
		ax1.set_xlim([-2.7, 2.7])
		ax1.set_xlabel("N-terminal amino acid")

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
