"""
Plots the differences between the relative expression levels of each tRNA
cistron that is expected from experimental data (Dong et al., 1996) and the
actual expression levels of each tRNA cistron after accounting for operon
structure.
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils.fitting import normalize


CONDITION = 'basal'
NUMERICAL_ZERO = 1e-20
LABEL_BOUNDARY = 2


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load data from sim_data
		transcription = sim_data.process.transcription
		all_cistron_ids = transcription.cistron_data['id']
		all_cistron_names = [
			sim_data.common_names.get_common_name(cistron_id)
			for cistron_id in all_cistron_ids
			]
		cistron_is_tRNA = transcription.cistron_data['is_tRNA']
		rna_includes_tRNA = transcription.rna_data['includes_tRNA']
		tRNA_cistron_ids = all_cistron_ids[cistron_is_tRNA]

		# Get boolean array of relevant cistrons that belong to at least one
		# polycistronic transcript
		polycistronic_cistron_indexes = []
		for rna_id in sim_data.process.transcription.rna_data['id']:
			cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes(rna_id)
			if len(cistron_indexes) > 1:
				polycistronic_cistron_indexes.extend(cistron_indexes)
		is_polycistronic = np.zeros(len(all_cistron_ids), bool)
		if len(polycistronic_cistron_indexes) > 0:
			is_polycistronic[np.array(
				list(set(polycistronic_cistron_indexes)))] = True
		is_polycistronic = is_polycistronic[cistron_is_tRNA]

		# Get expected expression levels of each tRNA cistron
		doubling_time = sim_data.condition_to_doubling_time[CONDITION]
		tRNA_distribution = sim_data.mass.get_trna_distribution(doubling_time)
		tRNA_id_to_dist = {
			trna_id: dist for (trna_id, dist)
			in zip(tRNA_distribution['id'], tRNA_distribution['molar_ratio_to_16SrRNA'])
			}
		expected_tRNA_cistron_exp = np.zeros(len(tRNA_cistron_ids))
		for i, tRNA_id in enumerate(tRNA_cistron_ids):
			expected_tRNA_cistron_exp[i] = tRNA_id_to_dist[tRNA_id]
		expected_tRNA_cistron_exp = normalize(expected_tRNA_cistron_exp)

		# Get actual expression levels of each tRNA cistron
		rna_exp = transcription.rna_expression[CONDITION]
		actual_tRNA_cistron_exp = transcription.tRNA_cistron_tu_mapping_matrix.dot(
			rna_exp[rna_includes_tRNA])
		actual_tRNA_cistron_exp = normalize(actual_tRNA_cistron_exp)

		# Find cistrons with more than a 10-fold difference between actual vs
		# expected expression
		exp_diff = np.log10(actual_tRNA_cistron_exp + NUMERICAL_ZERO) - np.log10(expected_tRNA_cistron_exp + NUMERICAL_ZERO)
		largest_diff_indexes = np.where(np.abs(exp_diff) > np.log10(LABEL_BOUNDARY))[0]
		largest_diff_indexes = [index for index in largest_diff_indexes]

		plt.figure(figsize=(9, 9))
		ls = np.logspace(-3, -1, 2)
		plt.plot(ls, LABEL_BOUNDARY * ls, c='#dddddd', ls='--')
		plt.plot(ls, 1 / LABEL_BOUNDARY * ls, c='#dddddd', ls='--')

		plt.scatter(
			expected_tRNA_cistron_exp, actual_tRNA_cistron_exp,
			c='#dddddd', s=2, label='monocistronic')

		# Highlight genes that are polycistronic
		plt.scatter(
			expected_tRNA_cistron_exp[is_polycistronic],
			actual_tRNA_cistron_exp[is_polycistronic],
			c='#333333', s=4, label='polycistronic')

		# Label tRNAs that have more than a 10-fold difference
		for index in largest_diff_indexes:
			plt.text(
				expected_tRNA_cistron_exp[index],
				actual_tRNA_cistron_exp[index] * 1.03,
				all_cistron_names[index],
				ha='center', va='bottom', fontsize=5)

		plt.xlabel('Expression expected from Dong et al. (1996)')
		plt.ylabel('Actual expression after applying operon structure')
		plt.xlim([1e-3, 5e-2])
		plt.ylim([1e-3, 5e-2])
		plt.xscale('log')
		plt.yscale('log')
		plt.legend()

		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
