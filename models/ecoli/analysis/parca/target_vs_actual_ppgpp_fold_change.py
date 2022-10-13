"""
Analysis that compares the target fold changes (experimentally measured FCs)
against the actual fold changes that are computed for ppGpp-regulated genes.
Discrepancies between operon structures and FC data can result in differences
between these two values.
"""

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

		transcription = sim_data.process.transcription

		# Get fraction of ppgpp-bound RNAPs in rich media
		ppgpp_aa = sim_data.growth_rate_parameters.get_ppGpp_conc(
			sim_data.condition_to_doubling_time['with_aa'])
		f_ppgpp_aa = transcription.fraction_rnap_bound_ppgpp(ppgpp_aa)

		# Get target fold changes for each gene given in experimental data
		cistron_id_to_idx = {
			cistron: i for i, cistron
			in enumerate(transcription.cistron_data['id'])}
		ppgpp_regulated_genes = transcription.ppgpp_regulated_genes
		ppgpp_regulated_gene_indexes = np.array([
			cistron_id_to_idx[cistron_id] for cistron_id in ppgpp_regulated_genes
			])
		target_fold_changes = transcription.ppgpp_fold_changes

		# Get ppgpp-bound and unbound expression levels for each RNA
		exp_ppgpp_rna = transcription.exp_ppgpp
		exp_free_rna = transcription.exp_free

		# Get ppgpp-bound and unbound net expression levels for each gene
		cistron_tu_mapping_matrix = transcription.cistron_tu_mapping_matrix
		exp_ppgpp_gene = cistron_tu_mapping_matrix.dot(exp_ppgpp_rna)
		exp_aa_gene = cistron_tu_mapping_matrix.dot(
			(1 - f_ppgpp_aa) * exp_free_rna + f_ppgpp_aa * exp_ppgpp_rna)

		# Get actual fold changes after applying operon structures
		with np.errstate(divide='ignore', invalid='ignore'):
			all_actual_fold_changes = np.log2(exp_ppgpp_gene / exp_aa_gene)
		actual_fold_changes = all_actual_fold_changes[ppgpp_regulated_gene_indexes]

		plt.figure(figsize=(8, 8))
		plt.scatter(target_fold_changes, actual_fold_changes, s=5)
		plt.xlabel('Target FCs')
		plt.ylabel('Actual FCs')
		plt.xlim([-4, 4])
		plt.ylim([-4, 4])
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
