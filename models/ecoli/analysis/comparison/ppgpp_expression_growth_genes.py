"""
Compares the ppGpp-regulated expression levels of growth genes (genes encoding
for RNAP and ribosome subunits).
"""

from typing import Tuple

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.analysis.analysis_tools import exportFigure
# noinspection PyUnresolvedReferences
from wholecell.io.tablereader import TableReader


CONDITION = 'basal'


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plot_out_dir, plot_out_filename, input_sim_dir, unused, metadata):
		# From fw_queue, reference_sim_dir has operons="off"; input_sim_dir has
		# operons="on".
		# manual/analysisComparison.py can compare any two sim dirs.
		# sim_data1.operons_on and sim_data2.operons_on indicate operons on/off.

		# noinspection PyUnusedLocal
		_, sim_data1, _ = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		_, sim_data2, _ = self.setup(input_sim_dir)

		# Get mask for genes that encode for ribosomal proteins or RNAPs
		is_rnap = sim_data2.process.transcription.cistron_data['is_RNAP']
		is_ribosomal_protein = sim_data2.process.transcription.cistron_data['is_ribosomal_protein']
		growth_genes_mask = np.logical_or(is_rnap, is_ribosomal_protein)
		plotted_cistron_ids = sim_data2.process.transcription.cistron_data['id'][growth_genes_mask]

		def get_exp(sim_data):
			transcription = sim_data.process.transcription

			# Get fraction of ppgpp-bound RNAPs in this condition
			ppgpp_conc = sim_data.growth_rate_parameters.get_ppGpp_conc(
				sim_data.condition_to_doubling_time[CONDITION])

			# Get basal synthesis probabilities for each RNA in this condition
			basal_synth_prob, factor = transcription.synth_prob_from_ppgpp(
				ppgpp_conc, sim_data.process.replication.get_average_copy_number)
			ppgpp_scale = basal_synth_prob.copy()
			ppgpp_scale[ppgpp_scale == 0] = 1

			# Add changes in probabilites from TFs
			p_promoter_bound = np.array([
				sim_data.pPromoterBound[CONDITION][tf]
				for tf in sim_data.process.transcription_regulation.tf_ids
				])
			delta_prob_matrix = sim_data.process.transcription_regulation.get_delta_prob_matrix(
				ppgpp=True)
			synth_prob = basal_synth_prob + ppgpp_scale*(
					delta_prob_matrix @ p_promoter_bound)
			synth_prob[synth_prob < 0] = 0
			synth_prob /= synth_prob.sum()

			# Covert to gene-level expression
			exp_rna = synth_prob / factor
			exp_gene = transcription.cistron_tu_mapping_matrix.dot(exp_rna)

			return exp_gene / exp_gene.sum()

		exp_gene1 = get_exp(sim_data1)[growth_genes_mask]
		exp_gene2 = get_exp(sim_data2)[growth_genes_mask]

		lower_exp_indexes = np.where(exp_gene2 < exp_gene1)[0]

		plt.figure(figsize=(8, 8))
		plt.scatter(exp_gene1, exp_gene2, s=10)
		plt.plot([1e-5, 1e-3], [1e-5, 1e-3], ls='--', lw=3, c='#dddddd')
		plt.xlabel('Gene expression levels, old sims')
		plt.ylabel('Gene expression levels, new sims')
		plt.xlim([1e-5, 1e-3])
		plt.ylim([1e-5, 1e-3])
		plt.xscale('log')
		plt.yscale('log')

		for i in lower_exp_indexes:
			plt.text(1.05 * exp_gene1[i], 1.05 * exp_gene2[i], plotted_cistron_ids[i])

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename + f'_{CONDITION}', metadata)
		plt.close('all')


	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
