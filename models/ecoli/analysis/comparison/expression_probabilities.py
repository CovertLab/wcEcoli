"""
Compare expression probabilities.
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


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):

		# noinspection PyUnusedLocal
		ap1, sim_data1, validation_data1 = self.setup(reference_sim_dir)
		# noinspection PyUnusedLocal
		ap2, sim_data2, validation_data2 = self.setup(input_sim_dir)

		# sim_data classes used
		t_reg1 = sim_data1.process.transcription_regulation
		transcription1 = sim_data1.process.transcription
		replication1 = sim_data1.process.replication

		t_reg2 = sim_data2.process.transcription_regulation
		transcription2 = sim_data2.process.transcription
		replication2 = sim_data2.process.replication

		# Normal expression
		tf_to_gene_id1 = t_reg1.tf_to_gene_id
		tf_ids1 = [tf_to_gene_id1[tf] for tf in t_reg1.tf_ids]
		basal_prob1 = t_reg1.basal_prob
		rna_ids1 = transcription1.rna_data['id']
		gene_symbols1 = [
			sim_data1.common_names.get_common_name(rna_id[:-3])
			for rna_id in rna_ids1]

		tf_to_gene_id2 = t_reg2.tf_to_gene_id
		tf_ids2 = [tf_to_gene_id2[tf] for tf in t_reg2.tf_ids]
		basal_prob2 = t_reg2.basal_prob
		rna_ids2 = transcription2.rna_data['id']
		gene_symbols2 = [
			sim_data2.common_names.get_common_name(rna_id[:-3])
			for rna_id in rna_ids2]

		# TODO: Make more elegant/generalizable
		# basal_prob2 = basal_prob2[:-1]
		# gene_symbols2 = gene_symbols2[:-1]
		basal_prob1 = np.append(basal_prob1, 0.0)
		gene_symbols1.append(gene_symbols2[-1])
		assert gene_symbols1 == gene_symbols2

		abs_diff = basal_prob2 - basal_prob1
		rel_diff = (basal_prob2 - basal_prob1) / basal_prob1

		mag_rel_diff = np.absolute(rel_diff)
		descending_mag_rel_diff_indices = np.argsort(-mag_rel_diff)

		quantiles = [np.nanquantile(abs_diff, 0),
			np.nanquantile(abs_diff, .25),
			np.nanquantile(abs_diff, .5),
			np.nanquantile(abs_diff, .75),
			np.nanquantile(abs_diff, 1)]

		rel_quantiles = [np.nanquantile(rel_diff, 0),
			np.nanquantile(rel_diff, .25),
			np.nanquantile(rel_diff, .5),
			np.nanquantile(rel_diff, .75),
			np.nanquantile(rel_diff, 1)]

		mag_rel_quantiles = [np.nanquantile(mag_rel_diff, 0),
			np.nanquantile(mag_rel_diff, .25),
			np.nanquantile(mag_rel_diff, .5),
			np.nanquantile(mag_rel_diff, .75),
			np.nanquantile(mag_rel_diff, 1)]

		print("Comparing basal probabilities")
		print("Quantiles (0 25 50 75 100): ", quantiles)
		print("Relative Quantiles (0 25 50 75 100): ", rel_quantiles)
		print("Magnitude Relative Quantiles (0 25 50 75 100): ", mag_rel_quantiles)
		print("\n")
		#
		# np.set_printoptions(threshold=np.inf)
		# print(np.array(gene_symbols1)[descending_mag_rel_diff_indices])
		#
		# for j in range(len(descending_mag_rel_diff_indices)):
		# 	if j < 200:
		# 		i = descending_mag_rel_diff_indices[j]
		# 		print("\n")
		# 		print(i)
		# 		print(gene_symbols1[i])
		# 		print("Rel. Diff: ", rel_diff[i])
		# 		print("Abs. Diff: ", abs_diff[i])
		# 		print("basal prob ref: ", basal_prob1[i])
		# 		print("basal prob input: ", basal_prob2[i])





		RNA_data = sim_data2.process.transcription.rna_data.struct_array
		RNAP_subunit_mask = RNA_data['includes_RNAP']
		rRNA_mask = RNA_data['is_rRNA']
		ribosomal_protein_mask = RNA_data['includes_ribosomal_protein']

		fig, ax = plt.subplots()
		plt.scatter(basal_prob1[np.where(ribosomal_protein_mask)[0]],
					basal_prob2[np.where(ribosomal_protein_mask)[0]],
					color="red", label="Ribosomal Protein")
		plt.scatter(basal_prob1[np.where(RNAP_subunit_mask)[0]],
					basal_prob2[np.where(RNAP_subunit_mask)[0]],
					color = "blue",label = "RNAP Subunit")
		plt.scatter(basal_prob1[np.where(rRNA_mask)[0]],
					basal_prob2[np.where(rRNA_mask)[0]],
					color="green", label="rRNA")
		plt.title("Basal Probability Comparison")
		plt.xlabel("Basal Probabaility Reference")
		plt.ylabel("Basal Probability Input")
		plt.legend()
		ax.axline((0, 0), slope=1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


		print("\n\nRNAP Subunit Basal Prob Comparisons")
		for i in np.where(RNAP_subunit_mask)[0]:
			print("\n")
			print(i)
			print(gene_symbols1[i])
			print("Rel. Diff: ", rel_diff[i])
			print("Abs. Diff: ", abs_diff[i])
			print("basal prob ref: ", basal_prob1[i])
			print("basal prob input: ", basal_prob2[i])

		print("\n\nrRNA Basal Prob Comparisons")
		for i in np.where(rRNA_mask)[0]:
			print("\n")
			print(i)
			print(gene_symbols1[i])
			print("Rel. Diff: ", rel_diff[i])
			print("Abs. Diff: ", abs_diff[i])
			print("basal prob ref: ", basal_prob1[i])
			print("basal prob input: ", basal_prob2[i])

		print("\n\nRibosomal Protein Basal Prob Comparisons")
		for i in np.where(ribosomal_protein_mask)[0]:
			print("\n")
			print(i)
			print(gene_symbols1[i])
			print("Rel. Diff: ", rel_diff[i])
			print("Abs. Diff: ", abs_diff[i])
			print("basal prob ref: ", basal_prob1[i])
			print("basal prob input: ", basal_prob2[i])

	def setup(self, inputDir: str) -> Tuple[
			AnalysisPaths, SimulationDataEcoli, ValidationDataEcoli]:
		"""Return objects used for analyzing multiple sims."""
		ap = AnalysisPaths(inputDir, variant_plot=True)
		sim_data = self.read_sim_data_file(inputDir)
		validation_data = self.read_validation_data_file(inputDir)
		return ap, sim_data, validation_data


if __name__ == "__main__":
	Plot().cli()
