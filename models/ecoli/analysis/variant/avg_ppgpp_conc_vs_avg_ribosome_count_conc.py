"""
Plot a scatterplot of the average ppgpps concentration vs. average ribosome count.
This plot is intended to be run on simulations where
the new gene option was enabled.
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.io.tablereader import TableReader


# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variants = self.ap.get_variants()
		min_variant = min(variants)

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in new_gene_mRNA_ids]
		if len(new_gene_mRNA_ids) == 0:
			print("This plot is intended to be run on simulations where the"
				  " new gene option was enabled, but no new gene mRNAs were "
				  "found.")
			return
		if len(new_gene_monomer_ids) == 0:
			print("This plot is intended to be run on simulations where the "
				  "new gene option was enabled, but no new gene proteins "
				  "were "
				  "found.")
			return
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
			'number of new gene monomers and mRNAs should be equal'

		# Data extraction
		plot_variant_mask = np.full(len(variants), True)
		trl_eff_values = np.zeros(len(variants))
		expression_factors = np.zeros(len(variants))
		avg_ribosome_count = np.zeros(len(variants))
		avg_ribosome_conc = np.zeros(len(variants))
		ppgpp_concentration = np.zeros(len(variants))
		n_total_gens = self.ap.n_generation
		variant_name = metadata["variant"]

		# Loop through all variant indexes
		for i, variant in enumerate(variants):

			all_cells = self.ap.get_cells(
				variant=[variant],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				plot_variant_mask[i] = False
				continue

			# Get new gene parameters for this variant index
			condition_index = variant // 1000
			index_remainder = variant - condition_index * 1000

			if variant_name == "new_gene_param_sampling_internal_shift":
				from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift import get_sampled_new_gene_expression_factor_and_translation_efficiency
				np.random.seed(index_remainder)
				expression_factors[i], trl_eff_values[i] = get_sampled_new_gene_expression_factor_and_translation_efficiency(
					index_remainder)

			elif variant_name == "new_gene_param_sampling_internal_shift_narrow":
				from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift_narrow import get_sampled_new_gene_expression_factor_and_translation_efficiency
				expression_factors[i], trl_eff_values[i] = get_sampled_new_gene_expression_factor_and_translation_efficiency(
					index_remainder)

			else:
				print(variant_name + " is not a valid variant name for this plot")
				return

			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')

			# get average ribosome count
			active_ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				remove_first=True,
				ignore_exception=True)[:,ribosome_index]
			avg_ribosome_count[i] = np.mean(active_ribosome_counts)

			# calculate the ribosome concentration
			counts_to_molar = read_stacked_columns(
				all_cells, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True, ignore_exception=True).squeeze()
			avg_ribosome_conc[i] = np.mean(active_ribosome_counts * counts_to_molar)

			# get average concentration of ppgpp
			avg_ppgpp_concentration = read_stacked_columns(
				all_cells, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True, fun=lambda x: np.mean(x)).squeeze()
			ppgpp_concentration[i] = np.mean(avg_ppgpp_concentration)


		# ppgpp concentration vs average ribosome counts plot
		plt.subplot(2, 1, 1)
		plt.xlabel("ppgpp Concentration")
		plt.ylabel("Average ribosome counts")
		plt.scatter(
			ppgpp_concentration[plot_variant_mask],
			avg_ribosome_count[plot_variant_mask])
		plt.plot(
			ppgpp_concentration[plot_variant_mask],
			np.poly1d(np.polyfit(ppgpp_concentration[plot_variant_mask],
								 avg_ribosome_count[plot_variant_mask], 1))(
				ppgpp_concentration[plot_variant_mask]))

		# ppgpp concentration vs average ribosome concentration plot
		plt.subplot(2, 1, 2)
		plt.xlabel("ppgpp Concentration")
		plt.ylabel("Average ribosome concentration")
		plt.scatter(
			ppgpp_concentration[plot_variant_mask],
			avg_ribosome_conc[plot_variant_mask])
		plt.plot(
			ppgpp_concentration[plot_variant_mask],
			np.poly1d(np.polyfit(ppgpp_concentration[plot_variant_mask],
								 avg_ribosome_conc[plot_variant_mask], 1))(
				ppgpp_concentration[plot_variant_mask]))

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()