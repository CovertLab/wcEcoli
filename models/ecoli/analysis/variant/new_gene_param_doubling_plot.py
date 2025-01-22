"""
Plot a scatterplot of the average new gene protein counts vs. doubling time,
with the option to color the points by expression factor or translation
efficiency. This plot is intended to be run on simulations where the new gene
option was enabled.
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift import (
	NEW_GENE_EXPRESSION_FACTOR_CONTROL, NEW_GENE_EXPRESSION_FACTOR_MIN, NEW_GENE_EXPRESSION_FACTOR_MAX,
	NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL, NEW_GENE_TRANSLATION_EFFICIENCY_MIN,
	NEW_GENE_TRANSLATION_EFFICIENCY_MAX)

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

COLOR_BY = "same"  # ["same", "expression_factor", "translation_efficiency"]

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# Determine new gene ids
		variants = self.ap.get_variants()

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
		doubling_times = []
		avg_ng_monomer = []
		trl_eff_values = []
		expression_factors = []
		colors = []
		n_total_gens = self.ap.n_generation

		min_variant = min(variants)

		for variant in variants:
			np.random.seed(variant)
			if variant == 0:
				expression_factors.append(NEW_GENE_EXPRESSION_FACTOR_CONTROL)
				trl_eff_values.append(NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL)
			else:
				expression_factors.append(np.random.uniform(NEW_GENE_EXPRESSION_FACTOR_MIN,
													  NEW_GENE_EXPRESSION_FACTOR_MAX))
				trl_eff_values.append(10 ** np.random.uniform(NEW_GENE_TRANSLATION_EFFICIENCY_MIN,
														NEW_GENE_TRANSLATION_EFFICIENCY_MAX))

			all_cells = self.ap.get_cells(
				variant=[variant],
				generation= np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
					continue

			# Doubling times
			dt = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times.append(np.mean(dt))
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(
					simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
					enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id)
					for monomer_id in new_gene_monomer_ids]

			avg_new_gene_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				fun=lambda x: np.mean(x[:,new_gene_monomer_indexes],axis=0))

			avg_ng_monomer.append(np.mean(avg_new_gene_monomer_counts))

			if COLOR_BY == "translation_efficiency":
				colors.append(trl_eff_values[variant]/10)
			elif COLOR_BY == "expression_factor":
				colors.append(expression_factors[variant] / 10)

		plt.figure()
		plt.xlabel("Protein Counts")
		plt.ylabel("Doubling Time")

		if COLOR_BY != "same":
			plt.scatter(avg_ng_monomer, doubling_times, c=colors, cmap='coolwarm')
			plt.colorbar(orientation='horizontal', label='expression factor / 10')
		else:
			plt.scatter(avg_ng_monomer, doubling_times)

		plt.plot(avg_ng_monomer,
				 np.poly1d(np.polyfit(avg_ng_monomer, doubling_times, 1))(avg_ng_monomer))

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_" + COLOR_BY, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
