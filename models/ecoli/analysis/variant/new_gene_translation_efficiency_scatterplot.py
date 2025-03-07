"""
Plot translation efficiency of new gene vs
- new gene count
- doubling time
- Average number of ribosomes
- Average number of RNA polymerases
- Average ppGpp concentration
- average cell mass
"""

from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
import numpy as np

from wholecell.analysis.analysis_tools import exportFigure, \
	read_stacked_columns, read_stacked_bulk_molecules, \
	stacked_cell_identification
from wholecell.analysis.plotting_tools import heatmap
from unum.units import fg

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift import (
	NEW_GENE_EXPRESSION_FACTOR_CONTROL, NEW_GENE_EXPRESSION_FACTOR_MIN,
	NEW_GENE_EXPRESSION_FACTOR_MAX, NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL,
	NEW_GENE_TRANSLATION_EFFICIENCY_MIN, NEW_GENE_TRANSLATION_EFFICIENCY_MAX)
from wholecell.io.tablereader import TableReader

import os.path
import pickle


IGNORE_FIRST_N_GENS = 0

"""
PLOTS_LIST = ["translation_efficiency_vs_doubling_times",
				"translation_efficiency_vs_cell_mass",
				"translation_efficiency_vs_ppgpp_concentration",
				"translation_efficiency_vs_rnap_counts",
				"translation_efficiency_vs_ribosome_counts",
				"translation_efficiency_vs_new_gene_monomer_counts"
				]
"""
class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		variants = self.ap.get_variants()

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

		avg_doubling_time = np.zeros(len(variants))
		avg_cell_mass = np.zeros(len(variants))
		ppgpp_concentration = np.zeros(len(variants))
		avg_active_rnap = np.zeros(len(variants))
		avg_ribosomes_count = np.zeros(len(variants))
		avg_ng_monomer = np.zeros(len(variants))


		n_total_gens = self.ap.n_generation
		variant_name = metadata["variant"]
		min_variant = min(variants)

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
				from models.ecoli.sim.variants.new_gene_param_sampling_internal_shift import \
					get_sampled_new_gene_expression_factor_and_translation_efficiency
				np.random.seed(index_remainder)
				expression_factors[i], trl_eff_values[
					i] = get_sampled_new_gene_expression_factor_and_translation_efficiency(
					index_remainder)

			elif variant_name == "new_gene_param_sampling_internal_shift_narrow":
				params_to_use = metadata["params_to_use"]
				if variant == 0:
					expression_factors[i] = 0
					trl_eff_values[i] = 0

				else:
					expression_factors[i] = float(params_to_use[str(variant)]["expression_factor"])
					trl_eff_values[i] = float(params_to_use[str(variant)]["translation_efficiency"])

			else:
				print(variant_name + " is not a valid variant name for this plot")
				return

			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				monomer_counts_reader = TableReader(os.path.join(
					simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
									enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id)
											for monomer_id in new_gene_monomer_ids]
				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				active_rnap_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_RNAP')
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')


			# get average doubling time
			dt = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			avg_doubling_time[i] = np.mean(dt)

			# get average cell mass
			mean_cell_mass = read_stacked_columns(
				all_cells, 'Mass', 'cellMass',
				remove_first=True,
				fun=np.mean).reshape(-1)
			avg_cell_mass[i] = np.mean(mean_cell_mass)


			# get ppgpp concentration
			avg_ppgpp_concentration = read_stacked_columns(
				all_cells, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True, fun=lambda x: np.mean(x)).squeeze()
			ppgpp_concentration[i] = np.mean(avg_ppgpp_concentration)

			# get average RNA polymerase count
			active_rnap_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				fun=lambda x: np.mean(x[:, active_rnap_index], axis=0))
			avg_active_rnap[i] = np.mean(active_rnap_counts)

			# get average ribosomes count
			ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				fun=lambda x: np.mean(x[:, ribosome_index], axis=0))
			avg_ribosomes_count[i] = np.mean(ribosome_counts)

			# get average new gene monomer counts
			avg_new_gene_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				fun=lambda x: np.mean(x[:, new_gene_monomer_indexes],axis=0))
			avg_ng_monomer[i] = np.mean(avg_new_gene_monomer_counts)

		# plotting
			grid_spec = GridSpec(6, 1)
			plt.figure(figsize=(5, 15))

			# # translation_efficiency_vs_doubling_times
			ax = plt.subplot(grid_spec[0, 0])
			ax.scatter(
				avg_doubling_time[plot_variant_mask],
				trl_eff_values[plot_variant_mask])
			ax.plot(
				avg_doubling_time[plot_variant_mask],
				np.poly1d(np.polyfit(avg_doubling_time[plot_variant_mask],
									 trl_eff_values[plot_variant_mask], 1))(
					avg_doubling_time[plot_variant_mask]))
			ax.set_xlabel("average doubling time")
			ax.set_ylabel("translation efficiency")
			self.remove_border(ax=ax, bottom=True)

			# translation_efficiency_vs_cell_mass
			ax = plt.subplot(grid_spec[1, 0])

			ax.scatter(
				avg_cell_mass[plot_variant_mask],
				trl_eff_values[plot_variant_mask])
			ax.plot(
				avg_cell_mass[plot_variant_mask],
				np.poly1d(np.polyfit(avg_cell_mass[plot_variant_mask],
									 trl_eff_values[plot_variant_mask], 1))(
					avg_cell_mass[plot_variant_mask]))
			ax.set_xlabel("average cell mass")
			ax.set_xlabel("translation efficiency")
			self.remove_border(ax=ax, bottom=True)

			# translation_efficiency_vs_ppgpp_concentration 3
			ax = plt.subplot(grid_spec[2, 0])
			ax.scatter(
				ppgpp_concentration[plot_variant_mask],
				trl_eff_values[plot_variant_mask])
			ax.plot(
				ppgpp_concentration[plot_variant_mask],
				np.poly1d(np.polyfit(ppgpp_concentration[plot_variant_mask],
									 trl_eff_values[plot_variant_mask], 1))(
					ppgpp_concentration[plot_variant_mask]))

			ax.set_xlabel("average ppgpp concentration")
			ax.set_ylabel("translation efficiency")
			self.remove_border(ax=ax, bottom=True)

			# translation_efficiency_vs_rnap_counts 4
			ax = plt.subplot(grid_spec[3, 0])
			ax.scatter(
				avg_active_rnap[plot_variant_mask],
				trl_eff_values[plot_variant_mask])
			ax.plot(
				avg_active_rnap[plot_variant_mask],
				np.poly1d(np.polyfit(avg_active_rnap[plot_variant_mask],
									 trl_eff_values[plot_variant_mask], 1))(
					avg_active_rnap[plot_variant_mask]))

			ax.set_xlabel("average active RNA polymerase")
			ax.set_ylabel("translation efficiency")
			self.remove_border(ax=ax, bottom=True)

			# translation_efficiency_vs_ribosome_counts 5
			ax = plt.subplot(grid_spec[4, 0])
			ax.scatter(
				avg_ribosomes_count[plot_variant_mask],
				trl_eff_values[plot_variant_mask])
			ax.plot(
				avg_ribosomes_count[plot_variant_mask],
				np.poly1d(np.polyfit(avg_ribosomes_count[plot_variant_mask],
									 trl_eff_values[plot_variant_mask], 1))(
					avg_ribosomes_count[plot_variant_mask]))

			ax.set_xlabel("average active ribosome counts")
			ax.set_ylabel("translation efficiency")
			self.remove_border(ax=ax, bottom=True)

			# translation_efficiency_vs_new_gene_monomer_counts 6
			ax = plt.subplot(grid_spec[5, 0])
			ax.scatter(
				avg_ng_monomer[plot_variant_mask],
				trl_eff_values[plot_variant_mask])
			ax.plot(
				avg_ng_monomer[plot_variant_mask],
				np.poly1d(np.polyfit(avg_ng_monomer[plot_variant_mask],
									 trl_eff_values[plot_variant_mask], 1))(
					avg_ng_monomer[plot_variant_mask]))

			ax.set_xlabel("average new gene monomer counts")
			ax.set_ylabel("translation efficiency")
			self.remove_border(ax=ax, bottom=True)

			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + "IGN_FIRST" + str(IGNORE_FIRST_N_GENS) + "GEN", metadata)
			plt.close('all')

if __name__ == '__main__':
	Plot().cli()