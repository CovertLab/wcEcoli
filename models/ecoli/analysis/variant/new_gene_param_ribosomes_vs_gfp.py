"""

"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

import pickle
import os
from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure,\
	read_stacked_columns, stacked_cell_threshold_mask
from wholecell.analysis.plotting_tools import DEFAULT_MATPLOTLIB_COLORS\
	as COLORS, labeled_indexable_hist, labeled_indexable_scatter

FONT_SIZE = 9
DOUBLING_TIME_BOUNDS_MINUTES = [0, 180]
N_BINS = 36

MAX_VARIANT = 10 # do not include any variant >= this index
MAX_CELL_INDEX = 8 # do not include any generation >= this index

# Set True to plot completion rate
PLOT_COMPLETION_RATES = True

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

exclude_timeout_cells = 0

MAX_CELL_LENGTH = 180
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# Determine new gene ids


		variants = self.ap.get_variants()

		NEW_GENE_EXPRESSION_FACTOR_CONTROL = 0
		NEW_GENE_EXPRESSION_FACTOR_MIN = 7
		NEW_GENE_EXPRESSION_FACTOR_MAX = 10
		NEW_GENE_TRANSLATION_EFFICIENCY_CONTROL = 0
		NEW_GENE_TRANSLATION_EFFICIENCY_MIN = np.log10(0.01)
		NEW_GENE_TRANSLATION_EFFICIENCY_MAX = np.log10(10)


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
		ribosome_count = []
		generations = {}
		new_gene_mRNA_counts = [{} for id_ in new_gene_mRNA_ids]
		new_gene_monomer_counts = [{} for id_ in new_gene_monomer_ids]
		n_total_gens = self.ap.n_generation


		plt.figure()

		plt.xlabel("Number of Ribosomes")
		plt.ylabel("Protein Counts")
		plt.ylim(0, 1.2e6)

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


			all_cells = self.ap.get_cells(variant=[variant], generation= np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
											  only_successful=True)


			if len(all_cells) == 0:
					continue

			exclude_timeout_cell_mask = stacked_cell_threshold_mask(
				all_cells, 'Main', 'time', MAX_CELL_LENGTH,
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			all_cells_gens = np.array([int(os.path.basename(os.path.dirname(
					cell_path))[-6:]) for cell_path in all_cells])
			generations[variant] = all_cells_gens

			# Doubling times
			dt = read_stacked_columns(all_cells, 'Main', 'time',
										  fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times.append(np.mean(dt[exclude_timeout_cell_mask]))
			if variant == min_variant:
				sim_dir = all_cells[0]
				simOutDir = os.path.join(sim_dir, 'simOut')

				uniqueMoleculeCounts = TableReader(
					os.path.join(simOutDir, "UniqueMoleculeCounts"))
				ribosome_index = uniqueMoleculeCounts.readAttribute(
					"uniqueMoleculeIds").index('active_ribosome')


				# Extract mRNA indexes for each new gene
				mRNA_counts_reader = TableReader(os.path.join(simOutDir,
															  'RNACounts'))
				mRNA_idx_dict = {rna[:-3]: i for i, rna in
					enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
				new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id)
										 for mRNA_id in new_gene_mRNA_ids]

				# Extract protein indexes for each new gene
				monomer_counts_reader = TableReader(os.path.join(
					simOutDir, "MonomerCounts"))
				monomer_idx_dict = {monomer: i for i, monomer in
					enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
				new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id)
					for monomer_id in new_gene_monomer_ids]

			avg_new_gene_mRNA_counts = read_stacked_columns(all_cells,
				'RNACounts', 'mRNA_counts',
				fun=lambda x: np.mean(x[:,new_gene_mRNA_indexes],axis=0))
			avg_new_gene_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				fun=lambda x: np.mean(x[:,new_gene_monomer_indexes],axis=0))
			ribosome_counts = read_stacked_columns(
				all_cells, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', fun=lambda x: np.mean(x[:,ribosome_index],axis=0))

			avg_ng_monomer.append(np.mean(avg_new_gene_monomer_counts[exclude_timeout_cell_mask]))
			ribosome_count.append(np.mean(ribosome_counts[exclude_timeout_cell_mask]))

			for i in range(len(new_gene_mRNA_ids)):
				new_gene_mRNA_counts[i][variant] = \
					np.log10(avg_new_gene_mRNA_counts[:,i] + 1)
				new_gene_monomer_counts[i][variant] = \
					np.log10(avg_new_gene_monomer_counts[:,i] + 1)

		plt.scatter(ribosome_count, avg_ng_monomer)
		#plt.plot(ribosome_count, np.poly1d(np.polyfit(ribosome_count, avg_ng_monomer, 1))(ribosome_count))

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()