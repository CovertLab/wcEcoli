"""
Save the average new gene protein counts and doubling time,
This plot is intended to be run on simulations where the new gene
option was enabled.
"""

import numpy as np
from matplotlib import pyplot as plt

import pickle
import os

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader

# Remove first N gens from plot
IGNORE_FIRST_N_GENS = 16

INTERACTIVE_PLOT = True

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
				  "were found.")
			return
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
			'number of new gene monomers and mRNAs should be equal'

		# Data extraction
		plot_variant_mask = np.full(len(variants), True)
		doubling_times = np.zeros(len(variants))
		avg_ng_monomer = np.zeros(len(variants))
		trl_eff_values = np.zeros(len(variants))
		expression_factors = np.zeros(len(variants))
		n_total_gens = self.ap.n_generation

		min_variant = min(variants)

		for i, variant in enumerate(variants):

			all_cells = self.ap.get_cells(
				variant=[variant],
				generation= np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				plot_variant_mask[i] = False
				continue

			# Doubling times
			dt = read_stacked_columns(
				all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[i] = np.mean(dt)
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

			avg_ng_monomer[i] = np.mean(avg_new_gene_monomer_counts)

		# Save variant index, doubling time, and new gene protein counts
		# to file
		all_avg_monomer_counts = np.vstack((variants, doubling_times,
			avg_ng_monomer)).T
		header = np.array(['Variant Index', 'Doubling Time (min)',
			'New Gene Protein Counts'])
		all_avg_monomer_counts = np.vstack((header, all_avg_monomer_counts))
		np.savetxt(os.path.join(plotOutDir, plotOutFileName + ".csv"),
			all_avg_monomer_counts, delimiter=",", fmt="%s")

		plt.figure()
		plt.xlabel("New Gene Protein Counts")
		plt.ylabel("Doubling Time")

		plt.scatter(
			avg_ng_monomer[plot_variant_mask],
			doubling_times[plot_variant_mask])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		if INTERACTIVE_PLOT:
			import plotly.express as px

			fig = px.scatter(
				x=avg_ng_monomer[plot_variant_mask],
				y=doubling_times[plot_variant_mask],
				hover_name=np.array(variants)[plot_variant_mask],
				labels={'x': 'New Gene Protein Counts', 'y': 'Doubling Time'},
				hover_data={
					'Variant Index': np.array(variants)[plot_variant_mask]})

			fig.write_html(os.path.join(
				plotOutDir, plotOutFileName + ".html"))
			fig.show()

if __name__ == "__main__":
	Plot().cli()