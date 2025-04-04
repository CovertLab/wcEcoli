"""
Computes and plots average init probs for proteins involved in ribosomes or
RNA polymerases
"""

import os
import pickle

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd

from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader

from models.ecoli.analysis import cohortAnalysisPlot

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		if "init_sims" in metadata:
			num_seeds = metadata['init_sims']
		else:
			num_seeds = metadata['total_init_sims']

		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False

		RNAP_component_monomer_ids = sim_data.molecule_groups.RNAP_subunits
		ribosome_component_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
		monomer_ids_to_plot = RNAP_component_monomer_ids + ribosome_component_monomer_ids

		total_plots = len(monomer_ids_to_plot)
		avg_protein_init_probs = np.zeros((total_plots, num_seeds))

		for seed in range(num_seeds):
			print("Seed: ", seed)

			plt.figure(figsize=(8.5, total_plots * 3))
			plot_num = 1

			cell_paths = self.ap.get_cells(seed=[seed])

			if len(cell_paths) == 0:
				continue

			if seed == 0:
				# Get monomer indexes for each monomer_id_to_plot
				monomer_reader = TableReader(
					os.path.join(cell_paths[0], 'simOut', 'MonomerCounts'))
				monomer_ids = np.array(monomer_reader.readAttribute('monomerIds'))
				monomer_id_to_idx_dict = {
					monomer_id: i for i, monomer_id in enumerate(monomer_ids)}

			for i, monomer_id in enumerate(monomer_ids_to_plot):

				monomer_idx = monomer_id_to_idx_dict[monomer_id]

				if i == 0:
					ax1 = plt.subplot(total_plots, 1, plot_num)
				else:
					plt.subplot(total_plots, 1, plot_num, sharex=ax1)

				time = read_stacked_columns(
					cell_paths, 'Main', 'time',
					ignore_exception=True)
				protein_init_prob = read_stacked_columns(
					cell_paths, 'RibosomeData',
					'target_prob_translation_per_transcript',
					ignore_exception=True)[:, monomer_idx]

				plt.plot(time / 60.0, protein_init_prob)

				# Add a line for the average
				avg_protein_init_prob = np.mean(protein_init_prob, axis=0)
				avg_protein_init_probs[i, seed] = avg_protein_init_prob
				plt.axhline(y=avg_protein_init_prob, color='r', linestyle='--')

				plt.xlabel("Time (min)", fontsize="small")
				plt.ylabel(
					f"Protein Init Prob: {monomer_id}",fontsize="small")
				plt.title(
					f"Protein Init Prob: {monomer_id}", fontsize="small")

				plot_num += 1

			plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
			exportFigure(
				plt, plotOutDir, plotOutFileName + f"_seed_{seed}", metadata)
			plt.close("all")

		# for each monomer_id, print the averages for each seed
		for i, monomer_id in enumerate(monomer_ids_to_plot):
			print(f"\n\nAvg init probs for {monomer_id}:")
			print(avg_protein_init_probs[i, :])

			print("Average of the averages: ")
			print(np.mean(avg_protein_init_probs[i, :]))

		# Save a csv of protein id, average init probs for each monomer_id
		avg_protein_init_probs_df = pd.DataFrame(
			avg_protein_init_probs.T, columns=monomer_ids_to_plot)
		avg_protein_init_probs_df.to_csv(
			os.path.join(plotOutDir, "avg_protein_init_probs.csv"),
			index=False)

		avg_of_avg_init_probs = np.mean(avg_protein_init_probs, axis=1)
		avg_of_avg_init_probs = np.reshape(avg_of_avg_init_probs, (len(avg_of_avg_init_probs), 1))
		avg_of_avg_init_probs_df = pd.DataFrame(
			avg_of_avg_init_probs.T, columns=monomer_ids_to_plot)
		avg_of_avg_init_probs_df.to_csv(
			os.path.join(plotOutDir, "avg_of_avg_init_probs.csv"),
			index=False)


if __name__ == "__main__":
	Plot().cli()
