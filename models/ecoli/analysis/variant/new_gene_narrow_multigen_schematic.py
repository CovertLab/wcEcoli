"""
Analyze how production machinery counts change over time and across at most two
variants.
"""

import pickle
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (
	exportFigure, read_stacked_columns, read_stacked_bulk_molecules)
from wholecell.io.tablereader import TableReader

poster_colors = {
    "light_gray": (0.75, 0.75, 0.75),
    "poster_green": (66/255, 170/255, 154/255),
    "poster_blue": (27/255, 132/255, 198/255),
    "poster_purple": (188/255, 140/255, 191/255),
    "poster_gold": (221/255, 203/255, 119/255),
    "poster_light_blue": (136/255, 205/255, 240/255),
    "poster_red": (202/255, 0/255, 32/255),
}

START_GEN_INDEX = 7
END_GEN_INDEX = 11 # Not inclusive
GEN_RANGE = np.arange(START_GEN_INDEX, END_GEN_INDEX)

SELECTED_SEED_INDEXES = [
			0, 1, 2, 3, 4, 5, 6, 7]

GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: poster_colors["poster_green"],
	9: poster_colors["poster_blue"],
	10: poster_colors["poster_purple"],
	}
VARIANT_B_COLOR = poster_colors["light_gray"]

# FOR SHERLOCK
VARIANT_INDEX_A = 23
VARIANT_INDEX_B = 0 # Optional, set to -1 if only one variant is desired

# FOR LOCAL DEVELOPMENT
VARIANT_INDEX_A = 21
VARIANT_INDEX_B = 0 # Optional, set to -1 if only one variant is desired

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		plotOutDir = os.path.join(plotOutDir, plotOutFileName)

		for seed_index in SELECTED_SEED_INDEXES:
			all_cells_A = self.ap.get_cells(
				variant=[VARIANT_INDEX_A],
				seed=[seed_index],
				generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
				only_successful=True)
			if len(all_cells_A) == 0:
				continue
			if VARIANT_INDEX_B != -1:
				all_cells_B = self.ap.get_cells(
					variant=[VARIANT_INDEX_B],
					seed=[seed_index],
					generation=np.arange(START_GEN_INDEX, END_GEN_INDEX),
					only_successful=True)
				if len(all_cells_B) == 0:
					continue

			total_plots = 15 # TODO: Modularize and get rid of this magic number
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.figure(figsize=(14, total_plots * 3))

			plot_num = 1
			ax1 = plt.subplot(total_plots, 1, plot_num)

			# Load generation labels
			dt_A = read_stacked_columns(
				all_cells_A, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			num_time_steps_A = read_stacked_columns(
				all_cells_A, 'Main', 'time',
				fun=lambda x: len(x)).squeeze()
			gen_labels_A = np.repeat(np.arange(len(dt_A)), num_time_steps_A)
			unique_gen_labels_A = np.unique(gen_labels_A)
			gen_start_indexes_A = np.array(
				[gen_labels_A.tolist().index(i) for i in unique_gen_labels_A])
			gen_end_indexes_A = np.concatenate((
				np.array(gen_start_indexes_A[1:] - 1), np.array([len(gen_labels_A) - 1])))
			time_A = read_stacked_columns(
				all_cells_A, 'Main', 'time', ignore_exception=True)

			if VARIANT_INDEX_B != -1:
				dt_B = read_stacked_columns(
					all_cells_B, 'Main', 'time',
					fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
				num_time_steps_B = read_stacked_columns(
					all_cells_B, 'Main', 'time',
					fun=lambda x: len(x)).squeeze()
				gen_labels_B = np.repeat(np.arange(len(dt_B)), num_time_steps_B)
				unique_gen_labels_B = np.unique(gen_labels_B)
				gen_start_indexes_B = np.array(
					[gen_labels_B.tolist().index(i) for i in unique_gen_labels_B])
				gen_end_indexes_B = np.concatenate((
					np.array(gen_start_indexes_B[1:] - 1), np.array([len(gen_labels_B) - 1])))
				time_B = read_stacked_columns(
					all_cells_B, 'Main', 'time', ignore_exception=True)

			# Determine new gene ids
			with open(simDataFile, 'rb') as f:
				sim_data = pickle.load(f)
			mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
			monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
			new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
			mRNA_monomer_id_dict = dict(zip(
				monomer_sim_data['cistron_id'], monomer_sim_data['id']))
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

			# Extract mRNA indexes for each new gene
			sim_dir = all_cells_A[0]
			simOutDir = os.path.join(sim_dir, 'simOut')
			mRNA_counts_reader = TableReader(os.path.join(
				simOutDir, 'RNACounts'))
			mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(
				mRNA_counts_reader.readAttribute('mRNA_ids'))}
			new_gene_mRNA_indexes = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
			mRNA_counts_reader.close()

			# Extract protein indexes for each new gene
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in
				enumerate(monomer_counts_reader.readAttribute('monomerIds'))}
			new_gene_monomer_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]
			monomer_counts_reader.close()

			# New Gene mRNA Counts
			new_gene_mRNA_counts_A = read_stacked_columns(
				all_cells_A, 'RNACounts',
		'mRNA_counts')[:, new_gene_mRNA_indexes]
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = new_gene_mRNA_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1, :]
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				new_gene_mRNA_counts_B = read_stacked_columns(
					all_cells_B, 'RNACounts',
				'mRNA_counts')[:, new_gene_mRNA_indexes]
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = new_gene_mRNA_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1, :]
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('New Gene mRNA Counts')
			plt.xlabel('Time (minutes)')
			handles, labels = plt.gca().get_legend_handles_labels()
			by_label = dict(zip(labels, handles))
			plt.legend(
				by_label.values(), by_label.keys(),
				loc='upper left', fontsize=8)
			plot_num += 1

			# New Gene Protein Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			new_gene_protein_counts_A = read_stacked_columns(
				all_cells_A, 'MonomerCounts',
				'monomerCounts')[:, new_gene_monomer_indexes]
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = new_gene_protein_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				new_gene_protein_counts_B = read_stacked_columns(
					all_cells_B, 'MonomerCounts',
					'monomerCounts')[:, new_gene_monomer_indexes]
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = new_gene_protein_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('New Gene Protein Counts')
			plt.xlabel('Time (minutes)')
			handles, labels = plt.gca().get_legend_handles_labels()
			by_label = dict(zip(labels, handles))
			plt.legend(
				by_label.values(), by_label.keys(),
				loc='upper left', fontsize=8)
			plot_num += 1


			# Save figure
			print(f"\nSeed: {seed_index}")
			print("Total number of plots made: ", plot_num - 1)
			plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
			if VARIANT_INDEX_B != -1:
				variants_str = f"{VARIANT_INDEX_A}_vs_{VARIANT_INDEX_B}"
			else:
				variants_str = f"{VARIANT_INDEX_A}"
			exportFigure(
				plt, plotOutDir,
				plotOutFileName + f"_seed_{seed_index}_var_{variants_str}",
				metadata)
			plt.close("all")


if __name__ == "__main__":
	Plot().cli()
