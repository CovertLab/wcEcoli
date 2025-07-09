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

DEFAULT_COLOR = poster_colors["light_gray"]

INDIVIDUAL_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: poster_colors["poster_green"],
	9: poster_colors["poster_blue"],
	10: poster_colors["poster_purple"],
	}
VARIANT_B_COLOR = poster_colors["light_gray"]

GFP_INDUCED_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: poster_colors["poster_green"],
	9: DEFAULT_COLOR,
	10: DEFAULT_COLOR
}

RNAP_RIBO_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: poster_colors["poster_blue"],
	9: DEFAULT_COLOR,
	10: DEFAULT_COLOR
}

CELL_MASS_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: DEFAULT_COLOR,
	9: poster_colors["poster_gold"],
	10: DEFAULT_COLOR
	}

FORKS_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: DEFAULT_COLOR,
	9: DEFAULT_COLOR,
	10: poster_colors["poster_light_blue"]
}

DT_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: DEFAULT_COLOR,
	9: DEFAULT_COLOR,
	10: poster_colors["poster_purple"]
}

RRNA_GEN_TO_COLOR = {
	7: poster_colors["light_gray"],
	8: DEFAULT_COLOR,
	9: DEFAULT_COLOR,
	10: poster_colors["poster_red"]
}

# FOR SHERLOCK
VARIANT_INDEX_A = 1
VARIANT_INDEX_B = -1 # Optional, set to -1 if only one variant is desired

# # FOR LOCAL DEVELOPMENT
# VARIANT_INDEX_A = 21
# VARIANT_INDEX_B = -1 # Optional, set to -1 if only one variant is desired

std_linewidth = 3


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

			total_plots = 8 # TODO: Modularize and get rid of this magic number
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.rcParams['font.size'] = 20
			fig= plt.figure(figsize=(18, total_plots * 3))

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
			initial_time_absolute_A = time_A[0]
			time_A = time_A - initial_time_absolute_A

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
				GEN_TO_COLOR = GFP_INDUCED_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',
					linewidth=std_linewidth)
			plt.ylabel(r'$\it{gfp}$ mRNA Counts')
			ax1.spines["bottom"].set_position(("outward", 10))
			ax1.spines["left"].set_position(("outward", 10))
			ax1.spines["bottom"].set_visible(False)
			ax1.get_xaxis().set_visible(False)
			max_y = 1500
			ax1.set_ylim([0, max_y])
			ax1.set_yticks([0, max_y / 2, max_y])
			plot_num += 1

			# New Gene Protein Counts
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			new_gene_protein_counts_A = read_stacked_columns(
				all_cells_A, 'MonomerCounts',
				'monomerCounts')[:, new_gene_monomer_indexes]
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = new_gene_protein_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				GEN_TO_COLOR = GFP_INDUCED_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',linewidth=std_linewidth)
			plt.ylabel('GFP Counts')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			max_y = 1.5e6
			ax.set_ylim([0, max_y])
			ax.set_yticks([0, max_y / 2, max_y])
			from matplotlib.ticker import ScalarFormatter
			ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
			ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
			plot_num += 1

			# RNAP Counts
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts_A,) = read_stacked_bulk_molecules(
				all_cells_A, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts_A = read_stacked_columns(
				all_cells_A, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			total_rnap_counts_A = inactive_rnap_counts_A + active_rnap_counts_A
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = total_rnap_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				GEN_TO_COLOR = RNAP_RIBO_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}', linewidth=std_linewidth)
			plt.ylabel('RNA Polymerase Counts')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			max_y = 7500
			ax.set_ylim([0, max_y])
			ax.set_yticks([0, max_y / 2, max_y])
			plot_num += 1

			# Ribosome Counts
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s_A, complex_counts_50s_A) = read_stacked_bulk_molecules(
				all_cells_A, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts_A = np.minimum(
				complex_counts_30s_A, complex_counts_50s_A)
			# Active
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			unique_molecule_counts_table.close()
			active_ribosome_counts_A = read_stacked_columns(
				all_cells_A, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# Total
			total_ribosome_counts_A = inactive_ribosome_counts_A + active_ribosome_counts_A
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = total_ribosome_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				GEN_TO_COLOR = RNAP_RIBO_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}', linewidth=std_linewidth)
			plt.ylabel('Ribosome Counts')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			max_y = 32000
			ax.set_ylim([0, max_y])
			ax.set_yticks([0, max_y / 2, max_y])
			plot_num += 1

			# Cell Mass
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			cell_mass_A = read_stacked_columns(
				all_cells_A, 'Mass', 'cellMass', ignore_exception=True).squeeze()
			critical_initiation_mass_A = read_stacked_columns(
				all_cells_A, 'ReplicationData', 'criticalInitiationMass',
				ignore_exception=True).squeeze()
			plt.axhline(y=np.mean(critical_initiation_mass_A), color='k', linestyle='--', alpha=0.5)
			plt.axhline(y=2 * np.mean(critical_initiation_mass_A), color='k', linestyle='--', alpha=0.5)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = cell_mass_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				GEN_TO_COLOR = CELL_MASS_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}', linewidth=std_linewidth)
			plt.ylabel('Cell Mass (fg)')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			max_y = 2500
			min_y = 750
			ax.set_ylim([min_y, max_y])
			ax.set_yticks([min_y, (max_y - min_y) / 2 + min_y, max_y])
			plot_num += 1

			# Pairs of Forks (Chromosome Replication Initiation)
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				# Because of dynamic sizing, read_stacked_columns does not work for fork_coordinates
				sim_dir_A = all_cells_A[rel_gen_index]
				simOutDirA = os.path.join(sim_dir_A, 'simOut')
				time_data = TableReader(os.path.join(simOutDirA, "Main")).readColumn(
					"time")
				time_data = time_data - initial_time_absolute_A
				replication_data_file_A = TableReader(os.path.join(simOutDirA, "ReplicationData"))
				fork_coordinates_A = replication_data_file_A.readColumn("fork_coordinates")
				pairs_of_forks_A = np.logical_not(np.isnan(fork_coordinates_A)).sum(axis=1) / 2
				replication_data_file_A.close()
				counts_data = pairs_of_forks_A
				GEN_TO_COLOR = FORKS_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}', linewidth=std_linewidth)
			plt.ylabel('Rounds of Active\nChromosome Replication')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			max_y = 2
			ax.set_ylim([0, max_y])
			ax.set_yticks([0, max_y / 2, max_y])
			plot_num += 1

			# Doubling Time
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			dt_to_plot_A = np.repeat(dt_A, num_time_steps_A)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = dt_to_plot_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				GEN_TO_COLOR = DT_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}', linewidth=std_linewidth)
			plt.ylabel('Doubling Time (min)')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			max_y = 62
			min_y = 46
			ax.set_ylim([min_y, max_y])
			ax.set_yticks([min_y, (max_y - min_y) / 2 + min_y, max_y])
			plot_num += 1

			# rRNA Copy Number
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Get rna_ids attribute from RnaSynthProb table in reference cell path
			reference_cell_path = all_cells_A[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')
			# Get indexes of rRNAs in RnaSynthProb table
			transcription = sim_data.process.transcription
			rna_id_to_is_rRNA = {
				rna['id']: rna['is_rRNA'] for rna in transcription.rna_data
				}
			rrna_indexes_rna_synth_prob = np.array([
				i for (i, rna_id) in enumerate(rna_ids)
				if rna_id_to_is_rRNA[rna_id]
				])
			rna_synth_prob_reader.close()
			rrna_copy_numbers_A = read_stacked_columns(
				all_cells_A, 'RnaSynthProb', 'promoter_copy_number',
				fun=lambda x: x[:, rrna_indexes_rna_synth_prob])
			rrna_copy_numbers_A = rrna_copy_numbers_A.sum(axis=1)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = rrna_copy_numbers_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				GEN_TO_COLOR = RRNA_GEN_TO_COLOR
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}', linewidth=std_linewidth)
			plt.ylabel('rRNA Operon Copy Number')
			plt.xlabel('Time (minutes)')
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			max_y = 29
			min_y = 7
			ax.set_ylim([min_y, max_y])
			ax.set_yticks([min_y, (max_y - min_y) / 2 + min_y, max_y])
			max_x = (np.ceil(time_A[-1] / 60 / 20) * 20).item()
			min_x = 0
			ax.set_xlim([min_x, max_x])
			ax.set_xticks([
				min_x, (max_x - min_x) / 4 + min_x,
				2 * (max_x - min_x) / 4 + min_x,
				3 * (max_x - min_x) / 4 + min_x,
				max_x])
			plot_num += 1

			# Disable line clipping
			for ax in fig.get_axes():
				for artist in ax.get_children():
					try:
						artist.set_clip_on(False)
					except AttributeError:
						pass
			# Save figure
			fig.subplots_adjust(hspace=0.5)
			print(f"\nSeed: {seed_index}")
			print("Total number of plots made: ", plot_num - 1)
			# plt.subplots_adjust(hspace=0.7, top=0.95, bottom=0.05)
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
