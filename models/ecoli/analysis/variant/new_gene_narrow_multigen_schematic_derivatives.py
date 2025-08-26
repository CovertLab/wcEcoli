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

START_GEN_INDEX = 0
END_GEN_INDEX = 24 # Not inclusive
GEN_RANGE = np.arange(START_GEN_INDEX, END_GEN_INDEX)

SELECTED_SEED_INDEXES = [
			0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

#GEN_TO_COLOR = {
#	7: poster_colors["light_gray"],
#	8: poster_colors["poster_green"],
#	9: poster_colors["poster_blue"],
#	10: poster_colors["poster_purple"],
#}

GEN_TO_COLOR = {
	0: (0, 0, 0),
	1: (0, 0, 0),
	2: (0, 0, 0),
	3: (0, 0, 0),
	4: (0, 0, 0),
	5: (0, 0, 0),
	6: (0, 0, 0),
	7: (0, 0, 0),
	8: poster_colors["poster_blue"],
	9: poster_colors["poster_blue"],
	10: poster_colors["poster_blue"],
	11: poster_colors["poster_blue"],
	12: poster_colors["poster_blue"],
	13: poster_colors["poster_blue"],
	14: poster_colors["poster_blue"],
	15: poster_colors["poster_blue"],
	16: poster_colors["poster_green"],
	17: poster_colors["poster_green"],
	18: poster_colors["poster_green"],
	19: poster_colors["poster_green"],
	20: poster_colors["poster_green"],
	21: poster_colors["poster_green"],
	22: poster_colors["poster_green"],
	23: poster_colors["poster_green"],
	}

VARIANT_B_COLOR = poster_colors["light_gray"]

# FOR SHERLOCK
VARIANT_INDEX_A = 2
VARIANT_INDEX_B = -1 # Optional, set to -1 if only one variant is desired

# # FOR LOCAL DEVELOPMENT
# VARIANT_INDEX_A = 21
# VARIANT_INDEX_B = 0 # Optional, set to -1 if only one variant is desired

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
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
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
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv. New Gene mRNA Counts')
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
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
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
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv.New Gene Protein Counts')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Limiting RNAP subunit mRNA counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Mapping from cistron ids to indexes
			sim_dir = all_cells_A[0]
			simOutDir = os.path.join(sim_dir, 'simOut')
			mRNA_counts_reader = TableReader(os.path.join(
				simOutDir, 'RNACounts'))
			mRNA_cistron_idx_dict = {rna: i for i, rna in enumerate(
				mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
			mRNA_counts_reader.close()
			limiting_RNAP_subunit_mRNA_cistron_id = "EG10895_RNA"
			limiting_RNAP_subunit_mRNA_cistron_index = mRNA_cistron_idx_dict.get(
				limiting_RNAP_subunit_mRNA_cistron_id)
			limiting_RNAP_subunit_mRNA_counts_A = read_stacked_columns(
				all_cells_A, 'RNACounts',
				'mRNA_cistron_counts')[:, limiting_RNAP_subunit_mRNA_cistron_index]
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = limiting_RNAP_subunit_mRNA_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				limiting_RNAP_subunit_mRNA_counts_B = read_stacked_columns(
					all_cells_B, 'RNACounts',
					'mRNA_cistron_counts')[:, limiting_RNAP_subunit_mRNA_cistron_index]
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = limiting_RNAP_subunit_mRNA_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv.Limiting RNAP Subunit (rpoC) mRNA Counts', fontsize=8)
			plt.xlabel('Time (minutes)')
			plt.ylim(-0.03, 0.03)
			plot_num += 1

			# Limiting RNAP subunit monomer counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			limiting_RNAP_subunit_monomer_id = ["RPOC-MONOMER[c]"]
			limiting_RNAP_subunit_monomer_index = monomer_idx_dict.get(
				limiting_RNAP_subunit_monomer_id[0])
			limiting_RNAP_subunit_monomer_counts_A = read_stacked_columns(
				all_cells_A, 'MonomerCounts',
				'monomerCounts')[:, limiting_RNAP_subunit_monomer_index]
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = limiting_RNAP_subunit_monomer_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				limiting_RNAP_subunit_monomer_counts_B = read_stacked_columns(
					all_cells_B, 'MonomerCounts',
					'monomerCounts')[:, limiting_RNAP_subunit_monomer_index]
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = limiting_RNAP_subunit_monomer_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv.Limiting RNAP Subunit (RpoC) Monomer Counts', fontsize=8)
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# RNAP Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
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
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				(inactive_rnap_counts_B,) = read_stacked_bulk_molecules(
					all_cells_B, (rnap_id,), ignore_exception=True)
				active_rnap_counts_B = read_stacked_columns(
					all_cells_B, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts',
					ignore_exception=True)[:, active_rnap_index]
				total_rnap_counts_B = inactive_rnap_counts_B + active_rnap_counts_B
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = total_rnap_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv.RNAP Counts')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Total Transcriptional Output - NTPs Used
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			count_NTPs_used_A = read_stacked_columns(
				all_cells_A, "TranscriptElongationListener", "countNTPsUSed",
				ignore_exception=True)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = count_NTPs_used_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:]
				# Skip first 10 data points to avoid artifacts from division 
				counts_data = counts_data[10:]
				time_data = time_data[10:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				count_NTPs_used_B = read_stacked_columns(
					all_cells_B, "TranscriptElongationListener", "countNTPsUSed",
					ignore_exception=True)
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = count_NTPs_used_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv. Total Transcriptional Output (NTPs Used)', fontsize=8)
			plt.xlabel('Time (minutes)')
			plt.ylim(-50,50)
			plot_num += 1

			# Limiting rRNA Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			limiting_rRNA_cistron_id = sim_data.molecule_groups.s50_23s_rRNA[0]
			limiting_ribosomal_subunit_id = sim_data.molecule_ids.s50_full_complex
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			unique_molecule_counts_table.close()
			(limiting_rRNA_counts_A) = read_stacked_bulk_molecules(
				all_cells_A, ([limiting_rRNA_cistron_id],), ignore_exception=True)
			(limiting_ribosomal_subunit_counts_A) = read_stacked_bulk_molecules(
				all_cells_A, ([limiting_ribosomal_subunit_id],), ignore_exception=True)
			full_ribosome_counts_A = (read_stacked_columns(
				all_cells_A, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, ribosome_index])
			limiting_rRNA_counts_A = limiting_rRNA_counts_A[0]
			limiting_rRNA_counts_A += limiting_ribosomal_subunit_counts_A[0]
			limiting_rRNA_counts_A += full_ribosome_counts_A
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = limiting_rRNA_counts_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:]
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved 
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				(limiting_rRNA_counts_B) = read_stacked_bulk_molecules(
					all_cells_B, ([limiting_rRNA_cistron_id],), ignore_exception=True)
				(limiting_ribosomal_subunit_counts_B) = read_stacked_bulk_molecules(
					all_cells_B, ([limiting_ribosomal_subunit_id],), ignore_exception=True)
				full_ribosome_counts_B = read_stacked_columns(
					all_cells_B, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts',
					ignore_exception=True)[:, ribosome_index]
				limiting_rRNA_counts_B = limiting_rRNA_counts_B[0]
				limiting_rRNA_counts_B += limiting_ribosomal_subunit_counts_B[0]
				limiting_rRNA_counts_B += full_ribosome_counts_B
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = limiting_rRNA_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv. Limiting rRNA (23s) Counts', fontsize=8)
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Ribosome Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
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
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:] 
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				(complex_counts_30s_B, complex_counts_50s_B) = read_stacked_bulk_molecules(
					all_cells_B, (complex_id_30s, complex_id_50s), ignore_exception=True)
				inactive_ribosome_counts_B = np.minimum(
					complex_counts_30s_B, complex_counts_50s_B)
				active_ribosome_counts_B = read_stacked_columns(
					all_cells_B, 'UniqueMoleculeCounts',
					'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
				total_ribosome_counts_B = inactive_ribosome_counts_B + active_ribosome_counts_B
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = total_ribosome_counts_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv. Ribosome Counts')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Total Translational Output - counts of translated amino acids
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			translated_AAs_A = read_stacked_columns(
				all_cells_A, 'RibosomeData', 'actualElongations',
				ignore_exception=True).squeeze()
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = translated_AAs_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:]
				# Skip first 10 data points to avoid artifacts from division 
				counts_data = counts_data[10:]
				time_data = time_data[10:]  
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				translated_AAs_B = read_stacked_columns(
					all_cells_B, 'RibosomeData', 'actualElongations',
					ignore_exception=True).squeeze()
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = translated_AAs_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv. Total Translational Output (AAs)')
			plt.xlabel('Time (minutes)')
			plt.ylim(-300,300)
			plot_num += 1

			# Cell Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			cell_mass_A = read_stacked_columns(
				all_cells_A, 'Mass', 'cellMass', ignore_exception=True).squeeze()
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = cell_mass_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				# Use the derivatives for plotting
				time_data = time_data.squeeze()
				counts_data = counts_data.squeeze()
				counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
				time_data = time_data[1:]
				# Skip first 10 data points to avoid artifacts from division
				counts_data = counts_data[10:]
				time_data = time_data[10:]
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				counts_data = counts_data_convolved 
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				cell_mass_B = read_stacked_columns(
					all_cells_B, 'Mass', 'cellMass', ignore_exception=True).squeeze()
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = cell_mass_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					# Use the derivatives for plotting
					time_data = time_data.squeeze()
					counts_data = counts_data.squeeze()
					counts_data = (counts_data[1:] - counts_data[:-1]) / (time_data[1:] - time_data[:-1])
					time_data = time_data[1:] 
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Deriv. Cell Mass (fg)')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Instantaneous Doubling Time
			ax = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			instantaneous_growth_rate_A = read_stacked_columns(
				all_cells_A, 'Mass', 'instantaneous_growth_rate',
				ignore_exception=True).squeeze()
			instantaneous_dt_A = np.log(2) / instantaneous_growth_rate_A / 60.0
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = instantaneous_dt_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				time_data = time_data[100:]
				counts_data = counts_data[100:]
				moving_window = min(301, len(counts_data))
				convolution_array = (np.ones(moving_window) / moving_window)
				pad = moving_window // 2
				padded = np.pad(counts_data, pad_width=pad, mode='edge')
				counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
				plt.plot(
					time_data / 60., counts_data_convolved,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}')
			if VARIANT_INDEX_B != -1:
				instantaneous_growth_rate_B = read_stacked_columns(
					all_cells_B, 'Mass', 'instantaneous_growth_rate',
					ignore_exception=True).squeeze()
				instantaneous_dt_B = np.log(2) / instantaneous_growth_rate_B / 60.0
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = instantaneous_dt_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					time_data = time_data[100:]
					counts_data = counts_data[100:]
					moving_window = min(301, len(counts_data))
					convolution_array = (np.ones(moving_window) / moving_window)
					pad = moving_window // 2
					padded = np.pad(counts_data, pad_width=pad, mode='edge')
					counts_data_convolved = np.convolve(padded, convolution_array, mode='valid')
					plt.plot(
						time_data / 60., counts_data_convolved,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}')
			plt.ylabel('Instantaneous Doubling Time (min)')
			# max_y = 80
			# min_y = 40
			# ax.set_ylim([min_y, max_y])
			# ax.set_yticks([min_y, (max_y - min_y) / 2 + min_y, max_y])
			plot_num += 1

			# Critical Chromosome Replication Initiation Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			critical_initiation_mass_A = read_stacked_columns(
				all_cells_A, 'ReplicationData', 'criticalInitiationMass',
				ignore_exception=True).squeeze()
			critical_mass_equivalents_A = (
				cell_mass_A / critical_initiation_mass_A)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = critical_mass_equivalents_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				critical_initiation_mass_B = read_stacked_columns(
					all_cells_B, 'ReplicationData', 'criticalInitiationMass',
					ignore_exception=True).squeeze()
				critical_mass_equivalents_B = (
					cell_mass_B / critical_initiation_mass_B)
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = critical_mass_equivalents_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.axhline(y=1, color='k', linestyle='--')
			plt.axhline(y=2, color='k', linestyle='--')
			plt.ylabel('Factors of Critical Initiation Mass')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Pairs of Forks (Chromosome Replication Initiation)
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				# Because of dynamic sizing, read_stacked_columns does not work for fork_coordinates
				sim_dir_A = all_cells_A[rel_gen_index]
				simOutDirA = os.path.join(sim_dir_A, 'simOut')
				time_data = TableReader(os.path.join(simOutDirA, "Main")).readColumn(
					"time")
				replication_data_file_A = TableReader(os.path.join(simOutDirA, "ReplicationData"))
				fork_coordinates_A = replication_data_file_A.readColumn("fork_coordinates")
				pairs_of_forks_A = np.logical_not(np.isnan(fork_coordinates_A)).sum(axis=1) / 2
				replication_data_file_A.close()
				counts_data = pairs_of_forks_A
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					# Because of dynamic sizing, read_stacked_columns does not work for fork_coordinates
					sim_dir_B = all_cells_B[rel_gen_index]
					simOutDirB = os.path.join(sim_dir_B, 'simOut')
					time_data = TableReader(os.path.join(simOutDirB, "Main")).readColumn(
						"time")
					replication_data_file_B = TableReader(os.path.join(simOutDirB, "ReplicationData"))
					fork_coordinates_B = replication_data_file_B.readColumn("fork_coordinates")
					pairs_of_forks_B = np.logical_not(np.isnan(fork_coordinates_B)).sum(axis=1) / 2
					replication_data_file_B.close()
					counts_data = pairs_of_forks_B
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Pairs of Forks')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# Doubling Time
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			dt_to_plot_A = np.repeat(dt_A, num_time_steps_A)
			for gen in GEN_RANGE:
				rel_gen_index = gen - START_GEN_INDEX
				time_data = time_A[gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				counts_data = dt_to_plot_A[
					gen_start_indexes_A[rel_gen_index]:gen_end_indexes_A[rel_gen_index] + 1]
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				dt_to_plot_B = np.repeat(dt_B, num_time_steps_B)
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = dt_to_plot_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Doubling Time (min)')
			plt.xlabel('Time (minutes)')
			plot_num += 1

			# rRNA Copy Number
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
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
				plt.plot(
					time_data / 60., counts_data,
					color=GEN_TO_COLOR[gen],
					label = f'Variant {VARIANT_INDEX_A}',)
			if VARIANT_INDEX_B != -1:
				rrna_copy_numbers_B = read_stacked_columns(
					all_cells_B, 'RnaSynthProb', 'promoter_copy_number',
					fun=lambda x: x[:, rrna_indexes_rna_synth_prob])
				rrna_copy_numbers_B = rrna_copy_numbers_B.sum(axis=1)
				for gen in GEN_RANGE:
					rel_gen_index = gen - START_GEN_INDEX
					time_data = time_B[gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					counts_data = rrna_copy_numbers_B[
						gen_start_indexes_B[rel_gen_index]:gen_end_indexes_B[rel_gen_index] + 1]
					plt.plot(
						time_data / 60., counts_data,
						color=VARIANT_B_COLOR,
						label = f'Variant {VARIANT_INDEX_B}',)
			plt.ylabel('Summed rRNA Copy Number')
			plt.xlabel('Time (minutes)')
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
