"""
Plot mRNA and protein counts for new genes across multiple generations, as well
as plots to analyze the impact of new gene expression, including growth rate,
RNAP and ribosome counts, and ppGpp concentration.
"""
# TODO: update file header comment

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib as mpl
# noinspection PyUnresolvedReferences
import numpy as np
from numpy import inf

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# START_GEN = 0
# END_GEN = 3

START_GEN = 16
END_GEN = 24

LINE_COLOR = (66/255, 170/255, 154/255)
LINE_COLOR2 = (27/255, 132/255, 198/255)

POSTER_VARIANT_COLORS = [ (136/255, 205/255, 240/255),
						  (188/255, 140/255, 191/255),
						  (66/255, 170/255, 154/255),
						  (221/255, 203/255, 119/255),
						  (27/255, 132/255, 198/255)]

LINE_WIDTH = 0.7



class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def get_mRNA_ids_from_monomer_ids(self, sim_data, target_monomer_ids):
		"""
		Map monomer ids back to the mRNA ids that they were translated from.

		Args:
			target_monomer_ids: ids of the monomers to map to mRNA ids

		Returns: set of mRNA ids
		"""
		# Map protein ids to cistron ids
		monomer_ids = sim_data.process.translation.monomer_data['id']
		cistron_ids = sim_data.process.translation.monomer_data[
			'cistron_id']
		monomer_to_cistron_id_dict = {
			monomer_id: cistron_ids[i] for i, monomer_id in
			enumerate(monomer_ids)}
		target_cistron_ids = [
			monomer_to_cistron_id_dict.get(RNAP_monomer_id) for
			RNAP_monomer_id in target_monomer_ids]
		# Map cistron ids to RNA indexes
		target_RNA_indexes = [
			sim_data.process.transcription.cistron_id_to_rna_indexes(
				RNAP_cistron_id) for RNAP_cistron_id in
			target_cistron_ids]
		# Map RNA indexes to RNA ids
		RNA_ids = sim_data.process.transcription.rna_data['id']
		target_RNA_ids = set()
		for i in range(len(target_RNA_indexes)):
			for index in target_RNA_indexes[i]:
				target_RNA_ids.add(RNA_ids[index])
		return target_RNA_ids

	def get_mRNA_indexes_from_monomer_ids(self, sim_data, all_cells, target_monomer_ids, index_type):
		"""
		Retrieve new gene indexes of a given type.

		Args:
			all_cells: Paths to all cells to read data from (directories should
				contain a simOut/ subdirectory), typically the return from
				AnalysisPaths.get_cells()
			target_monomer_ids: ids of the monomers to map to mRNA indexes
			index_type: Type of indexes to extract, currently supported options
				are 'mRNA' and 'monomer'

		Returns:
			List of requested indexes
		"""
		sim_dir = all_cells[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		if index_type == 'mRNA':
			# Map protein ids to RNA ids
			target_RNA_ids = self.get_mRNA_ids_from_monomer_ids(
				sim_data, target_monomer_ids)
			# Get index of those RNA ids in the output
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_idx_dict = {
				rna: i for i, rna in
				enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}
			output_indexes = [
				mRNA_idx_dict.get(mRNA_id) for mRNA_id in target_RNA_ids]

		elif index_type == 'monomer':
			# Get index of those monomer ids in the output
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, "MonomerCounts"))
			monomer_idx_dict = {
				monomer: i for i, monomer in enumerate(
					monomer_counts_reader.readAttribute('monomerIds'))}
			output_indexes = [
				monomer_idx_dict.get(monomer_id) for monomer_id in
				target_monomer_ids]

		else:
			raise Exception(
				"Index type " + index_type +
				" has no instructions for data extraction.")

		return output_indexes

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		cell_paths = self.ap.get_cells(
			variant=[0], seed=[0],
			generation=np.arange(START_GEN, END_GEN))
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

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

		# Extract mRNA indexes for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_ids'))}
		new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
								 new_gene_mRNA_ids]
		mRNA_counts_reader.close()

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}
		new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									monomer_id in new_gene_monomer_ids]
		monomer_counts_reader.close()

		# Extract RNA indexes for each new gene
		rnap_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))
		RNA_idx_dict = {
			rna[:-3]: i for i, rna in
			enumerate(rnap_reader.readAttribute('rnaIds'))}
		new_gene_RNA_indexes = [
			RNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]
		rnap_reader.close()

		# Load data
		# time = read_stacked_columns(
		# 	cell_paths, 'Main', 'time', ignore_exception=True)
		# time_no_first = read_stacked_columns(
		# 	cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
		#
		# time = time - time[0]
		# time_no_first = time_no_first - time_no_first[0]

		# (new_gene_monomer_counts,) = read_stacked_bulk_molecules(
		# 	cell_paths, new_gene_monomer_ids, ignore_exception=True)
		# all_mRNA_stacked_counts = read_stacked_columns(
		# 	cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)
		# new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]

		plot_suffix = ""

		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False
		# mpl.spines["bottom"].set_position(("outward", 10))
		# mpl.spines["left"].set_position(("outward", 10))

		# Plotting
		max_x = 300
		standard_xlim = (0, max_x)

		# # Mass
		# VARIANTS_TO_PLOT = [0, 1, 2]
		# SEEDS_TO_PLOT = [0, 0, 0]
		# assert len(VARIANTS_TO_PLOT) == len(SEEDS_TO_PLOT)
		# plt.figure(figsize=(6, 3))
		# total_plots = len(VARIANTS_TO_PLOT)
		# plot_name = "mass"
		# plot_num = 1
		# for i in range(len(VARIANTS_TO_PLOT)):
		# 	curr_var = VARIANTS_TO_PLOT[i]
		# 	curr_seed = SEEDS_TO_PLOT[i]
		# 	print(curr_var, curr_seed)
		# 	if i == 0:
		# 		ax1 = plt.subplot(total_plots, 1, plot_num)
		# 	if i == 1:
		# 		ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
		# 	if i == 2:
		# 		ax3 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
		#
		# 	cell_paths = self.ap.get_cells(
		# 		variant=[curr_var], seed=[curr_seed],
		# 		generation=np.arange(START_GEN, END_GEN))
		#
		# 	time = read_stacked_columns(
		# 		cell_paths, 'Main', 'time', ignore_exception=True)
		# 	time_no_first = read_stacked_columns(
		# 		cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
		# 	time = time - time[0]
		# 	time_no_first = time_no_first - time_no_first[0]
		#
		# 	mass = read_stacked_columns(
		# 		cell_paths, "Mass", "cellMass", ignore_exception=True)
		#
		# 	plt.plot(time / 60., mass, color=LINE_COLOR, clip_on=False)
		# 	plt.xlim(standard_xlim)
		#
		# 	if i == 0:
		# 		ax1.spines["bottom"].set_position(("outward", 10))
		# 		ax1.spines["left"].set_position(("outward", 10))
		# 		ax1.spines["bottom"].set_visible(False)
		# 		ax1.get_xaxis().set_visible(False)
		# 		max_y = 2600
		# 		ax1.set_ylim([0, max_y])
		# 		ax1.set_yticks([0, max_y /2, max_y])
		# 	elif i == 1:
		# 		ax2.spines["bottom"].set_position(("outward", 10))
		# 		ax2.spines["left"].set_position(("outward", 10))
		# 		ax2.spines["bottom"].set_visible(False)
		# 		ax2.get_xaxis().set_visible(False)
		# 		max_y = 2600
		# 		ax2.set_ylim([0, max_y])
		# 		ax2.set_yticks([0, max_y /2, max_y])
		# 	else:
		# 		ax3.spines["bottom"].set_position(("outward", 10))
		# 		ax3.spines["left"].set_position(("outward", 10))
		# 		max_y = 2600
		# 		ax3.set_ylim([0, max_y])
		# 		ax3.set_yticks([0, max_y /2, max_y])
		# 		ax3.set_xticks([0, 150, 300])
		#
		# 	plot_num += 1
		#
		# # plt.xlabel("Time (min)")
		# # plt.ylabel("Cell Mass (fg)", fontsize="small")
		# plt.tight_layout()
		# exportFigure(
		# 	plt, plotOutDir,
		# 	plotOutFileName + plot_suffix + "_" + plot_name + "_"
		# 		+ str(START_GEN) + "_" + str(END_GEN), metadata)
		# plt.close("all")


		# mRNA Counts
		VARIANTS_TO_PLOT = [0, 2, 7]
		SEEDS_TO_PLOT = [1, 3, 2]
		assert len(VARIANTS_TO_PLOT) == len(SEEDS_TO_PLOT)
		plt.figure(figsize=(6, 3))
		total_plots = len(VARIANTS_TO_PLOT)

		plot_name = "mRNA_counts"
		plot_num = 1
		for i in range(len(VARIANTS_TO_PLOT)):
			curr_var = VARIANTS_TO_PLOT[i]
			curr_seed = SEEDS_TO_PLOT[i]
			print(curr_var, curr_seed)
			if i == 0:
				ax1 = plt.subplot(total_plots, 1, plot_num)
			if i == 1:
				ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if i == 2:
				ax3 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			cell_paths = self.ap.get_cells(
				variant=[curr_var], seed=[curr_seed],
				generation=np.arange(START_GEN, END_GEN))

			time = read_stacked_columns(
				cell_paths, 'Main', 'time', ignore_exception=True)
			time_no_first = read_stacked_columns(
				cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
			time = time - time[0]
			time_no_first = time_no_first - time_no_first[0]

			all_mRNA_stacked_counts = read_stacked_columns(
					cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)
			new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]

			plt.plot(
				time[time / 60. <= max_x] / 60.,
				new_gene_mRNA_counts[time / 60. <= max_x],
				color=LINE_COLOR, clip_on=False, linewidth=LINE_WIDTH)
			plt.xlim(standard_xlim)

			if i == 0:
				ax1.spines["bottom"].set_position(("outward", 10))
				ax1.spines["left"].set_position(("outward", 10))
				ax1.spines["bottom"].set_visible(False)
				ax1.get_xaxis().set_visible(False)
				max_y = 100
				ax1.set_ylim([0, max_y])
				ax1.set_yticks([0, max_y /2, max_y])
			elif i == 1:
				ax2.spines["bottom"].set_position(("outward", 10))
				ax2.spines["left"].set_position(("outward", 10))
				ax2.spines["bottom"].set_visible(False)
				ax2.get_xaxis().set_visible(False)
				max_y = 300
				ax2.set_ylim([0, max_y])
				ax2.set_yticks([0, max_y /2, max_y])
			else:
				ax3.spines["bottom"].set_position(("outward", 10))
				ax3.spines["left"].set_position(("outward", 10))
				max_y = 900
				ax3.set_ylim([0, max_y])
				ax3.set_yticks([0, max_y /2, max_y])
				ax3.set_xticks([0, 150, 300])

			plot_num += 1
		plt.tight_layout()
		exportFigure(
			plt, plotOutDir,
			plotOutFileName + plot_suffix + "_" + plot_name + "_"
				+ str(START_GEN) + "_" + str(END_GEN), metadata)
		plt.close("all")


		# Protein Counts
		# VARIANTS_TO_PLOT = [1, 2]
		# SEEDS_TO_PLOT = [0, 0]
		VARIANTS_TO_PLOT = [7, 6]
		SEEDS_TO_PLOT = [2, 0]
		assert len(VARIANTS_TO_PLOT) == len(SEEDS_TO_PLOT)
		plt.figure(figsize=(6, 3))
		total_plots = 3

		plot_name = "protein_counts"
		plot_num = 1
		for i in range(len(VARIANTS_TO_PLOT)):
			curr_var = VARIANTS_TO_PLOT[i]
			curr_seed = SEEDS_TO_PLOT[i]
			print(curr_var, curr_seed)
			if i == 0:
				ax1 = plt.subplot(total_plots, 1, plot_num)
			if i == 1:
				ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			cell_paths = self.ap.get_cells(
				variant=[curr_var], seed=[curr_seed],
				generation=np.arange(START_GEN, END_GEN))

			time = read_stacked_columns(
				cell_paths, 'Main', 'time', ignore_exception=True)
			time = time - time[0]

			(new_gene_monomer_counts,) = read_stacked_bulk_molecules(
				cell_paths, new_gene_monomer_ids, ignore_exception=True)

			plt.plot(
				time[time / 60. <= max_x] / 60.,
				new_gene_monomer_counts[(time / 60. <= max_x).squeeze()],
				color=LINE_COLOR, clip_on=False, linewidth=LINE_WIDTH)
			plt.xlim(standard_xlim)

			if i == 0:
				ax1.spines["bottom"].set_position(("outward", 10))
				ax1.spines["left"].set_position(("outward", 10))
				ax1.spines["bottom"].set_visible(False)
				ax1.get_xaxis().set_visible(False)
				max_y = 750000
				ax1.set_ylim([0, max_y])
				ax1.set_yticks([0, max_y /2, max_y])
				# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
				ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
			elif i == 1:
				ax2.spines["bottom"].set_position(("outward", 10))
				ax2.spines["left"].set_position(("outward", 10))
				max_y = 1500000
				ax2.set_ylim([0, max_y])
				ax2.set_yticks([0, max_y /2, max_y])
				ax2.set_xticks([0, 150, 300])
				ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
			plot_num += 1
		plt.tight_layout()
		exportFigure(
			plt, plotOutDir,
			plotOutFileName + plot_suffix + "_" + plot_name + "_"
				+ str(START_GEN) + "_" + str(END_GEN), metadata)
		plt.close("all")

		VARIANTS_TO_PLOT = [0, 6]
		SEEDS_TO_PLOT = [1, 0]
		assert len(VARIANTS_TO_PLOT) == len(SEEDS_TO_PLOT)
		plt.figure(figsize=(6, 3))
		total_plots = 3
		plot_name = "ppGpp"
		plot_num = 1
		for i in range(len(VARIANTS_TO_PLOT)):
			curr_var = VARIANTS_TO_PLOT[i]
			curr_seed = SEEDS_TO_PLOT[i]
			print(curr_var, curr_seed)
			if i == 0:
				ax1 = plt.subplot(total_plots, 1, plot_num)
			if i == 1:
				ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			cell_paths = self.ap.get_cells(
				variant=[curr_var], seed=[curr_seed],
				generation=np.arange(START_GEN, END_GEN))

			time_no_first = read_stacked_columns(
				cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
			time_no_first = time_no_first - time_no_first[0]

			ppGpp_concentration = read_stacked_columns(
				cell_paths, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)

			plt.plot(
				time_no_first[time_no_first / 60. <= max_x] / 60.,
				ppGpp_concentration[time_no_first / 60. <= max_x],
				color=LINE_COLOR2,
				clip_on=False, linewidth=LINE_WIDTH)
			plt.xlim(standard_xlim)

			if i == 0:
				ax1.spines["bottom"].set_position(("outward", 10))
				ax1.spines["left"].set_position(("outward", 10))
				ax1.spines["bottom"].set_visible(False)
				ax1.get_xaxis().set_visible(False)
				max_y = 140
				ax1.set_ylim([0, max_y])
				ax1.set_yticks([0, max_y / 2, max_y])
			elif i == 1:
				ax2.spines["bottom"].set_position(("outward", 10))
				ax2.spines["left"].set_position(("outward", 10))
				max_y = 140
				ax2.set_ylim([0, max_y])
				ax2.set_yticks([0, max_y / 2, max_y])
				ax2.set_xticks([0, 150, 300])
			plot_num += 1
		plt.tight_layout()
		exportFigure(
			plt, plotOutDir,
			plotOutFileName + plot_suffix + "_" + plot_name + "_0_6_"
			+ str(START_GEN) + "_" + str(END_GEN), metadata)
		plt.close("all")

		# Growth Rate
		plt.figure(figsize=(6, 3))
		total_plots = 3

		plot_name = "growth_rate"
		plot_num = 1
		for i in range(len(VARIANTS_TO_PLOT)):
			curr_var = VARIANTS_TO_PLOT[i]
			curr_seed = SEEDS_TO_PLOT[i]
			print(curr_var, curr_seed)
			if i == 0:
				ax1 = plt.subplot(total_plots, 1, plot_num)
			if i == 1:
				ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			cell_paths = self.ap.get_cells(
				variant=[curr_var], seed=[curr_seed],
				generation=np.arange(START_GEN, END_GEN))

			time = read_stacked_columns(
				cell_paths, 'Main', 'time', ignore_exception=True)
			time = time - time[0]

			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(150, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')

			plt.plot(
				time[time / 60. <= max_x] / 60.,
				growth_rate_convolved[(time / 60. <= max_x).squeeze()],
				color=LINE_COLOR2,
				clip_on=False, linewidth=LINE_WIDTH)
			plt.xlim(standard_xlim)

			if i == 0:
				ax1.spines["bottom"].set_position(("outward", 10))
				ax1.spines["left"].set_position(("outward", 10))
				ax1.spines["bottom"].set_visible(False)
				ax1.get_xaxis().set_visible(False)
				max_y = 0.0005
				ax1.set_ylim([0, max_y])
				ax1.set_yticks([0, max_y / 2, max_y])
			elif i == 1:
				ax2.spines["bottom"].set_position(("outward", 10))
				ax2.spines["left"].set_position(("outward", 10))
				max_y = 0.0005
				ax2.set_ylim([0, max_y])
				ax2.set_yticks([0, max_y / 2, max_y])
				ax2.set_xticks([0, 150, 300])
			plot_num += 1
		plt.tight_layout()
		exportFigure(
			plt, plotOutDir,
			plotOutFileName + plot_suffix + "_" + plot_name + "_0_6_"
			+ str(START_GEN) + "_" + str(END_GEN), metadata)
		plt.close("all")

if __name__ == '__main__':
	Plot().cli()