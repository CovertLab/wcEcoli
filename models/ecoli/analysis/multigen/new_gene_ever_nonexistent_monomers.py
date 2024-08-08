"""
Plot ...
"""
# TODO: update file header comment

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
# noinspection PyUnresolvedReferences
import numpy as np
from numpy import inf

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.sim.variants.new_gene_internal_shift import determine_new_gene_ids_and_indices
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

# START_GEN = 0
# END_GEN = 3

START_GEN = 17
END_GEN = 24

LINE_COLOR = (66/255, 170/255, 154/255)

cistron_of_interest = "EG10945_RNA"
gene_name = "serB" # TODO: don't hardcode this
reaction_ids_of_interest = ['RXN0-5114[CCO-CYTOSOL]-3-P-SERINE/WATER//SER/Pi.38.'] # TODO: don't hardcode this
AA_of_interest = "SER[c]" # TODO: don't hardcode this

# 'RXN0-5114[CCO-CYTOSOL]-3-P-SERINE/WATER//SER/Pi.38.', 'RXN0-5114[CCO-PERI-BAC]-3-P-SERINE/WATER//SER/Pi.39.',

# EG10945 EG10945_RNA	PSERPHOSPHA-MONOMER

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
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
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells(generation=np.arange(START_GEN, END_GEN))
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array


		gene_of_interest_cistron_ids = [cistron_of_interest]

		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		gene_of_interest_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in gene_of_interest_cistron_ids]

		# TODO: make a option for this plot that you could run without new genes

		# Extract mRNA indexes for each gene of interest
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		gene_of_interest_cistron_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
								 gene_of_interest_cistron_ids]
		mRNA_counts_reader.close()

		# Extract protein indexes for each gene of interest
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}
		gene_of_interest_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									monomer_id in gene_of_interest_monomer_ids]
		monomer_counts_reader.close()

		# Extract reaction index for each reaction of interest
		fba_results_reader = TableReader(os.path.join(simOutDir, 'FBAResults'))
		reaction_ids = fba_results_reader.readAttribute('reactionIDs')
		reaction_idx_dict = {reaction: i for i, reaction in enumerate(reaction_ids)}
		reaction_of_interest_indexes = [
			reaction_idx_dict.get(reaction_id) for reaction_id in reaction_ids_of_interest]
		fba_results_reader.close()

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)
		time_no_first = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)

		time = time - time[0]
		time_no_first = time_no_first - time_no_first[0]

		all_monomer_stacked_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)
		gene_of_interest_monomer_counts = all_monomer_stacked_counts[:,gene_of_interest_monomer_indexes]
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts', ignore_exception=True)
		gene_of_interest_mRNA_counts = all_mRNA_stacked_counts[:,gene_of_interest_cistron_indexes]

		all_rection_fluxes = read_stacked_columns(
			cell_paths, 'FBAResults', 'reactionFluxes', ignore_exception=True)
		reaction_of_interest_flux = all_rection_fluxes[:,reaction_of_interest_indexes]

		(aa_of_interest_counts, ) = read_stacked_bulk_molecules(
			cell_paths, [AA_of_interest], ignore_exception=True)

		plot_suffixes = [""]
		standard_xlim = (0,1500)
		total_plots = 7 # TODO Modularize and get rid of this magic number

		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			plt.figure(figsize = (8.5, 11))
			plot_num = 1

			# Mass
			ax1 = plt.subplot(total_plots, 1, plot_num)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			plt.plot(time / 60., mass, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,4000)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Cell Mass (fg)", fontsize="small")
			plot_num += 1

			# TODO: transcription event

			# mRNA Cistron Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if plot_suffix == "":
				if len(gene_of_interest_cistron_ids) == 1:
					plt.plot(time / 60., gene_of_interest_mRNA_counts, color=LINE_COLOR)
				else:
					for r in range(len(gene_of_interest_cistron_ids)):
						plt.plot(time / 60., gene_of_interest_mRNA_counts[:,r],
								 label = gene_of_interest_cistron_ids[r], color=LINE_COLOR)
					plt.legend()
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# plot on log scale instead
				if len(gene_of_interest_cistron_ids) == 1:
					plt.plot(time / 60., np.log10(gene_of_interest_mRNA_counts + 1), color=LINE_COLOR)
				else:
					for r in range(len(gene_of_interest_cistron_ids)):
						plt.plot(time / 60., np.log10(gene_of_interest_mRNA_counts[:,r] + 1),
								 label = gene_of_interest_cistron_ids[r], color=LINE_COLOR)
					plt.legend()
				plt.ylim((-1,4.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA Cistron Counts: " + gene_name, fontsize="small")
			plot_num += 1

			# Protein Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if plot_suffix == "":
				if len(gene_of_interest_monomer_ids) == 1:
					plt.plot(time / 60., gene_of_interest_monomer_counts, color=LINE_COLOR)
				else:
					for m in range(len(gene_of_interest_monomer_ids)):
						plt.plot(time / 60., gene_of_interest_monomer_counts[:,m],
								 label = gene_of_interest_monomer_ids[m], color=LINE_COLOR)
					plt.legend()
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# plot on log scale instead
				if len(gene_of_interest_monomer_ids) == 1:
					plt.plot(time / 60., np.log10(gene_of_interest_monomer_counts + 1), color=LINE_COLOR)
				else:
					for m in range(len(gene_of_interest_monomer_ids)):
						plt.plot(time / 60., np.log10(gene_of_interest_monomer_counts[:,m] + 1),
								 label = gene_of_interest_monomer_ids[m], color=LINE_COLOR)
					plt.legend()
				plt.ylim((-1,7.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Monomer Counts: " + gene_name, fontsize="small")
			plot_num += 1

			# TODO: reaction flux
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.plot(time / 60., reaction_of_interest_flux, color=LINE_COLOR)
			plt.xlabel("Time (min)")
			plt.ylabel("Reaction flux (mmol/g DCW/hr)", fontsize="small")
			plot_num += 1


			# TODO: reaction product
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.plot(time / 60., aa_of_interest_counts, color=LINE_COLOR)
			plt.xlabel("Time (min)")
			plt.ylabel("" + AA_of_interest + " counts", fontsize="small")
			plot_num += 1

			# TODO: translation stalling

			# Growth Rate
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(150, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')
			plt.plot(time.flatten() / 60., growth_rate_convolved, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,.0004)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Growth Rate", fontsize="small")
			plot_num += 1

			# ppGpp
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			ppGpp_concentration = read_stacked_columns(
				cell_paths, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)
			plt.plot(time_no_first / 60., ppGpp_concentration, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,150))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("ppGpp Concentration ($\mu$M)", fontsize="small")
			plot_num += 1

			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(
				plt, plotOutDir, plotOutFileName + plot_suffix + "_"
				+ str(START_GEN) + "_" + str(END_GEN), metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
