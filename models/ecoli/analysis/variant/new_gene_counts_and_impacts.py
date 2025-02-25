"""
Plot mRNA and protein counts for new genes across multiple generations, as well
as plots to analyze the impact of new gene expression, including growth rate,
RNAP and ribosome counts, and ppGpp concentration.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

START_GEN = 0
END_GEN = 24

VARIANT_1 = 3
VARIANT_2 = 6

VARIANT_1_SEED = 9
VARIANT_2_SEED = 9

LINE_COLOR = (66/255, 170/255, 154/255)
LINE_COLOR2 = (152/255, 78/255, 163/255)

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
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells(
			variant=np.array([VARIANT_1]), generation=np.arange(START_GEN, END_GEN),
			seed=np.array([VARIANT_1_SEED]))
		cell_paths2 = self.ap.get_cells(
			variant=np.array([VARIANT_2]), generation=np.arange(START_GEN, END_GEN),
			seed=np.array([VARIANT_2_SEED]))
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
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)
		time_no_first = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
		time2 = read_stacked_columns(
			cell_paths2, 'Main', 'time', ignore_exception=True)
		time_no_first2 = read_stacked_columns(
			cell_paths2, 'Main', 'time', remove_first=True, ignore_exception=True)

		time = time - time[0]
		time_no_first = time_no_first - time_no_first[0]
		time2 = time2 - time2[0]
		time_no_first2 = time_no_first2 - time_no_first2[0]

		(new_gene_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, new_gene_monomer_ids, ignore_exception=True)
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)
		new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]

		(new_gene_monomer_counts2,) = read_stacked_bulk_molecules(
			cell_paths2, new_gene_monomer_ids, ignore_exception=True)
		all_mRNA_stacked_counts2 = read_stacked_columns(
			cell_paths2, 'RNACounts', 'mRNA_counts', ignore_exception=True)
		new_gene_mRNA_counts2 = all_mRNA_stacked_counts2[:,new_gene_mRNA_indexes]

		plot_suffixes = ["", "_standard_axes_y"]
		standard_xlim = (0,400)
		total_plots = 7 # TODO Modularize and get rid of this magic number

		mpl.rcParams['axes.spines.right'] = False
		mpl.rcParams['axes.spines.top'] = False

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			plt.figure(figsize = (8.5, 14))
			plot_num = 1

			# Growth Rate
			ax1 = plt.subplot(total_plots, 1, plot_num)

			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(150, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')
			plt.plot(time.flatten() / 60., growth_rate_convolved, color=LINE_COLOR)

			growth_rate2 = np.ravel(read_stacked_columns(
				cell_paths2, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window2 = min(150, len(growth_rate2))
			convolution_array2 = (np.ones(moving_window2) / moving_window2)
			growth_rate_convolved2 = np.convolve(
				convolution_array2, growth_rate2, mode='same')
			plt.plot(time2.flatten() / 60., growth_rate_convolved2, color=LINE_COLOR2)

			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,.0004)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Growth Rate", fontsize="small")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			plt.plot(time / 60., mass, color=LINE_COLOR, label = "Variant " + str(VARIANT_1))

			mass2 = read_stacked_columns(
				cell_paths2, "Mass", "cellMass", ignore_exception=True)
			plt.plot(time2 / 60., mass2, color=LINE_COLOR2, label = "Variant " + str(VARIANT_2))

			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,4000)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Cell Mass (fg)", fontsize="small")
			plt.legend()
			plot_num += 1

			# mRNA Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if plot_suffix == "":
				plt.plot(time / 60., new_gene_mRNA_counts, color=LINE_COLOR)
				plt.plot(time2 / 60., new_gene_mRNA_counts2, color=LINE_COLOR2)
				plt.ylabel("gfp mRNA Counts", fontsize="small")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.plot(time / 60., np.log10(new_gene_mRNA_counts + 1), color=LINE_COLOR)
				plt.plot(time2 / 60., np.log10(new_gene_mRNA_counts2 + 1), color=LINE_COLOR2)
				plt.ylim((-1,4.5))
				plt.ylabel("Log(gfp mRNA Counts + 1)", fontsize="small")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plot_num += 1

			# Protein Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if plot_suffix == "":
				plt.plot(time / 60., new_gene_monomer_counts, color=LINE_COLOR)
				plt.plot(time2 / 60., new_gene_monomer_counts2, color=LINE_COLOR2)
				plt.ylabel("GFP Counts", fontsize="small")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.plot(time / 60., np.log10(new_gene_monomer_counts + 1), color=LINE_COLOR)
				plt.plot(time2 / 60., np.log10(new_gene_monomer_counts2 + 1), color=LINE_COLOR2)
				plt.ylim((-1,7.5))
				plt.ylabel("Log(GFP Counts + 1)", fontsize="small")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")

			plot_num += 1

			# ppGpp
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)

			ppGpp_concentration = read_stacked_columns(
				cell_paths, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)
			plt.plot(time_no_first / 60., ppGpp_concentration, color=LINE_COLOR)

			ppGpp_concentration2 = read_stacked_columns(
				cell_paths2, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)
			plt.plot(time_no_first2 / 60., ppGpp_concentration2, color=LINE_COLOR2)

			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,150))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("ppGpp Concentration ($\mu$M)", fontsize="small")
			plot_num += 1

			# Total RNAP Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts,) = read_stacked_bulk_molecules(
				cell_paths, (rnap_id,), ignore_exception=True)
			(inactive_rnap_counts2,) = read_stacked_bulk_molecules(
				cell_paths2, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			active_rnap_counts2 = read_stacked_columns(
				cell_paths2, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			# Total
			total_rnap_counts = inactive_rnap_counts + active_rnap_counts
			total_rnap_counts2 = inactive_rnap_counts2 + active_rnap_counts2
			plt.plot(time / 60., total_rnap_counts, color=LINE_COLOR)
			plt.plot(time2 / 60., total_rnap_counts2, color=LINE_COLOR2)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,10000))
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Total RNAP Counts", fontsize="small")
			plot_num += 1

			# Total Ribosome Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
				cell_paths, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts = np.minimum(
				complex_counts_30s, complex_counts_50s)

			(complex_counts_30s2, complex_counts_50s2) = read_stacked_bulk_molecules(
				cell_paths2, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts2 = np.minimum(
				complex_counts_30s2, complex_counts_50s2)

			# Active
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			active_ribosome_counts2 = read_stacked_columns(
				cell_paths2, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# Total
			total_ribosome_counts = active_ribosome_counts + inactive_ribosome_counts
			total_ribosome_counts2 = active_ribosome_counts2 + inactive_ribosome_counts2
			plt.plot(time / 60., total_ribosome_counts, color=LINE_COLOR)
			plt.plot(time2 / 60., total_ribosome_counts2, color=LINE_COLOR2)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,40000))
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Total Ribosome Counts", fontsize="small")

			plot_num += 1

			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(
				plt, plotOutDir, plotOutFileName + plot_suffix + "_"
				+ "VAR_" + str(VARIANT_1) + "_SEED_" + str(VARIANT_1_SEED) + "_"
				+ "VAR_" + str(VARIANT_2) + "_SEED_" + str(VARIANT_2_SEED) + "_GENS_"
				+ str(START_GEN) + "_" + str(END_GEN), metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
