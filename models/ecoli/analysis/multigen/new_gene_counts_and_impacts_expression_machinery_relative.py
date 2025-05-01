"""
Plot mRNA and protein counts for new genes across multiple generations, as well
as plots to analyze the impact of new gene expression, including growth rate,
RNAP and ribosome counts, and ppGpp concentration.
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
	read_stacked_bulk_molecules, read_stacked_columns, read_bulk_molecule_counts)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

LINE_COLOR = (66/255, 170/255, 154/255)

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

		cell_paths = self.ap.get_cells()
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

		# TODO: make a option for this plot that you could run without new genes

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
		(new_gene_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, new_gene_monomer_ids, ignore_exception=True)
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)
		new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]
		new_gene_promoter_copy_numbers = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True)[:,new_gene_RNA_indexes]

		# plot_suffixes = ["", "_standard_axes_y", "_standard_axes_both"]
		plot_suffixes = [""]
		standard_xlim = (0,2000)
		total_plots = 100 # TODO Modularize and get rid of this magic number

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.figure(figsize = (12, total_plots*3))
			plot_num = 1
			ax1 = plt.subplot(total_plots, 1, plot_num)

			# Get time marker where GFP induced
			dt = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			num_time_steps = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: len(x)).squeeze()
			gen_labels = np.repeat(np.arange(len(dt)), num_time_steps)
			unique_gen_labels = np.unique(gen_labels)
			gen_start_index = np.array(
				[gen_labels.tolist().index(i) for i in unique_gen_labels])
			gen_end_index = np.concatenate((
				np.array(gen_start_index[1:] - 1), np.array([len(gen_labels) - 1])))
			dt_to_plot = np.repeat(dt, num_time_steps)
			doubling_time_index = np.where(dt_to_plot == dt[8])[0][0]
			analysis_gen_index = np.where(dt_to_plot == dt[16])[0][0]

			# Doubling time
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., dt_to_plot)
			plt.xlabel("Time (min)")
			plt.ylabel("Doubling Time (min)", fontsize="small")
			plt.title("Doubling Time")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., mass, color=LINE_COLOR)
			plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)", fontsize="x-small")
			plt.title("Cell Mass")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			avg_before_gfp = np.mean(mass[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(mass[gen_start_index[i]:gen_end_index[i]])
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			# plot horizontal line for average mass before gfp
			plt.axhline(
				avg_before_gfp, color='gray', linestyle='--', lw=0.5)
			plt.plot(time / 60., avg_per_gen_plot, color=LINE_COLOR)
			plt.ylim(bottom=0)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)", fontsize="x-small")
			plt.title("Cell Mass")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			avg_before_gfp = np.mean(mass[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(mass[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., avg_per_gen_percent_plot, color=LINE_COLOR)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)", fontsize="x-small")
			plt.title("Cell Mass")
			plot_num += 1

			print("Total number of plots made: ", plot_num)
			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			plt.close("all")


if __name__ == '__main__':
	Plot().cli()
