"""
Save the average, initial, and final values for cell level properties for each
generation as csv files.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import inf
import pandas as pd

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
			# mpl.rcParams['axes.spines.right'] = False
			# mpl.rcParams['axes.spines.top'] = False
			# plt.figure(figsize = (12, total_plots*3))
			plot_num = 1
			# ax1 = plt.subplot(total_plots, 1, plot_num)

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

			all_averages = np.zeros((total_plots, len(gen_labels))) - 1.0
			all_averages_percent = np.zeros((total_plots, len(gen_labels))) - 1.0
			all_averages_by_gen = np.zeros((total_plots, len(unique_gen_labels))) - 1.0
			all_averages_index_mapping = {}
			all_initial_values_by_gen = np.zeros((total_plots, len(unique_gen_labels))) - 1.0
			all_final_values_by_gen = np.zeros((total_plots, len(unique_gen_labels))) - 1.0

			# Doubling time
			avg_before_gfp = np.mean(dt[4:8])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = dt[i]
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["doubling_time"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("Percent Change Doubling Time (min)", fontsize="x-small")
			# plt.title("Doubling Time")
			plot_num += 1

			# Mass
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			avg_before_gfp = np.mean(mass[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(mass[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["mass"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = mass[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = mass[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Percent Change in Mass (fg)", fontsize="x-small")
			# plt.title("Cell Mass")
			plot_num += 1

			# Total Ribosome Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
				cell_paths, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts = np.minimum(
				complex_counts_30s, complex_counts_50s)
			# Active
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# Total
			total_ribosome_counts = active_ribosome_counts + inactive_ribosome_counts
			avg_before_gfp = np.mean(total_ribosome_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(total_ribosome_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["total_ribosome_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = total_ribosome_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = total_ribosome_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Total Ribosome Counts", fontsize="x-small")
			# plt.title("Total Ribosome Counts")
			plot_num += 1

			#  Ribosome Components Data Loading
			# Load IDs of ribosome components from sim_data
			s30_protein_ids = sim_data.molecule_groups.s30_proteins
			s30_16s_rRNA_ids = sim_data.molecule_groups.s30_16s_rRNA
			s30_full_complex_id = [sim_data.molecule_ids.s30_full_complex]
			s50_protein_ids = sim_data.molecule_groups.s50_proteins
			s50_23s_rRNA_ids = sim_data.molecule_groups.s50_23s_rRNA
			s50_5s_rRNA_ids = sim_data.molecule_groups.s50_5s_rRNA
			s50_full_complex_id = [sim_data.molecule_ids.s50_full_complex]
			# Get complexation stoichiometries of ribosomal proteins
			complexation = sim_data.process.complexation
			s30_monomers = complexation.get_monomers(s30_full_complex_id[0])
			s50_monomers = complexation.get_monomers(s50_full_complex_id[0])
			s30_subunit_id_to_stoich = {
				subunit_id: stoich for (subunit_id, stoich)
				in zip(s30_monomers['subunitIds'], s30_monomers['subunitStoich'])
				}
			s50_subunit_id_to_stoich = {
				subunit_id: stoich for (subunit_id, stoich)
				in zip(s50_monomers['subunitIds'], s50_monomers['subunitStoich'])
				}
			s30_protein_stoich = np.array([
				s30_subunit_id_to_stoich[subunit_id]
				for subunit_id in s30_protein_ids
				])
			s50_protein_stoich = np.array([
				s50_subunit_id_to_stoich[subunit_id]
				for subunit_id in s50_protein_ids
				])

			# Read free counts of rRNAs
			time = read_stacked_columns(cell_paths, 'Main', 'time')
			(s30_16s_rRNA_counts, s50_23s_rRNA_counts, s50_5s_rRNA_counts
			) = read_stacked_bulk_molecules(
				cell_paths,
				(s30_16s_rRNA_ids, s50_23s_rRNA_ids, s50_5s_rRNA_ids))

			# Read counts of free 30S and 50S subunits
			(s30_full_complex_counts, s50_full_complex_counts
			) = read_stacked_bulk_molecules(
				cell_paths, (s30_full_complex_id, s50_full_complex_id))

			# Read counts of active ribosomes
			simOutDir = os.path.join(cell_paths[0], 'simOut')
			unique_molecule_counts_reader = TableReader(
				os.path.join(simOutDir, 'UniqueMoleculeCounts'))
			active_ribosome_index = unique_molecule_counts_reader.readAttribute(
				'uniqueMoleculeIds').index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				)[:, active_ribosome_index]

			# Read counts of free and complexed ribosomal proteins (this includes
			# ribosomal proteins in inactive ribosomal subunits and active
			# ribosomes)
			simOutDir = os.path.join(cell_paths[0], 'simOut')
			monomer_counts_reader = TableReader(
				os.path.join(simOutDir, 'MonomerCounts'))
			monomer_ids = monomer_counts_reader.readAttribute('monomerIds')
			monomer_id_to_index = {
				monomer_id: i for (i, monomer_id) in enumerate(monomer_ids)
				}
			s30_protein_indexes = np.array([
				monomer_id_to_index[protein_id] for protein_id in s30_protein_ids
				])
			s50_protein_indexes = np.array([
				monomer_id_to_index[protein_id] for protein_id in s50_protein_ids
				])
			monomer_counts = read_stacked_columns(
				cell_paths, 'MonomerCounts', 'monomerCounts',
				)
			# Divide by complexation stoichiometric constant
			s30_protein_counts = monomer_counts[:, s30_protein_indexes] / s30_protein_stoich
			s50_protein_counts = monomer_counts[:, s50_protein_indexes] / s50_protein_stoich

			# Calculate total counts of all components
			s30_limiting_protein_counts = s30_protein_counts.min(axis=1)
			s30_16s_rRNA_total_counts = (
					s30_16s_rRNA_counts.sum(axis=1)
					+ s30_full_complex_counts + active_ribosome_counts)
			s50_limiting_protein_counts = s50_protein_counts.min(axis=1)
			s50_23s_rRNA_total_counts = (
					s50_23s_rRNA_counts.sum(axis=1)
					+ s50_full_complex_counts + active_ribosome_counts)
			s50_5s_rRNA_total_counts = (
					s50_5s_rRNA_counts.sum(axis=1)
					+ s50_full_complex_counts + active_ribosome_counts)

			# Calculate total counts of both subunits
			s30_total_counts = s30_full_complex_counts + active_ribosome_counts
			s50_total_counts = s50_full_complex_counts + active_ribosome_counts


			# Plot 50s Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s50_total_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s50_total_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["50s_total_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s50_total_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s50_total_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("50S Total Counts", fontsize="x-small")
			# plt.title("50S Total Counts")
			plot_num += 1

			# Plot 50s 23s rRNA Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s50_23s_rRNA_total_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s50_23s_rRNA_total_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["50s_23s_rRNA_total_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s50_23s_rRNA_total_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s50_23s_rRNA_total_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("50S 23S rRNA Counts", fontsize="x-small")
			# plt.title("50S 23S rRNA Counts")
			plot_num += 1

			# Plot 50s 5s rRNA Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s50_5s_rRNA_total_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s50_5s_rRNA_total_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["50s_5s_rRNA_total_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s50_5s_rRNA_total_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s50_5s_rRNA_total_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("50S 5S rRNA Counts", fontsize="x-small")
			# plt.title("50S 5S rRNA Counts")
			plot_num += 1

			# Plot 50s Limiting Protein Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s50_limiting_protein_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s50_limiting_protein_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["50s_limiting_protein_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s50_limiting_protein_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s50_limiting_protein_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("50S Limiting Protein Counts", fontsize="x-small")
			# plt.title("50S Limiting Protein Counts")
			plot_num += 1

			# Plot 30s Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s30_total_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels)) - 1.0
			avg_per_gen_percent = np.zeros(len(unique_gen_labels)) - 1.0
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s30_total_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["30s_total_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s30_total_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s30_total_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("30S Total Counts", fontsize="x-small")
			# plt.title("30S Total Counts")
			plot_num += 1

			# Plot 30s 16s rRNA Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s30_16s_rRNA_total_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s30_16s_rRNA_total_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["30s_16s_rRNA_total_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s30_16s_rRNA_total_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s30_16s_rRNA_total_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Time (min)")
			# plt.ylabel("30S 16S rRNA Counts", fontsize="x-small")
			# plt.title("30S 16S rRNA Counts")
			plot_num += 1

			# Plot 30s Limiting Protein Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(s30_limiting_protein_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(s30_limiting_protein_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["30s_limiting_protein_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = s30_limiting_protein_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = s30_limiting_protein_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("30S Limiting Protein Counts", fontsize="x-small")
			# plt.title("30S Limiting Protein Counts")
			plot_num += 1

			# rRNA copy number
			GENE_ID_TO_RRNA_OPERON_ID = {
				'EG30084': 'rrnA',
				'EG30085': 'rrnB',
				'EG30086': 'rrnC',
				'EG30087': 'rrnD',
				'EG30088': 'rrnE',
				'EG30089': 'rrnG',
				'EG30090': 'rrnH',
				}
			# Get gene_ids attribute from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			gene_ids = rna_synth_prob_reader.readAttribute('gene_ids')
			# Get indexes of 16S genes (first gene in each operon)
			rrna_gene_indexes = np.array([
				gene_ids.index(key) for key in GENE_ID_TO_RRNA_OPERON_ID.keys()
				])
			# Get copy numbers of 16S genes
			rrna_gene_copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'gene_copy_number',
				ignore_exception=True, fun=lambda x: x[:, rrna_gene_indexes])
			# Make a separate plot for each rrna gene
			for i, gene_id in enumerate(GENE_ID_TO_RRNA_OPERON_ID.keys()):
				# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				avg_before_gfp = np.mean(
					rrna_gene_copy_numbers[gen_start_index[4]:gen_start_index[8], i])
				avg_per_gen = np.zeros(len(unique_gen_labels))
				avg_per_gen_percent = np.zeros(len(unique_gen_labels))
				for j in range(len(unique_gen_labels)):
					avg_per_gen[j] = np.mean(
						rrna_gene_copy_numbers[gen_start_index[j]:gen_end_index[j], i])
					avg_per_gen_percent[j] = (
						avg_per_gen[j] - avg_before_gfp) / avg_before_gfp * 100.0
				avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
				all_averages[plot_num, :] = avg_per_gen_plot
				avg_per_gen_percent_plot = np.repeat(
					avg_per_gen_percent, num_time_steps)
				all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
				all_averages_index_mapping[
					GENE_ID_TO_RRNA_OPERON_ID[gene_id] + "_copy_num"] = plot_num
				all_averages_by_gen[plot_num, :] = avg_per_gen
				all_initial_values_by_gen[plot_num, :] = (
					rrna_gene_copy_numbers[gen_start_index, i].squeeze())
				all_final_values_by_gen[plot_num, :] = (
					rrna_gene_copy_numbers[gen_end_index, i].squeeze())
				x_data = unique_gen_labels
				y_data = avg_per_gen
				# plt.axvline(
				# 	(8), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.axvline(
				# 	(16), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.plot(x_data, y_data)
				# plt.xlabel("Generation")
				# plt.ylabel("Copy Number", fontsize="x-small")
				# plt.title("Copy Number of " + GENE_ID_TO_RRNA_OPERON_ID[gene_id])
				plot_num += 1

			# rRNA RNA Synth Prob
			TU_ID_TO_RRNA_OPERON_ID = {
				'TU0-1181[c]': 'rrnA',
				'TU0-1182[c]': 'rrnB',
				'TU0-1183[c]': 'rrnC',
				'TU0-1191[c]': 'rrnD',
				'TU0-1186[c]': 'rrnE',
				'TU0-1187[c]': 'rrnG',
				'TU0-1189[c]': 'rrnH',
				}
			# Get rna ID attributes from reference cell path
			reference_cell_path = cell_paths[0]
			sim_out_dir = os.path.join(reference_cell_path, 'simOut')
			rnap_data_reader = TableReader(
				os.path.join(sim_out_dir, 'RnapData'))
			rnap_data_rna_ids = rnap_data_reader.readAttribute('rnaIds')
			rna_synth_prob_reader = TableReader(
				os.path.join(sim_out_dir, 'RnaSynthProb'))
			rna_synth_prob_rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')
			# Get indexes of rRNA operons in both listeners
			rnap_data_rrna_indexes = np.array([
				rnap_data_rna_ids.index(key)
				for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])
			rna_synth_prob_rrna_indexes = np.array([
				rna_synth_prob_rna_ids.index(key)
				for key in TU_ID_TO_RRNA_OPERON_ID.keys()
				])
			rrna_rna_target_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
				ignore_exception=True)[:, rna_synth_prob_rrna_indexes]
			rrna_actual_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
				ignore_exception=True)[:, rna_synth_prob_rrna_indexes]
			# Target prob
			for i, gene_id in enumerate(GENE_ID_TO_RRNA_OPERON_ID.keys()):
				# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				avg_before_gfp = np.mean(
					rrna_rna_target_synth_prob[gen_start_index[4]:gen_start_index[8], i])
				avg_per_gen = np.zeros(len(unique_gen_labels))
				avg_per_gen_percent = np.zeros(len(unique_gen_labels))
				for j in range(len(unique_gen_labels)):
					avg_per_gen[j] = np.mean(
						rrna_rna_target_synth_prob[gen_start_index[j]:gen_end_index[j], i])
					avg_per_gen_percent[j] = (
						avg_per_gen[j] - avg_before_gfp) / avg_before_gfp * 100.0
				avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
				all_averages[plot_num, :] = avg_per_gen_plot
				avg_per_gen_percent_plot = np.repeat(
					avg_per_gen_percent, num_time_steps)
				all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
				all_averages_index_mapping[
					GENE_ID_TO_RRNA_OPERON_ID[gene_id] + "_target_synth_prob"] = plot_num
				all_averages_by_gen[plot_num, :] = avg_per_gen
				all_initial_values_by_gen[plot_num, :] = (
					rrna_rna_target_synth_prob[gen_start_index + 1, i].squeeze()) # RNASynthProb is always 0 at first time step
				all_final_values_by_gen[plot_num, :] = (
					rrna_rna_target_synth_prob[gen_end_index, i].squeeze())
				x_data = unique_gen_labels
				y_data = avg_per_gen
				# plt.axvline(
				# 	(8), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.axvline(
				# 	(16), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.plot(x_data, y_data)
				# plt.xlabel("Generation")
				# plt.ylabel("Target RNA Synthesis Probability", fontsize="x-small")
				# plt.title("Target RNA Synthesis Probability of " + GENE_ID_TO_RRNA_OPERON_ID[gene_id])
				plot_num += 1
			# Actual prob
			for i, gene_id in enumerate(GENE_ID_TO_RRNA_OPERON_ID.keys()):
				# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				avg_before_gfp = np.mean(
					rrna_actual_synth_prob[gen_start_index[4]:gen_start_index[8], i])
				avg_per_gen = np.zeros(len(unique_gen_labels))
				avg_per_gen_percent = np.zeros(len(unique_gen_labels))
				for j in range(len(unique_gen_labels)):
					avg_per_gen[j] = np.mean(
						rrna_actual_synth_prob[gen_start_index[j]:gen_end_index[j], i])
					avg_per_gen_percent[j] = (
						avg_per_gen[j] - avg_before_gfp) / avg_before_gfp * 100.0
				avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
				all_averages[plot_num, :] = avg_per_gen_plot
				avg_per_gen_percent_plot = np.repeat(
					avg_per_gen_percent, num_time_steps)
				all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
				all_averages_index_mapping[
					GENE_ID_TO_RRNA_OPERON_ID[gene_id] + "_actual_synth_prob"] = plot_num
				all_averages_by_gen[plot_num, :] = avg_per_gen
				all_initial_values_by_gen[plot_num, :] = (
					rrna_actual_synth_prob[gen_start_index + 1, i].squeeze()) # RNASynthProb is always 0 at first time step
				all_final_values_by_gen[plot_num, :] = (
					rrna_actual_synth_prob[gen_end_index, i].squeeze())
				x_data = unique_gen_labels
				y_data = avg_per_gen
				# plt.axvline(
				# 	(8), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.axvline(
				# 	(16), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.plot(x_data, y_data)
				# plt.xlabel("Generation")
				# plt.ylabel("Actual RNA Synthesis Probability", fontsize="x-small")
				# plt.title("Actual RNA Synthesis Probability of " + GENE_ID_TO_RRNA_OPERON_ID[gene_id])
				plot_num += 1

			# Total RNAP Counts
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts,) = read_stacked_bulk_molecules(
				cell_paths, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			# Total
			total_rnap_counts = inactive_rnap_counts + active_rnap_counts
			avg_before_gfp = np.mean(total_rnap_counts[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(total_rnap_counts[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping["total_rnap_counts"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = total_rnap_counts[gen_start_index].squeeze()
			all_final_values_by_gen[plot_num, :] = total_rnap_counts[gen_end_index].squeeze()
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			plot_num += 1

			# Plot protein counts for RNAP subunits
			RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			rnap_subunit_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, RNAP_subunit_monomer_ids, "monomer")
			rnap_subunit_protein_counts = read_stacked_columns(
				cell_paths, "MonomerCounts", "monomerCounts",
				ignore_exception=True)[:, rnap_subunit_monomer_indexes]
			for r in range(len(rnap_subunit_monomer_indexes)):
				# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				avg_before_gfp = np.mean(
					rnap_subunit_protein_counts[gen_start_index[4]:gen_start_index[8], r])
				avg_per_gen = np.zeros(len(unique_gen_labels))
				avg_per_gen_percent = np.zeros(len(unique_gen_labels))
				for i in range(len(unique_gen_labels)):
					avg_per_gen[i] = np.mean(
						rnap_subunit_protein_counts[gen_start_index[i]:gen_end_index[i], r])
					avg_per_gen_percent[i] = (
						avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
				avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
				all_averages[plot_num, :] = avg_per_gen_plot
				avg_per_gen_percent_plot = np.repeat(
					avg_per_gen_percent, num_time_steps)
				all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
				all_averages_index_mapping[
					RNAP_subunit_monomer_ids[r] + "_protein_counts"] = plot_num
				all_averages_by_gen[plot_num, :] = avg_per_gen
				all_initial_values_by_gen[plot_num, :] = (
					rnap_subunit_protein_counts[gen_start_index, r].squeeze())
				all_final_values_by_gen[plot_num, :] = (
					rnap_subunit_protein_counts[gen_end_index, r].squeeze())
				x_data = unique_gen_labels
				y_data = avg_per_gen
				# plt.axvline(
				# 	(8), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.axvline(
				# 	(16), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.plot(x_data, y_data)
				# plt.xlabel("Generation")
				# plt.ylabel("Protein Counts: " + RNAP_subunit_monomer_ids[r], fontsize="x-small")
				# plt.title("RNAP Subunit Protein Counts")
				plot_num += 1

			# RNAP subunit cistron counts
			monomer_ids = sim_data.process.translation.monomer_data['id']
			cistron_ids = sim_data.process.translation.monomer_data[
				'cistron_id']
			monomer_to_cistron_id_dict = {
				monomer_id: cistron_ids[i] for i, monomer_id in
				enumerate(monomer_ids)}
			target_cistron_ids = [
				monomer_to_cistron_id_dict.get(RNAP_monomer_id) for
				RNAP_monomer_id in RNAP_subunit_monomer_ids]
			# get cistron indexes from cistron ids
			sim_dir = cell_paths[0]
			simOutDir = os.path.join(sim_dir, 'simOut')
			mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
			mRNA_cistron_idx_dict = {
				rna: i for i, rna in
				enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
			rnap_subunit_cistron_indexes = [
				mRNA_cistron_idx_dict[target_cistron_id] for target_cistron_id in
				target_cistron_ids]
			rnap_subunit_mRNA_cistron_counts = read_stacked_columns(
				cell_paths, "RNACounts", "mRNA_cistron_counts",
				ignore_exception=True)[:, rnap_subunit_cistron_indexes]
			for r in range(len(rnap_subunit_cistron_indexes)):
				# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				avg_before_gfp = np.mean(
					rnap_subunit_mRNA_cistron_counts[gen_start_index[4]:gen_start_index[8], r])
				avg_per_gen = np.zeros(len(unique_gen_labels))
				avg_per_gen_percent = np.zeros(len(unique_gen_labels))
				for i in range(len(unique_gen_labels)):
					avg_per_gen[i] = np.mean(
						rnap_subunit_mRNA_cistron_counts[gen_start_index[i]:gen_end_index[i], r])
					avg_per_gen_percent[i] = (
						avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
				avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
				all_averages[plot_num, :] = avg_per_gen_plot
				avg_per_gen_percent_plot = np.repeat(
					avg_per_gen_percent, num_time_steps)
				all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
				all_averages_index_mapping[
					target_cistron_ids[r] + "_mRNA_cistron_counts"] = plot_num
				all_averages_by_gen[plot_num, :] = avg_per_gen
				all_initial_values_by_gen[plot_num, :] = (
					rnap_subunit_mRNA_cistron_counts[gen_start_index, r].squeeze())
				all_final_values_by_gen[plot_num, :] = (
					rnap_subunit_mRNA_cistron_counts[gen_end_index, r].squeeze())
				x_data = unique_gen_labels
				y_data = avg_per_gen
				# plt.axvline(
				# 	(8), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.axvline(
				# 	(16), color='gray',
				# 	linestyle='--', lw=0.5)
				# plt.plot(x_data, y_data)
				# plt.xlabel("Generation")
				# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# 	plt.ylim((-1,4.5))
				# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				# 	plt.xlim(standard_xlim)
				# plt.xlabel("Time (min)")
				# plt.ylabel("mRNA Counts: " + target_cistron_ids[r], fontsize="x-small")
				# plt.title("RNAP Subunit mRNA Counts")
				plot_num += 1


			# Active RNAP Portion Allocation
			# Active RNAP Counts
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			# New Gene RNAP Portion
			new_gene_mRNA_indexes = new_gene_mRNA_indexes
			new_gene_rnap_counts = read_stacked_columns(
				cell_paths, "RNACounts", "partial_mRNA_counts",
				ignore_exception=True)[:, new_gene_mRNA_indexes].flatten()
			new_gene_rnap_portion = new_gene_rnap_counts / active_rnap_counts
			# rRNA RNAP Portion
			rrna_rnap_counts = np.sum(read_stacked_columns(
				cell_paths, "RNACounts", "partial_rRNA_counts",
				ignore_exception=True), axis = 1).flatten()
			rrna_rnap_portion = rrna_rnap_counts / active_rnap_counts
			# RNAP Subunit RNAP Portion
			RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			rnap_subunit_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, RNAP_subunit_monomer_ids, "mRNA")
			rnap_subunit_rnap_counts = np.sum(read_stacked_columns(
				cell_paths, "RNACounts", "partial_mRNA_counts",
				ignore_exception=True)[:, rnap_subunit_mRNA_indexes], axis=1).flatten()
			rnap_subunit_rnap_portion = rnap_subunit_rnap_counts / active_rnap_counts
			# Ribosomal Proteins RNAP Portion
			ribosomal_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
			ribosomal_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, ribosomal_monomer_ids, "mRNA")
			ribosomal_rnap_counts = np.sum(read_stacked_columns(
				cell_paths, "RNACounts", "partial_mRNA_counts",
				ignore_exception=True)[:, ribosomal_mRNA_indexes], axis=1).flatten()
			ribosomal_rnap_portion = ribosomal_rnap_counts / active_rnap_counts
			# Plot
			# RNAP subunit RNAP portion
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(
				rnap_subunit_rnap_portion[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(
					rnap_subunit_rnap_portion[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (
					avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(
				avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping[
				"RNAP_subunit_rnap_portion"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = (
				rnap_subunit_rnap_portion[gen_start_index].squeeze())
			all_final_values_by_gen[plot_num, :] = (
				rnap_subunit_rnap_portion[gen_end_index].squeeze())
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("RNAP Subunit Portion of Active RNAPs", fontsize='x-small')
			# plt.title("Allocation of Active RNAPs")
			# plt.legend(fontsize="x-small")
			plot_num += 1

			# Ribosomal RNAP Portion
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(
				ribosomal_rnap_portion[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(
					ribosomal_rnap_portion[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (
					avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(
				avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping[
				"ribosomal_rnap_portion"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = (
				ribosomal_rnap_portion[gen_start_index].squeeze())
			all_final_values_by_gen[plot_num, :] = (
				ribosomal_rnap_portion[gen_end_index].squeeze())
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("Ribo Prot Portion of Active RNAPs", fontsize='x-small')
			# plt.title("Allocation of Active RNAPs")
			# plt.legend(fontsize="x-small")
			plot_num += 1

			# rRNA RNAP Portion
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(
				rrna_rnap_portion[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(
					rrna_rnap_portion[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (
					avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(
				avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping[
				"rrna_rnap_portion"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = (
				rrna_rnap_portion[gen_start_index].squeeze())
			all_final_values_by_gen[plot_num, :] = (
				rrna_rnap_portion[gen_end_index].squeeze())
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("rRNA Portion of Active RNAPs", fontsize='x-small')
			# plt.title("Allocation of Active RNAPs")
			# plt.legend(fontsize="x-small")
			plot_num += 1

			# Active Ribosome Portion Allocation
			# Active Ribosome Counts
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# New Gene Ribosome Portion
			new_gene_ribosome_counts = read_stacked_columns(
				cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
				ignore_exception=True)[:, new_gene_monomer_indexes].flatten()
			new_gene_ribosome_portion = new_gene_ribosome_counts / active_ribosome_counts
			# RNAP Subunits Ribosome Portion
			RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			rnap_subunit_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, RNAP_subunit_monomer_ids, "monomer")
			rnap_subunit_ribosome_counts = np.sum(read_stacked_columns(
				cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
				ignore_exception=True)[:, rnap_subunit_monomer_indexes], axis=1).flatten()
			rnap_subunit_ribosome_portion = rnap_subunit_ribosome_counts / active_ribosome_counts
			# Ribosomal Proteins Ribosome Portion
			ribosomal_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
			ribosomal_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, ribosomal_monomer_ids, "monomer")
			ribosomal_ribosome_counts = np.sum(read_stacked_columns(
				cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
				ignore_exception=True)[:, ribosomal_monomer_indexes], axis=1).flatten()
			ribosomal_ribosome_portion = ribosomal_ribosome_counts / active_ribosome_counts

			# Plot RNAP Subunit active ribosome portion
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(
				rnap_subunit_ribosome_portion[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(
					rnap_subunit_ribosome_portion[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(
				avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping[
				"rnap_subunit_ribosome_portion"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = (
				rnap_subunit_ribosome_portion[gen_start_index].squeeze())
			all_final_values_by_gen[plot_num, :] = (
				rnap_subunit_ribosome_portion[gen_end_index].squeeze())
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("RNAP Subunit Portion of Active Ribosomes", fontsize='x-small')
			# plt.title("Allocation of Active Ribosomes")
			plot_num += 1

			# Plot ribosomal protein active ribosome portion
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			avg_before_gfp = np.mean(
				ribosomal_ribosome_portion[gen_start_index[4]:gen_start_index[8]])
			avg_per_gen = np.zeros(len(unique_gen_labels))
			avg_per_gen_percent = np.zeros(len(unique_gen_labels))
			for i in range(len(unique_gen_labels)):
				avg_per_gen[i] = np.mean(
					ribosomal_ribosome_portion[gen_start_index[i]:gen_end_index[i]])
				avg_per_gen_percent[i] = (avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
			avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
			all_averages[plot_num, :] = avg_per_gen_plot
			avg_per_gen_percent_plot = np.repeat(
				avg_per_gen_percent, num_time_steps)
			all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
			all_averages_index_mapping[
				"ribosomal_ribosome_portion"] = plot_num
			all_averages_by_gen[plot_num, :] = avg_per_gen
			all_initial_values_by_gen[plot_num, :] = (
				ribosomal_ribosome_portion[gen_start_index].squeeze())
			all_final_values_by_gen[plot_num, :] = (
				ribosomal_ribosome_portion[gen_end_index].squeeze())
			x_data = unique_gen_labels
			y_data = avg_per_gen
			# plt.axvline(
			# 	(8), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.axvline(
			# 	(16), color='gray',
			# 	linestyle='--', lw=0.5)
			# plt.plot(x_data, y_data)
			# plt.xlabel("Generation")
			# plt.ylabel("Ribo Prot Portion of Active Ribosomes", fontsize='x-small')
			# plt.title("Allocation of Active Ribosomes")
			plot_num += 1

			# Cell mass breakdown
			mass_columns = [
				"dryMass", "dnaMass", "mRnaMass", "proteinMass", "rRnaMass",
				"tRnaMass", "membrane_mass", "waterMass", "smallMoleculeMass"]
			mass_labels = [
				"mass_dry_fg", "mass_dna_fg", "mass_mrna_fg", "mass_protein_fg",
				"mass_rrna_fg", "mass_trna_fg", "mass_membrane_fg",
				"mass_water_fg", "mass_small_molecule_fg"]
			for k, mass_column in enumerate(mass_columns):
				mass_label = mass_labels[k]
				mass_data = read_stacked_columns(
					cell_paths, 'Mass', mass_column)
				avg_before_gfp = np.mean(
					mass_data[gen_start_index[4]:gen_start_index[8]])
				avg_per_gen = np.zeros(len(unique_gen_labels))
				avg_per_gen_percent = np.zeros(len(unique_gen_labels))
				for i in range(len(unique_gen_labels)):
					avg_per_gen[i] = np.mean(
						mass_data[gen_start_index[i]:gen_end_index[i]])
					avg_per_gen_percent[i] = (
						avg_per_gen[i] - avg_before_gfp) / avg_before_gfp * 100.0
				avg_per_gen_plot = np.repeat(avg_per_gen, num_time_steps)
				all_averages[plot_num, :] = avg_per_gen_plot
				avg_per_gen_percent_plot = np.repeat(
					avg_per_gen_percent, num_time_steps)
				all_averages_percent[plot_num, :] = avg_per_gen_percent_plot
				all_averages_index_mapping[mass_label] = plot_num
				all_averages_by_gen[plot_num, :] = avg_per_gen
				all_initial_values_by_gen[plot_num, :] = (
					mass_data[gen_start_index].squeeze())
				all_final_values_by_gen[plot_num, :] = (
					mass_data[gen_end_index].squeeze())
				x_data = unique_gen_labels
				y_data = avg_per_gen
				plot_num += 1

			# TODO: RNAP subunit synth prob

			# TODO: RNAP subunit copy numbers

			print("Total number of plots made: ", plot_num)
			# plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			# exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			# plt.close("all")

			# Save csv files of average, initial, and final values for each generation and data keyword
			sim_seed = self.ap.get_cell_seed(cell_paths[0])
			import pandas as pd
			column_names = ["seed", "data_keyword"]
			for gen_num in unique_gen_labels:
				column_names.append(gen_num)

			all_averages_df = pd.DataFrame(columns=column_names)
			for data_keyword in all_averages_index_mapping.keys():
				data_index = all_averages_index_mapping[data_keyword]
				data_averages = all_averages_by_gen[data_index, :]
				new_row = [sim_seed, data_keyword] + list(data_averages)
				all_averages_df.loc[len(all_averages_df)] = new_row
			all_averages_df.to_csv(
				os.path.join(plotOutDir, f"all_averages_by_gen_for_seed_{sim_seed}.csv"),
				index=False)

			all_initial_values_df = pd.DataFrame(columns=column_names)
			for data_keyword in all_averages_index_mapping.keys():
				if data_keyword == "doubling_time":
					continue
				data_index = all_averages_index_mapping[data_keyword]
				data_initial_values = all_initial_values_by_gen[data_index, :]
				new_row = [sim_seed, data_keyword] + list(data_initial_values)
				all_initial_values_df.loc[len(all_initial_values_df)] = new_row
			all_initial_values_df.to_csv(
				os.path.join(plotOutDir, f"all_initial_values_by_gen_for_seed_{sim_seed}.csv"),
				index=False)

			all_final_values_df = pd.DataFrame(columns=column_names)
			for data_keyword in all_averages_index_mapping.keys():
				if data_keyword == "doubling_time":
					continue
				data_index = all_averages_index_mapping[data_keyword]
				data_final_values = all_final_values_by_gen[data_index, :]
				new_row = [sim_seed, data_keyword] + list(data_final_values)
				all_final_values_df.loc[len(all_final_values_df)] = new_row
			all_final_values_df.to_csv(
				os.path.join(plotOutDir, f"all_final_values_by_gen_for_seed_{sim_seed}.csv"),
				index=False)

			# print(all_averages_index_mapping.keys())
			#
			# comparisons = {
			# 	"production_machinery": [
			# 		"total_rnap_counts", "total_ribosome_counts",
			# 		],
			# 	"ribosome subunits": [
			# 		"50s_total_counts", "30s_total_counts",
			# 		],
			# 	"30s_subunit": [
			# 		'30s_16s_rRNA_total_counts', '30s_limiting_protein_counts'
			# 		],
			# 	"50s_subunit": [
			# 		'50s_23s_rRNA_total_counts', '50s_5s_rRNA_total_counts',
			# 		'50s_limiting_protein_counts'
			# 		],
			# 	"rnap_portions": [
			# 		'RNAP_subunit_rnap_portion', 'ribosomal_rnap_portion', 'rrna_rnap_portion',
			# 		],
			# 	"ribosome_portions": [
			# 		'rnap_subunit_ribosome_portion', 'ribosomal_ribosome_portion',
			# 		],
			# 	"rnap_subunit_protein_counts": [
			# 		'EG10893-MONOMER[c]_protein_counts', 'RPOC-MONOMER[c]_protein_counts',
			# 		'RPOB-MONOMER[c]_protein_counts',
			# 		],
			# 	"rnap_subunit_mRNA_counts": [
			# 		'EG10893_RNA_mRNA_cistron_counts', 'EG10895_RNA_mRNA_cistron_counts',
			# 		'EG10894_RNA_mRNA_cistron_counts',
			# 		],
			# 	"rnap_subunit_protein_and_mRNA_counts": [
			# 		'EG10893-MONOMER[c]_protein_counts', 'RPOC-MONOMER[c]_protein_counts',
			# 		'RPOB-MONOMER[c]_protein_counts', 'EG10893_RNA_mRNA_cistron_counts',
			# 		'EG10895_RNA_mRNA_cistron_counts', 'EG10894_RNA_mRNA_cistron_counts',
			# 		],
			# 	"rRNA_gene_copy_num": [
			# 		'rrnA_copy_num', 'rrnB_copy_num', 'rrnC_copy_num', 'rrnD_copy_num',
			# 		'rrnE_copy_num', 'rrnG_copy_num', 'rrnH_copy_num',
			# 		],
			# 	'rRNA_gene_target_synth_prob': [
			# 		'rrnA_target_synth_prob', 'rrnB_target_synth_prob', 'rrnC_target_synth_prob',
			# 		'rrnD_target_synth_prob', 'rrnE_target_synth_prob', 'rrnG_target_synth_prob',
			# 		'rrnH_target_synth_prob',
			# 		],
			# 	'rRNA_gene_actual_synth_prob': [
			# 		'rrnA_actual_synth_prob', 'rrnB_actual_synth_prob', 'rrnC_actual_synth_prob',
			# 		'rrnD_actual_synth_prob', 'rrnE_actual_synth_prob', 'rrnG_actual_synth_prob',
			# 		'rrnH_actual_synth_prob',
			# 		],
			# 	}
			# unique_plot_labels = set()
			# for labels in comparisons.values():
			# 	if not labels[0].startswith("rrn"):
			# 		unique_plot_labels.update(labels)
			# unique_plot_labels.update([
			# 	'rrnA_copy_num', 'rrnA_target_synth_prob',
			# 	'rrnA_actual_synth_prob'])
			# unique_plot_labels = list(unique_plot_labels)
			# comparisons["all"] = unique_plot_labels
			#
			# extended_colors = [
			# 	"#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
			# 	"#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100",
			# 	"#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73",
			# 	"#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000",
			# 	"#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
			# 	"#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02"
			# 	]
			#
			# for comparison, labels in comparisons.items():
			# 	plt.figure(figsize=(10, 7))
			# 	plt.axvline(
			# 		(8), color='gray',
			# 		linestyle='--', lw=0.5)
			# 	plt.axvline(
			# 		(16), color='gray',
			# 		linestyle='--', lw=0.5)
			#
			# 	for c, data_label in enumerate(labels):
			# 		plot_num_index = all_averages_index_mapping[data_label]
			# 		avg_per_gen = all_averages_by_gen[plot_num_index, :]
			#
			# 		# Make y data the normalized slope per gen
			# 		y_data = np.zeros(len(unique_gen_labels))
			# 		for j in range(len(unique_gen_labels)):
			# 			if j == 0:
			# 				y_data[j] = 0
			# 			else:
			# 				y_data[j] = (avg_per_gen[j]) / (
			# 					avg_per_gen[7])
			# 		plt.plot(
			# 			unique_gen_labels, y_data,
			# 			label=data_label, alpha=0.5, color=extended_colors[c])
			# 	plt.xlabel("Generation")
			# 	plt.ylabel("Avg of Gen / Avg of Gen 7")
			# 	plt.title(comparison + ": Normalized Averages",
			# 			  fontsize='x-small')
			# 	plt.legend(loc='upper right', fontsize="x-small")
			# 	exportFigure(plt, plotOutDir, plotOutFileName + "_" + comparison, metadata)
			# 	plt.close("all")




if __name__ == '__main__':
	Plot().cli()