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
			dt_to_plot = np.repeat(dt, num_time_steps)
			doubling_time_index = np.where(dt_to_plot == dt[8])[0][0]
			analysis_gen_index = np.where(dt_to_plot == dt[16])[0][0]

			# Growth Rate
			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(500, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time.flatten() / 60., growth_rate_convolved, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,.0004)
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Growth Rate", fontsize="x-small")
			plt.title("Growth Rate")
			plot_num += 1

			# Doubling time
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., dt_to_plot)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 80))
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
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,4000)
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)", fontsize="x-small")
			plt.title("Cell Mass")
			plot_num += 1

			# # Gene Promoter Copy Number
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# if len(new_gene_mRNA_ids) == 1:
			# 	plt.plot(time / 60., new_gene_promoter_copy_numbers,
			# 			 color=LINE_COLOR)
			# else:
			# 	for r in range(len(new_gene_mRNA_ids)):
			# 		plt.plot(time / 60., new_gene_promoter_copy_numbers[:,r],
			# 				 label = new_gene_mRNA_ids[r], color=LINE_COLOR)
			# 	plt.legend()
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	plt.ylim((0,6))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	plt.xlim(standard_xlim)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Gene Promoter Copy Number", fontsize="x-small")
			# plt.title("New Gene Promoter Copy Number")
			# plot_num += 1

			# mRNA Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			if plot_suffix == "":
				if len(new_gene_mRNA_ids) == 1:
					plt.plot(time / 60., new_gene_mRNA_counts, color=LINE_COLOR)
				else:
					for r in range(len(new_gene_mRNA_ids)):
						plt.plot(time / 60., new_gene_mRNA_counts[:,r],
								 label = new_gene_mRNA_ids[r],
								 color=LINE_COLOR)
					plt.legend()
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# plot on log scale instead
				if len(new_gene_mRNA_ids) == 1:
					plt.plot(time / 60., np.log10(new_gene_mRNA_counts + 1),
							 color=LINE_COLOR)
				else:
					for r in range(len(new_gene_mRNA_ids)):
						plt.plot(time / 60., np.log10(new_gene_mRNA_counts[:,r] + 1),
								 label = new_gene_mRNA_ids[r],
								 color=LINE_COLOR)
					plt.legend()
				plt.ylim((-1,4.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Log(mRNA Counts + 1)", fontsize="x-small")
			plt.title("New Gene mRNA Counts")
			plot_num += 1

			# Protein Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			if plot_suffix == "":
				if len(new_gene_monomer_ids) == 1:
					plt.plot(time / 60., new_gene_monomer_counts,
							 color=LINE_COLOR)
				else:
					for m in range(len(new_gene_monomer_ids)):
						plt.plot(time / 60., new_gene_monomer_counts[:,m],
								 label = new_gene_monomer_ids[m],
								 color=LINE_COLOR)
					plt.legend()
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# plot on log scale instead
				if len(new_gene_monomer_ids) == 1:
					plt.plot(time / 60., np.log10(new_gene_monomer_counts + 1),
							 color=LINE_COLOR)
				else:
					for m in range(len(new_gene_monomer_ids)):
						plt.plot(time / 60., np.log10(new_gene_monomer_counts[:,m] + 1),
								 label = new_gene_monomer_ids[m],
								 color=LINE_COLOR)
					plt.legend()
				plt.ylim((-1,7.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Log(Protein Counts + 1)", fontsize="x-small")
			plt.title("New Gene Protein Counts")
			plot_num += 1

			# Amino acid concentrations
			# aa_ids_of_interest = ['HIS[c]', 'THR[c]', 'TRP[c]']
			aa_ids_of_interest = sim_data.molecule_groups.amino_acids
			targets = np.array(
				[sim_data.process.metabolism.conc_dict[key].asNumber(units.mmol / units.L) for key
					in aa_ids_of_interest])
			aa_counts = read_stacked_bulk_molecules(
				cell_paths, aa_ids_of_interest, remove_first=True, ignore_exception=True)[0]

			counts_to_molar = read_stacked_columns(
				cell_paths, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True, ignore_exception=True)
			aa_conc = aa_counts * counts_to_molar

			for ii, aa in enumerate(aa_ids_of_interest):
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.axvline(
					(time[doubling_time_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.axvline(
					(time[analysis_gen_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.plot(time_no_first / 60., aa_conc[:, aa_ids_of_interest.index(aa)], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((0, 0.1))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.axhline(targets[ii], color='r', linestyle='--')
				plt.xlabel("Time (min)")
				plt.ylabel("Concentration (mmol/L)", fontsize='x-small')
				plt.title("Amino Acid Concentration: " + aa)
				plot_num += 1

			# ppGpp
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ppGpp_concentration = read_stacked_columns(
				cell_paths, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)
			plt.plot(time_no_first / 60., ppGpp_concentration, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,150))
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("ppGpp Concentration (uM)", fontsize="x-small")
			plt.title("ppGpp Concentration")
			plot_num += 1

			# Total RNAP Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
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
			plt.plot(time / 60., total_rnap_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,10000))
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Total RNAP Counts", fontsize="x-small")
			plt.title("Total RNAP Counts")
			plot_num += 1

			# Total Ribosome Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
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
			plt.plot(time / 60., total_ribosome_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,40000))
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Total Ribosome Counts", fontsize="x-small")
			plt.title("Total Ribosome Counts")
			plot_num += 1

			# # Glucose Comsumption Rate
			# # TODO: extend to other carbon sources
			# GLUCOSE_ID = "GLC[p]"
			# FLUX_UNITS = units.mmol / units.g / units.h
			# MASS_UNITS = units.fg
			# GROWTH_UNITS = units.fg / units.s
			# fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
			# external_molecule_ids = np.array(
			# 	fba_results.readAttribute("externalMoleculeIDs"))
			# fba_results.close()
			# if GLUCOSE_ID not in external_molecule_ids:
			# 	print("This plot only runs when glucose is the carbon source.")
			# 	return
			# glucose_idx = np.where(external_molecule_ids == GLUCOSE_ID)[0][0]
			# glucose_flux = FLUX_UNITS * read_stacked_columns(
			# 	cell_paths, "FBAResults", "externalExchangeFluxes",
			# 	ignore_exception=True, remove_first=True)[:, glucose_idx]
			# glucose_mw = sim_data.getter.get_mass(GLUCOSE_ID)
			# cell_dry_mass = MASS_UNITS * read_stacked_columns(
			# 	cell_paths, "Mass", "dryMass", ignore_exception=True,
			# 	remove_first=True).flatten()
			# glucose_mass_flux = glucose_flux * glucose_mw * cell_dry_mass
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# moving_window = min(300, len(glucose_mass_flux.asNumber()))
			# convolution_array = (np.ones(moving_window) / moving_window)
			# glucose_mass_flux_convolved = np.convolve(
			# 	convolution_array, glucose_mass_flux.asNumber(), mode='same')
			# plt.plot(time_no_first / 60., -glucose_mass_flux_convolved)
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	plt.ylim((0,1800))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	plt.xlim(standard_xlim)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Glucose Comsumption Rate (fg/h)", fontsize='x-small')
			# plt.title("Glucose Consumption Rate")
			# plot_num += 1

			# New Gene RNA Synthesis Probability
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			new_gene_target_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
				ignore_exception=True, remove_first=True)[:, new_gene_RNA_indexes]
			new_gene_actual_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
				ignore_exception=True, remove_first=True)[:, new_gene_RNA_indexes]

			if len(new_gene_mRNA_ids) == 1:
				plt.plot(time_no_first / 60., new_gene_target_rna_synth_prob,
						 label="Target")
				plt.plot(time_no_first / 60., new_gene_actual_rna_synth_prob,
						 label="Actual")
			else:
				for r in range(len(new_gene_mRNA_ids)):
					plt.plot(time_no_first / 60., new_gene_target_rna_synth_prob[:,r],
							 label = new_gene_mRNA_ids[r] + ": Target")
					plt.plot(time_no_first / 60., new_gene_actual_rna_synth_prob[:, r],
							 label=new_gene_mRNA_ids[r] + ": Actual")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1,0.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("RNA Synthesis Probability", fontsize='x-small')
			plt.title("New Gene RNA Synthesis Probability")
			plt.legend()
			plot_num += 1

			# New Gene Protein Initialization Probability
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			new_gene_target_protein_init_prob = read_stacked_columns(
				cell_paths, 'RibosomeData', 'target_prob_translation_per_transcript',
				ignore_exception=True, remove_first=True)[:, new_gene_monomer_indexes]
			new_gene_actual_protein_init_prob = read_stacked_columns(
				cell_paths, 'RibosomeData', 'actual_prob_translation_per_transcript',
				ignore_exception=True, remove_first=True)[:, new_gene_monomer_indexes]

			if len(new_gene_monomer_ids) == 1:
				plt.plot(time_no_first / 60., new_gene_target_protein_init_prob,
						 label="Target")
				plt.plot(time_no_first / 60., new_gene_actual_protein_init_prob,
						 label="Actual")
			else:
				for r in range(len(new_gene_monomer_ids)):
					plt.plot(time_no_first / 60., new_gene_target_protein_init_prob[:,r],
							 label = new_gene_monomer_ids[r] + ": Target")
					plt.plot(time_no_first/ 60., new_gene_actual_protein_init_prob[:, r],
							 label=new_gene_monomer_ids[r] + ": Actual")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 1.1))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Probability Translation Per Transcript", fontsize='x-small')
			plt.title("New Gene Protein Initialization Probability")
			plt.legend()
			plot_num += 1

			# mRNA Counts, Gene Promoter Copy Number, and RNA Synth Prob
			ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			ax2.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ax2.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ax3 = ax2.twinx()
			ax3.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ax3.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			# plot on log scale instead
			ax2.plot(time_no_first / 60., new_gene_target_rna_synth_prob,
					 label="Target")
			ax2.plot(time_no_first / 60., new_gene_actual_rna_synth_prob,
					 label="Actual")
			ax3.plot(time / 60., np.log10(new_gene_mRNA_counts + 1), label = "log10(mRNA counts + 1)", color = "cyan")
			ax3.plot(time / 60., 2 * new_gene_promoter_copy_numbers, label="Copy Number", color = "red")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				ax2.set_ylim(((-0.1, 0.5)))
				ax3.set_ylim((-1,4.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				ax2.set_xlim(standard_xlim)
			ax2.set_xlabel("Time (min)")
			plt.title("New Gene mRNA Counts, Copy Number, and RNA Synth Prob")
			plot_num += 1

			# mRNA Counts and RNAP Counts
			ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			ax2.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ax2.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ax3 = ax2.twinx()
			ax3.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			ax3.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			# plot on log scale instead
			ax2.plot(time / 60., total_rnap_counts)
			ax3.plot(time / 60., np.log10(new_gene_mRNA_counts + 1), label = "log10(mRNA counts + 1)", color = "cyan")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				ax2.set_ylim((0,10000))
				ax3.set_ylim((-1,4.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				ax2.set_xlim(standard_xlim)
			ax2.set_xlabel("Time (min)")
			plt.title("New Gene mRNA Counts and RNAP Counts")
			plot_num += 1

			# TODO: New Gene mRNA mass fraction

			# TODO: New Gene Protein Mass Fraction

			# Active RNAP Portion Allocation
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
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
			plt.plot(time / 60., rnap_subunit_rnap_portion, label="RNAP Subunit")
			plt.plot(time / 60., ribosomal_rnap_portion, label="Ribo. Prot.")
			plt.plot(time / 60., new_gene_rnap_portion, label = "New Gene")
			plt.plot(time / 60., rrna_rnap_portion, label="rRNA")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 1.1))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Portion of Active RNAPs", fontsize='x-small')
			plt.title("Allocation of Active RNAPs")
			plt.legend(fontsize="x-small")
			plot_num += 1

			# Active Ribosome Portion Allocation
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
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
			# Plot
			plt.plot(time / 60., rnap_subunit_ribosome_portion, label="RNAP Subunit")
			plt.plot(time / 60., ribosomal_ribosome_portion, label="Ribo. Prot.")
			plt.plot(time / 60., new_gene_ribosome_portion, label = "New Gene")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 1.1))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Portion of Active Ribosomes", fontsize='x-small')
			plt.title("Allocation of Active Ribosomes")
			plt.legend(fontsize="x-small")
			plot_num += 1

			# Get RNA Synth Prob of RNAP subunits
			rnap_subunit_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, RNAP_subunit_monomer_ids, "mRNA")
			rnap_subunit_mRNA_ids = list(self.get_mRNA_ids_from_monomer_ids(sim_data, RNAP_subunit_monomer_ids))
			rnap_subunit_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
				ignore_exception=True)[:, rnap_subunit_mRNA_indexes]
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			for ii in range(len(rnap_subunit_mRNA_indexes)):
				plt.plot(time / 60., rnap_subunit_rna_synth_prob[:,ii], label = rnap_subunit_mRNA_ids[ii])
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 0.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Target RNA Synthesis Probability", fontsize='x-small')
			plt.title("RNAP Subunit RNA Synthesis Probability")
			plt.legend(fontsize="x-small")
			plot_num += 1

			# Get Protein Init Prob of RNAP subunits
			rnap_subunit_protein_init_prob = read_stacked_columns(
				cell_paths, 'RibosomeData', 'target_prob_translation_per_transcript',
				ignore_exception=True)[:, rnap_subunit_monomer_indexes]
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			for ii in range(len(rnap_subunit_monomer_indexes)):
				plt.plot(time / 60., rnap_subunit_protein_init_prob[:,ii], label = RNAP_subunit_monomer_ids[ii])
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1, 1.1))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Target Probability Translation Per Transcript", fontsize='x-small')
			plt.title("RNAP Subunit Protein Initialization Probability")
			plt.legend(fontsize="x-small")
			plot_num += 1

			# Plot mRNA Counts for RNAP subunits

			# map rnap subunit monomer names to cistron indexes
			# Map protein ids to cistron ids
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
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.axvline(
					(time[doubling_time_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.axvline(
					(time[analysis_gen_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.plot(time / 60., rnap_subunit_mRNA_cistron_counts[:,r], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((-1,4.5))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("mRNA Counts: " + target_cistron_ids[r], fontsize="x-small")
				plt.title("RNAP Subunit mRNA Counts")
				plot_num += 1

			# Plot protein counts for RNAP subunits
			rnap_subunit_protein_counts = read_stacked_columns(
				cell_paths, "MonomerCounts", "monomerCounts",
				ignore_exception=True)[:, rnap_subunit_monomer_indexes]
			for r in range(len(rnap_subunit_monomer_indexes)):
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.axvline(
					(time[doubling_time_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.axvline(
					(time[analysis_gen_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.plot(time / 60., rnap_subunit_protein_counts[:,r], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((-1,7.5))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("Protein Counts: " + RNAP_subunit_monomer_ids[r], fontsize="x-small")
				plt.title("RNAP Subunit Protein Counts")
				plot_num += 1

			#  Ribosome Components
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

			# Plot 30s limiting protein counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s30_limiting_protein_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("30S Limiting Protein Counts", fontsize="x-small")
			plt.title("30S Limiting Protein Counts")
			plot_num += 1

			# Plot 30s 16s rRNA total counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s30_16s_rRNA_total_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("30S 16S rRNA Total Counts", fontsize="x-small")
			plt.title("30S 16S rRNA Total Counts")
			plot_num += 1

			# Plot 50s limiting protein counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s50_limiting_protein_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("50S Limiting Protein Counts", fontsize="x-small")
			plt.title("50S Limiting Protein Counts")
			plot_num += 1

			# Plot 50s 23s rRNA total counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s50_23s_rRNA_total_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("50S 23S rRNA Total Counts", fontsize="x-small")
			plt.title("50S 23S rRNA Total Counts")
			plot_num += 1

			# Plot 50s 5s rRNA total counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s50_5s_rRNA_total_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("50S 5S rRNA Total Counts", fontsize="x-small")
			plt.title("50S 5S rRNA Total Counts")
			plot_num += 1

			# Plot 30s total counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s30_total_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("30S Total Counts", fontsize="x-small")
			plt.title("30S Total Counts")
			plot_num += 1

			# Plot 50s total counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s50_total_counts, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("50S Total Counts", fontsize="x-small")
			plt.title("50S Total Counts")
			plot_num += 1

			# Plot 30s and 50s total counts together
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s30_total_counts, color="blue", label = "30S")
			plt.plot(time / 60., s50_total_counts, color="red", label = "50S")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Total Counts", fontsize="x-small")
			plt.title("30S and 50S Total Counts")
			plt.legend()
			plot_num += 1

			# Plot 30s limiting protein, total, and rRNA counts together
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s30_limiting_protein_counts, color="black", label = "Limiting Protein")
			plt.plot(time / 60., s30_16s_rRNA_total_counts, color="blue", label = "16S rRNA")
			plt.plot(time / 60., s30_total_counts, color="red", ls='--', label = "Total", lw=0.5)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Counts", fontsize="x-small")
			plt.title("30S Limiting Protein, Total, and 16S rRNA Counts")
			plt.legend()
			plot_num += 1

			# Plot 50s limiting protein, total, and rRNA counts together
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			plt.axvline(
				(time[doubling_time_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.axvline(
				(time[analysis_gen_index] / 60.), color='gray',
				linestyle='--', lw=0.5)
			plt.plot(time / 60., s50_limiting_protein_counts, color="black", label = "Limiting Protein")
			plt.plot(time / 60., s50_23s_rRNA_total_counts, color="blue", label = "23S rRNA")
			plt.plot(time / 60., s50_5s_rRNA_total_counts, color="green", label = "5S rRNA")
			plt.plot(time / 60., s50_total_counts, color="red", ls='--', label = "Total", lw=0.5)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0, 10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Counts", fontsize="x-small")
			plt.title("50S Limiting Protein, Total, and rRNA Counts")
			plt.legend()
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
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.axvline(
					(time[doubling_time_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.axvline(
					(time[analysis_gen_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.plot(time / 60., rrna_gene_copy_numbers[:, i], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((0, 10000))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("Copy Number", fontsize="x-small")
				plt.title("Copy Number of " + GENE_ID_TO_RRNA_OPERON_ID[gene_id])
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
			for i, gene_id in enumerate(GENE_ID_TO_RRNA_OPERON_ID.keys()):
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.axvline(
					(time[doubling_time_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.axvline(
					(time[analysis_gen_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.plot(time / 60., rrna_rna_target_synth_prob[:, i], color="gray", label = "Target")
				plt.plot(time / 60., rrna_actual_synth_prob[:, i], color=LINE_COLOR, label = "Actual")
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((-0.1, 0.5))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("RNA Synthesis Probability", fontsize="x-small")
				plt.title("RNA Synthesis Probability of " + GENE_ID_TO_RRNA_OPERON_ID[gene_id])
				plt.legend()
				plot_num += 1

			# rRNA occupancy
			# Get number of rRNA initiation events and copy numbers
			initiation_events = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True,
				fun=lambda x: x[:, rnap_data_rrna_indexes]
				)
			copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'promoter_copy_number',
				ignore_exception=True,
				fun=lambda x: x[:, rna_synth_prob_rrna_indexes])
			# Calculate total number of initiation events per copy
			initiation_events_per_copy_each_rrna = (
				initiation_events / copy_numbers
				)

			# Plot initiation events per copy for each rRNA operon
			for i, gene_id in enumerate(GENE_ID_TO_RRNA_OPERON_ID.keys()):
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.axvline(
					(time[doubling_time_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.axvline(
					(time[analysis_gen_index] / 60.), color='gray',
					linestyle='--', lw=0.5)
				plt.plot(time / 60., initiation_events_per_copy_each_rrna[:, i], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((-0.1, 0.5))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("Initiation Events Per Copy", fontsize="x-small")
				plt.title("Initiation Events Per Copy of " + GENE_ID_TO_RRNA_OPERON_ID[gene_id])
				plot_num += 1

			print("Total number of plots made: ", plot_num)
			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
