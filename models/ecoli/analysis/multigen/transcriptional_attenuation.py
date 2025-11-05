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

		plot_suffixes = [""]
		standard_xlim = (0,2000)
		total_plots = 45 # TODO Modularize and get rid of this magic number

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.figure(figsize = (8.5, 66))
			plot_num = 1
			ax1 = plt.subplot(total_plots, 1, plot_num)

			# Growth Rate
			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(500, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')
			plt.plot(time.flatten() / 60., growth_rate_convolved, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,.0004)
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Growth Rate", fontsize="small")
			plt.title("Growth Rate")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			plt.plot(time / 60., mass, color=LINE_COLOR)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,4000)
			else:
				plt.ylim(bottom=0)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)", fontsize="small")
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
			# plt.ylabel("Gene Promoter Copy Number", fontsize="small")
			# plt.title("New Gene Promoter Copy Number")
			# plot_num += 1

			# mRNA Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
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
			plt.ylabel("Log(mRNA Counts + 1)", fontsize="small")
			plt.title("New Gene mRNA Counts")
			plot_num += 1

			# Protein Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
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
			plt.ylabel("Log(Protein Counts + 1)", fontsize="small")
			plt.title("New Gene Protein Counts")
			plot_num += 1

			# Amino acid concentrations
			aa_ids_of_interest = ['HIS[c]', 'THR[c]', 'TRP[c]']
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
			plt.ylabel("ppGpp Concentration (uM)", fontsize="small")
			plt.title("ppGpp Concentration")
			plot_num += 1

			# Total RNAP Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
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
			plt.ylabel("Total RNAP Counts", fontsize="small")
			plt.title("Total RNAP Counts")
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
			plt.ylabel("Total Ribosome Counts", fontsize="small")
			plt.title("Total Ribosome Counts")
			plot_num += 1

			# Transcriptional Attenuation
			transcription_reader = TableReader(os.path.join(simOutDir, 'TranscriptElongationListener'))
			rna_ids = transcription_reader.readAttribute('rnaIds')
			attenuated_rnas = transcription_reader.readAttribute('attenuated_rnas')
			n_attenuated = len(attenuated_rnas)
			attenuated_mask = np.array([rna in set(attenuated_rnas) for rna in rna_ids])

			expected_probability = read_stacked_columns(
				cell_paths, 'TranscriptElongationListener', 'attenuation_probability')
			counts_attenuated = read_stacked_columns(
				cell_paths, 'TranscriptElongationListener', 'counts_attenuated')
			counts_synthesized = read_stacked_columns(
				cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized')[:, attenuated_mask]

			cum_attenuated = np.cumsum(counts_attenuated, axis=0)
			cum_synthesized = np.cumsum(counts_synthesized, axis=0)
			actual_probability = cum_attenuated / (cum_attenuated + cum_synthesized)

			# ids_of_interet = [
			# 	"TU3[c]", # trp
			# 	'TU00285[c]', # his
			# 	'TU00178[c]', # thr
			# ]

			for j, rna_id in enumerate(attenuated_rnas):
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time / 60.0, cum_attenuated[:, j], label = "Attenuated")
				plt.plot(time / 60.0, cum_synthesized[:, j], label = "Synthesized")
				plt.xlabel("Time (min)")
				plt.ylabel(f"Cumulative RNA Counts: {rna_id}")
				plt.legend(fontsize="x-small")
				plot_num += 1

				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time / 60.0, expected_probability[:, j], label = "Expected")
				plt.plot(time / 60.0, actual_probability[:, j], label = "Actual")
				plt.xlabel("Time (min)")
				plt.ylabel(f"Attenuation Probability: {rna_id}")
				plt.legend(fontsize="x-small")
				plot_num += 1

			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
