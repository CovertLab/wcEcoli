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

# OTHER_GENES_OF_INTEREST = [
# 	"EG11024 ", # trpA, TRYPSYN-APROTEIN[c]
# 	"EG11025", # trpB, TRYPSYN-BPROTEIN[c]
# 	"EG10447", # hisD, HISTDEHYD-MONOMER[c]
# 	"EG11000 ", # thrC, THRESYN-MONOMER[c]
# ]

OTHER_MONOMERS_OF_INTEREST = [
	"TRYPSYN-APROTEIN[c]",
	"TRYPSYN-BPROTEIN[c]",
	"HISTDEHYD-MONOMER[c]",
	"THRESYN-MONOMER[c]",
]

OTHER_GENE_COMMON_NAMES = [ # TODO: extract these using sim_data
	"trpA",
	"trpB",
	"hisD",
	"thrC",
]

OTHER_COMPLEXES_OF_INTEREST = ["TRYPSYN[c]"]

# TODO: also plot tryptophan synthase complex? ID: TRYPSYN

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
		mRNA_monomer_id_dict = dict(zip(
			monomer_sim_data['cistron_id'], monomer_sim_data['id']))
		monomer_to_mRNA_id_dict = dict(zip(
			monomer_sim_data['id'],monomer_sim_data['cistron_id']))
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

		other_gene_mRNA_ids = {}
		other_gene_mRNA_indexes = {}
		for monomer in OTHER_MONOMERS_OF_INTEREST:
			other_gene_mRNA_ids[monomer] = list(self.get_mRNA_ids_from_monomer_ids(
				sim_data, [monomer]))
			other_gene_mRNA_indexes[monomer] = self.get_mRNA_indexes_from_monomer_ids(
				sim_data, cell_paths, [monomer], 'mRNA')

		# determine cistron ids from otehr_monomers_of_interest
		other_gene_cistron_ids = [
			monomer_to_mRNA_id_dict.get(monomer_id) for monomer_id
			in OTHER_MONOMERS_OF_INTEREST]
		new_gene_cistron_ids = [
			monomer_to_mRNA_id_dict.get(monomer_id) for monomer_id in new_gene_monomer_ids]
		# Extract cistron indexes for each gene of interest
		cistron_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		cistron_idx_dict = {cistron: i for i, cistron in enumerate(
			cistron_counts_reader.readAttribute('mRNA_cistron_ids'))}
		other_gene_cistron_indexes = [
			cistron_idx_dict.get(cistron_id) for cistron_id in other_gene_cistron_ids]
		cistron_counts_reader.close()

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
		other_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									  monomer_id in OTHER_MONOMERS_OF_INTEREST]
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

		# Other genes of interest
		other_gene_monomer_counts = read_stacked_columns(
			cell_paths, 'MonomerCounts', 'monomerCounts', ignore_exception=True)[:,other_gene_monomer_indexes]
		all_mRNA_cistron_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts', ignore_exception=True)
		other_gene_cistron_counts = all_mRNA_cistron_stacked_counts[:,other_gene_cistron_indexes]

		plot_suffixes = ["", "_standard_axes_y", "_standard_axes_both"]
		standard_xlim = (0,2000)
		total_plots = 50 # TODO Modularize and get rid of this magic number

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			mpl.rcParams['axes.spines.right'] = False
			mpl.rcParams['axes.spines.top'] = False
			plt.figure(figsize = (8.5, 100))
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

			# New Gene RNA Synthesis Probability
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			new_gene_target_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
				ignore_exception=True, remove_first=True)[:, new_gene_RNA_indexes]
			new_gene_actual_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
				ignore_exception=True, remove_first=True)[:, new_gene_RNA_indexes]
			plt.plot(time_no_first / 60., new_gene_target_rna_synth_prob,
					 label="Target")
			plt.plot(time_no_first / 60., new_gene_actual_rna_synth_prob,
					 label="Actual")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1,0.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("RNA Synthesis Probability", fontsize='x-small')
			plt.title("New Gene RNA Synthesis Probability")
			plt.legend()
			plot_num += 1

			# Amino acid concentrations
			all_aa_ids = sim_data.molecule_groups.amino_acids
			# aa_ids_of_interest = ['HIS[c]', 'THR[c]', 'TRP[c]']
			aa_ids_of_interest = all_aa_ids
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

			# Other genes of interest
			for ii, gene in enumerate(OTHER_GENE_COMMON_NAMES):
				# RNA Synthesis Probability
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				for r in range(len(other_gene_mRNA_ids[OTHER_MONOMERS_OF_INTEREST[ii]])):
					other_gene_target_rna_synth_prob = read_stacked_columns(
						cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
						ignore_exception=True, remove_first=True)[:, other_gene_mRNA_indexes[OTHER_MONOMERS_OF_INTEREST[ii]][r]]
					other_gene_actual_rna_synth_prob = read_stacked_columns(
						cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
						ignore_exception=True, remove_first=True)[:, other_gene_mRNA_indexes[OTHER_MONOMERS_OF_INTEREST[ii]][r]]

					plt.plot(time_no_first / 60., other_gene_actual_rna_synth_prob,
							 label=other_gene_mRNA_ids[OTHER_MONOMERS_OF_INTEREST[ii]][r] + ": Actual")
					plt.plot(time_no_first / 60., other_gene_target_rna_synth_prob,
							 label=other_gene_mRNA_ids[OTHER_MONOMERS_OF_INTEREST[ii]][r] + ": Target")
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((-0.1,0.5))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("RNA Synthesis Probability", fontsize='x-small')
				plt.title("RNA Synthesis Probability: " + gene)
				plt.legend(loc='lower right')
				plot_num += 1

				# Cistron counts
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time / 60., other_gene_cistron_counts[:, ii], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((0, 1000))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("Cistron Counts", fontsize='small')
				plt.title("Cistron Counts: " + gene)
				plot_num += 1

				# Protein counts
				plt.subplot(total_plots, 1, plot_num, sharex=ax1)
				plt.plot(time / 60., other_gene_monomer_counts[:,ii], color=LINE_COLOR)
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
					plt.ylim((0, 1000))
				if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
					plt.xlim(standard_xlim)
				plt.xlabel("Time (min)")
				plt.ylabel("Protein Counts", fontsize='small')
				plt.title("Protein Counts: " + gene)
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
			#
			# # New Gene RNA Synthesis Probability
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# new_gene_target_rna_synth_prob = read_stacked_columns(
			# 	cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
			# 	ignore_exception=True, remove_first=True)[:, new_gene_RNA_indexes]
			# new_gene_actual_rna_synth_prob = read_stacked_columns(
			# 	cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
			# 	ignore_exception=True, remove_first=True)[:, new_gene_RNA_indexes]
			#
			# if len(new_gene_mRNA_ids) == 1:
			# 	plt.plot(time_no_first / 60., new_gene_target_rna_synth_prob,
			# 			 label="Target")
			# 	plt.plot(time_no_first / 60., new_gene_actual_rna_synth_prob,
			# 			 label="Actual")
			# else:
			# 	for r in range(len(new_gene_mRNA_ids)):
			# 		plt.plot(time_no_first / 60., new_gene_target_rna_synth_prob[:,r],
			# 				 label = new_gene_mRNA_ids[r] + ": Target")
			# 		plt.plot(time_no_first / 60., new_gene_actual_rna_synth_prob[:, r],
			# 				 label=new_gene_mRNA_ids[r] + ": Actual")
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	plt.ylim((-0.1,0.5))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	plt.xlim(standard_xlim)
			# plt.xlabel("Time (min)")
			# plt.ylabel("RNA Synthesis Probability", fontsize='x-small')
			# plt.title("New Gene RNA Synthesis Probability")
			# plt.legend()
			# plot_num += 1
			#
			# # New Gene Protein Initialization Probability
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# new_gene_target_protein_init_prob = read_stacked_columns(
			# 	cell_paths, 'RibosomeData', 'target_prob_translation_per_transcript',
			# 	ignore_exception=True, remove_first=True)[:, new_gene_monomer_indexes]
			# new_gene_actual_protein_init_prob = read_stacked_columns(
			# 	cell_paths, 'RibosomeData', 'actual_prob_translation_per_transcript',
			# 	ignore_exception=True, remove_first=True)[:, new_gene_monomer_indexes]
			#
			# if len(new_gene_monomer_ids) == 1:
			# 	plt.plot(time_no_first / 60., new_gene_target_protein_init_prob,
			# 			 label="Target")
			# 	plt.plot(time_no_first / 60., new_gene_actual_protein_init_prob,
			# 			 label="Actual")
			# else:
			# 	for r in range(len(new_gene_monomer_ids)):
			# 		plt.plot(time_no_first / 60., new_gene_target_protein_init_prob[:,r],
			# 				 label = new_gene_monomer_ids[r] + ": Target")
			# 		plt.plot(time_no_first/ 60., new_gene_actual_protein_init_prob[:, r],
			# 				 label=new_gene_monomer_ids[r] + ": Actual")
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	plt.ylim((-0.1, 1.1))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	plt.xlim(standard_xlim)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Probability Translation Per Transcript", fontsize='x-small')
			# plt.title("New Gene Protein Initialization Probability")
			# plt.legend()
			# plot_num += 1
			#
			# # mRNA Counts, Gene Promoter Copy Number, and RNA Synth Prob
			# ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# ax3 = ax2.twinx()
			# # plot on log scale instead
			# ax2.plot(time_no_first / 60., new_gene_target_rna_synth_prob,
			# 		 label="Target")
			# ax2.plot(time_no_first / 60., new_gene_actual_rna_synth_prob,
			# 		 label="Actual")
			# ax3.plot(time / 60., np.log10(new_gene_mRNA_counts + 1), label = "log10(mRNA counts + 1)", color = "cyan")
			# ax3.plot(time / 60., 2 * new_gene_promoter_copy_numbers, label="Copy Number", color = "red")
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	ax2.set_ylim(((-0.1, 0.5)))
			# 	ax3.set_ylim((-1,4.5))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	ax2.set_xlim(standard_xlim)
			# ax2.set_xlabel("Time (min)")
			# plt.title("New Gene mRNA Counts, Copy Number, and RNA Synth Prob")
			# plot_num += 1
			#
			# # mRNA Counts and RNAP Counts
			# ax2 = plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# ax3 = ax2.twinx()
			# # plot on log scale instead
			# ax2.plot(time / 60., total_rnap_counts)
			# ax3.plot(time / 60., np.log10(new_gene_mRNA_counts + 1), label = "log10(mRNA counts + 1)", color = "cyan")
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	ax2.set_ylim((0,10000))
			# 	ax3.set_ylim((-1,4.5))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	ax2.set_xlim(standard_xlim)
			# ax2.set_xlabel("Time (min)")
			# plt.title("New Gene mRNA Counts and RNAP Counts")
			# plot_num += 1
			#
			# # TODO: New Gene mRNA mass fraction
			#
			# # TODO: New Gene Protein Mass Fraction
			#
			# # Active RNAP Portion Allocation
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# # Active RNAP Counts
			# uniqueMoleculeCounts = TableReader(
			# 	os.path.join(simOutDir, "UniqueMoleculeCounts"))
			# active_rnap_index = uniqueMoleculeCounts.readAttribute(
			# 	"uniqueMoleculeIds").index('active_RNAP')
			# active_rnap_counts = read_stacked_columns(
			# 	cell_paths, 'UniqueMoleculeCounts',
			# 	'uniqueMoleculeCounts',
			# 	ignore_exception=True)[:, active_rnap_index]
			# # New Gene RNAP Portion
			# new_gene_mRNA_indexes = new_gene_mRNA_indexes
			# new_gene_rnap_counts = read_stacked_columns(
			# 	cell_paths, "RNACounts", "partial_mRNA_counts",
			# 	ignore_exception=True)[:, new_gene_mRNA_indexes].flatten()
			# new_gene_rnap_portion = new_gene_rnap_counts / active_rnap_counts
			# # rRNA RNAP Portion
			# rrna_rnap_counts = np.sum(read_stacked_columns(
			# 	cell_paths, "RNACounts", "partial_rRNA_counts",
			# 	ignore_exception=True), axis = 1).flatten()
			# rrna_rnap_portion = rrna_rnap_counts / active_rnap_counts
			# # RNAP Subunit RNAP Portion
			# RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			# rnap_subunit_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
			# 	sim_data, cell_paths, RNAP_subunit_monomer_ids, "mRNA")
			# rnap_subunit_rnap_counts = np.sum(read_stacked_columns(
			# 	cell_paths, "RNACounts", "partial_mRNA_counts",
			# 	ignore_exception=True)[:, rnap_subunit_mRNA_indexes], axis=1).flatten()
			# rnap_subunit_rnap_portion = rnap_subunit_rnap_counts / active_rnap_counts
			# # Ribosomal Proteins RNAP Portion
			# ribosomal_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
			# ribosomal_mRNA_indexes = self.get_mRNA_indexes_from_monomer_ids(
			# 	sim_data, cell_paths, ribosomal_monomer_ids, "mRNA")
			# ribosomal_rnap_counts = np.sum(read_stacked_columns(
			# 	cell_paths, "RNACounts", "partial_mRNA_counts",
			# 	ignore_exception=True)[:, ribosomal_mRNA_indexes], axis=1).flatten()
			# ribosomal_rnap_portion = ribosomal_rnap_counts / active_rnap_counts
			# # Plot
			# plt.plot(time / 60., rnap_subunit_rnap_portion, label="RNAP Subunit")
			# plt.plot(time / 60., ribosomal_rnap_portion, label="Ribo. Prot.")
			# plt.plot(time / 60., new_gene_rnap_portion, label = "New Gene")
			# plt.plot(time / 60., rrna_rnap_portion, label="rRNA")
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	plt.ylim((-0.1, 1.1))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	plt.xlim(standard_xlim)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Portion of Active RNAPs", fontsize='x-small')
			# plt.title("Allocation of Active RNAPs")
			# plt.legend(fontsize="x-small")
			# plot_num += 1
			#
			# # Active Ribosome Portion Allocation
			# plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# # Active Ribosome Counts
			# unique_molecule_counts_table = TableReader(
			# 	os.path.join(simOutDir, "UniqueMoleculeCounts"))
			# ribosome_index = unique_molecule_counts_table.readAttribute(
			# 	"uniqueMoleculeIds").index('active_ribosome')
			# active_ribosome_counts = read_stacked_columns(
			# 	cell_paths, 'UniqueMoleculeCounts',
			# 	'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# # New Gene Ribosome Portion
			# new_gene_ribosome_counts = read_stacked_columns(
			# 	cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
			# 	ignore_exception=True)[:, new_gene_monomer_indexes].flatten()
			# new_gene_ribosome_portion = new_gene_ribosome_counts / active_ribosome_counts
			# # RNAP Subunits Ribosome Portion
			# RNAP_subunit_monomer_ids = sim_data.molecule_groups.RNAP_subunits
			# rnap_subunit_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
			# 	sim_data, cell_paths, RNAP_subunit_monomer_ids, "monomer")
			# rnap_subunit_ribosome_counts = np.sum(read_stacked_columns(
			# 	cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
			# 	ignore_exception=True)[:, rnap_subunit_monomer_indexes], axis=1).flatten()
			# rnap_subunit_ribosome_portion = rnap_subunit_ribosome_counts / active_ribosome_counts
			# # Ribosomal Proteins Ribosome Portion
			# ribosomal_monomer_ids = sim_data.molecule_groups.ribosomal_proteins
			# ribosomal_monomer_indexes = self.get_mRNA_indexes_from_monomer_ids(
			# 	sim_data, cell_paths, ribosomal_monomer_ids, "monomer")
			# ribosomal_ribosome_counts = np.sum(read_stacked_columns(
			# 	cell_paths, "RibosomeData", "n_ribosomes_per_transcript",
			# 	ignore_exception=True)[:, ribosomal_monomer_indexes], axis=1).flatten()
			# ribosomal_ribosome_portion = ribosomal_ribosome_counts / active_ribosome_counts
			# # Plot
			# plt.plot(time / 60., rnap_subunit_ribosome_portion, label="RNAP Subunit")
			# plt.plot(time / 60., ribosomal_ribosome_portion, label="Ribo. Prot.")
			# plt.plot(time / 60., new_gene_ribosome_portion, label = "New Gene")
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
			# 	plt.ylim((-0.1, 1.1))
			# if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
			# 	plt.xlim(standard_xlim)
			# plt.xlabel("Time (min)")
			# plt.ylabel("Portion of Active Ribosomes", fontsize='x-small')
			# plt.title("Allocation of Active Ribosomes")
			# plt.legend(fontsize="x-small")
			# plot_num += 1

			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
