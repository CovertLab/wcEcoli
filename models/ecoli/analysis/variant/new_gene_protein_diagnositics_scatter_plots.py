"""
Generate scatter plots of new gene and protein diagnostics for each variant.
"""

import pickle
import os
import pandas as pd
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


""" USER INPUTS """

# Indicate the number of generations to be ignored at the start of each seed:
IGNORE_FIRST_N_GENS = 14 # 2 for local, 14 for Sherlock (w/ 24 total gens)

# Set this value to 1 to create this comparison graph, 0 otherwise:
plot_protein_count_comparisons = 1

# Set this to 1 if you wish to create the scatter plot on a log scale as well:
plot_protein_count_comparisons_log_scale = 1

# Set this number to the desired filter threshold in which the data must pass
# through to be plotted. For example, with a filter of 1, only proteins with
# average total monomer counts above 1 in both the control and experimental
# variant will be plotted. Can be equal to 0 or greater.
FILTER = 30

# Set this number equal to 1 to create comparison plots with the filtered data:
plot_comparisons_with_filtered_data = 1

""" END USER INPUTS """

class Plot(variantAnalysisPlot.VariantAnalysisPlot):

	def generate_data(self, simDataFile):
		"""
		Extracts protein count (PC) data from the simulation data file
		Args:
			simDataFile: simulation data file

		Returns:
			protein_counts: protein count (PC) data for native proteins
			 present on the E.coli chromosome for each variant
			  (the PC for each protein is averaged over all the generations)
			self.total_protein_counts: the original PCs and new gene (NG) PCs
			in one variable
			self.new_gene_monomer_ids: protein ids for new genes inserted into
			the E.coli genome
			self.original_gene_ids: protein ids for the original proteins on
			the E.coli genome
			self.all_monomer_ids: list of all the monomer ids (NG protein ids
			and native proteins' ids)
		"""
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = (
			sim_data.process.transcription.cistron_data.struct_array)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))

		# Extract new gene(s) data:
		new_gene_mRNA_ids = mRNA_sim_data[
			mRNA_sim_data['is_new_gene']]['id'].tolist()
		self.new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
									for mRNA_id in new_gene_mRNA_ids]

		# Extract data for all genes and the native genes (lacking new genes):
		self.all_monomer_ids = monomer_sim_data['id']
		self.original_monomer_ids = np.delete(self.all_monomer_ids, np.where(
			self.all_monomer_ids == self.new_gene_monomer_ids))

		# Initialize protein counts for all proteins (including new genes):
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(self.all_monomer_ids)}
		self.total_protein_counts = np.zeros((len(self.variant_pair),
											  len(self.all_monomer_ids)))

		# Initialize protein counts for the native proteins:
		protein_counts = np.zeros((len(self.variant_pair),
								   len(self.original_monomer_ids)))

		# Extract counts for proteins in each variant:
		for var_idx in range(len(self.variant_pair)):
			variant = self.variant_pair[var_idx]
			all_cells = (
				self.ap.get_cells(
					variant=[variant], generation=np.arange(IGNORE_FIRST_N_GENS,
					self.n_total_gens), only_successful=True))

			# Get the average total monomer counts over all cell durations:
			average_total_counts = (
				read_stacked_columns(all_cells,'MonomerCounts',
				'monomerCounts')).mean(axis=0)

			# Define the average total monomer counts for all proteins:
			self.total_protein_counts[var_idx] = average_total_counts

			# Extract the protein counts for the original/native proteins:
			old_gene_idxs = [monomer_idx_dict.get(monomer_id)
							 for monomer_id in self.original_monomer_ids]
			avg_native_monomer_counts = average_total_counts[old_gene_idxs]
			protein_counts[var_idx] = avg_native_monomer_counts

		# Return the protein counts for all proteins and the original proteins:
		self.total_protein_counts = np.array(self.total_protein_counts)
		protein_counts = np.array(protein_counts)

		return protein_counts

	def get_idxs(self, monomer_names):
		"""
		Obtain the indexes for specific proteins using their monomer IDs
		Args:
			monomer_names: a list of monomer names for desired proteins

		Returns: the indexes of the desired proteins within the list of all
		proteins in the chromosome (unfiltered proteins)
		"""
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(self.all_monomer_ids)}
		protein_idxs = [monomer_idx_dict.get(monomer_id)
						for monomer_id in monomer_names]
		protein_idxs = np.array(protein_idxs)

		return protein_idxs

	def get_ids(self, monomer_idx_dict, protein_idxs):
		"""
		Obtain the monomer IDs for proteins based from their indexes
		Args:
			monomer_idx_dict: a dictionary that maps protein names to their idx
			protein_idxs: an array of indices for proteins of interest

		Returns: the corresponding id for each respective protein index
		"""
		inv_monomer_idx_dict = {idx: i for i, idx in monomer_idx_dict.items()}
		protein_ids = [
			inv_monomer_idx_dict.get(monomer_id) for monomer_id in protein_idxs]
		return protein_ids

	def get_LogData(self, protein_idxs, interest_protein_counts, index_vals=[]):
		"""
		Covert monomer count data to a log10 scale
		Args:
			protein_idxs: an array of the indices for proteins to be plotted
             (this should be smaller than interest_protein_counts if the data
              has been filtered)
            interest_protein_counts: the full data structure of all proteins
            and their respective counts (no filter applied)
            (usually size variants by # of proteins), either filtered or not
            index_vals: if the protein idxs are not in sequential order (usually
            happens after filtering the data), include a list of the original
            indices for each protein in this variable.

		Returns: the log value of the total average counts for each protein
		"""
		# Initialize the array to store the log values of the protein counts:
		avg_log_interest_proteins = np.zeros((
			len(self.variant_pair), len(protein_idxs)))

		# Calculate the log values of the protein counts:
		for variant in range(len(self.variant_pair)):
			for idx in range(len(protein_idxs)):
				if len(index_vals) == 0:
					index = idx
				else:
					index = index_vals[idx]
				avg_log_interest_proteins[variant][idx] = \
					np.log10(interest_protein_counts[variant][index] + 1)

		return avg_log_interest_proteins

	def filter_data(self, nonfiltered_protein_counts,
					nonfiltered_monomer_idx_dict, filter_num=0):
		"""
		Filters the data to extract all proteins at least a specified number of
		protein counts in both variants.
		Args:
			nonfiltered_protein_counts: array of counts for all proteins
			nonfiltered_monomer_idx_dict: dictionary that maps ID to index
			filter_num: the minimum number of counts a protein must have in both
			variants to avoid being discarded from the new filtered data
			(by default this is set to zero).

		Returns: arrays of counts, IDs, and indexes for all proteins that
        successfully pass the user defined filter (or the default set to 0)
		"""
		# Extract all proteins with non-zero protein counts in both variants:
		nonzero_counts_var0_idxs = np.nonzero(nonfiltered_protein_counts[0])
		nonzero_counts_var1_idxs = np.nonzero(nonfiltered_protein_counts[1])
		shared_nonzero_PCs_idxs = np.intersect1d(nonzero_counts_var0_idxs,
												 nonzero_counts_var1_idxs)
		nonzero_PCs = nonfiltered_protein_counts[:, shared_nonzero_PCs_idxs]
		nonzero_idxs = shared_nonzero_PCs_idxs
		nonzero_ids = self.get_ids(nonfiltered_monomer_idx_dict,
									   shared_nonzero_PCs_idxs)

		if filter_num == 0:
			pass # if no filter, then return the data as is.
		else:
			# if there is a filter, reassign the variables accordingly:
			var0_filter_PCs_idxs = np.nonzero(nonzero_PCs[0] > filter_num)
			var1_filter_PCs_idxs = np.nonzero(nonzero_PCs[1] > filter_num)
			shared_filtered_PC_idxs = np.intersect1d(var0_filter_PCs_idxs,
													 var1_filter_PCs_idxs)
			# reassign variable values if the filter is applied:
			nonzero_PCs = nonzero_PCs[:, shared_filtered_PC_idxs]
			nonzero_ids = np.array(nonzero_ids)
			nonzero_ids = nonzero_ids[shared_filtered_PC_idxs]
			nonzero_idxs = shared_filtered_PC_idxs

		return nonzero_PCs, nonzero_ids, nonzero_idxs


	def generate_plot(self, var_num, all_protein_counts,
					  original_protein_counts, filtered=0):
		"""
		Generates a plot of the control variant's total monomer counts
		plotted against the experimental variant's total monomer counts
		Args:
			all_protein_counts: total protein counts (for new + native genes)
			original_protein_counts: total protein counts for native genes
			filtered: set to 1 if the data is filtered
			filter: set equal to the filter number used to filter the data

		Returns: an x-y style comparison plot.
		"""
		# Define protein counts for each variant:
		var0_x = all_protein_counts[0] # control variant, APCs
		var1_y = all_protein_counts[1] # experimental variant, APCs
		ori_var0_x = original_protein_counts[0] # control variant, NPCs
		ori_var1_y = original_protein_counts[1] # experimental variant, NPCs

		# Generate the plot:
		plt.figure(figsize=(10, 10))

		# find the location of the new gene's PCs within the dataframe:
		new_gene_PCs = var1_y[-1]

		# if the data is filtered, plot the filtered data:
		if filtered == 1:
			# Determine the maximum value for the x and y axes:
			max_val = [max(ori_var0_x), max(ori_var1_y)]
			max_val = max(max_val); max_val = round(max_val) + .5
			max_val = round(max_val); yxvals = range(0, max_val + 1)

			# plot the native proteins against each other with the new gene on
			# the plot but not included in linear calculations:
			plt.scatter(var0_x, var1_y, 1)
			plt.scatter(0, new_gene_PCs, 8, color="red", marker="*")

			# Calculate the linear fit: between the two variants:
			m, b = np.polyfit(var0_x, var1_y, 1)
			plt.plot(yxvals, m * yxvals + b, linewidth=.5, color='#bcbd22')
			legstr = ("linear fit: y = "+str(round(m, 3))+"x + "+str(round(b, 3)))

			# Plot an y=x line:
			plt.plot(yxvals, yxvals, linewidth=.5, linestyle="dashed",
					 color="grey", alpha=.5)

			# Plot the linear fit where the y-intercept is forced to be zero:
			x = ori_var0_x[:, np.newaxis]  # make a length array
			slope, _, _, _ = np.linalg.lstsq(x, ori_var1_y)
			plt.plot(yxvals, slope * yxvals,linewidth=.8, color='#FFA500')
			otherstr = "y = " + str(round(slope[0], 3)) + "x + 0"

			# Format the graph:
			plt.legend(["monomer count data","New Gene",legstr,"y=x",otherstr])
			plt.axis('square')
			plt.xlabel("variant 0 (control)")
			plt.ylabel(f"variant {var_num} (containing the new gene)")
			plt.title(f"The average total monomer counts for the {len(var0_x)}"
					  f" proteins in variant {var_num} \nplotted against the "
					  f"contorl variant (filter={FILTER})")
			plt.tight_layout()

		else:
			# Determine the maximum value for the x and y axes:
			max_val = [max(ori_var0_x), max(ori_var1_y)]
			max_val = max(max_val); max_val = round(max_val) + .5
			max_val = round(max_val); yxvals = range(0, max_val + 1)

			# plot the native proteins against each other with the new gene on
			# the plot but not included in linear calculations:
			plt.scatter(ori_var0_x, ori_var1_y, 1)
			plt.scatter(0, new_gene_PCs, 8, color="red", marker="*")

			# create a linear fit:
			m, b = np.polyfit(ori_var0_x, ori_var1_y, 1)
			plt.plot(yxvals, m * yxvals + b,linewidth=.5,color='#bcbd22')
			legstr = ("linear fit: y = " + str(round(m, 3)) +
						  "x + " + str(round(b, 3)))

			# Plot an y=x line:
			plt.plot(yxvals, yxvals, linewidth=.5, linestyle="dashed",
					 color="grey", alpha=.5)

			# Plot the linear fit where the y-intercept is forced to be zero:
			x = ori_var0_x[:, np.newaxis]  # make a length array
			slope, _, _, _ = np.linalg.lstsq(x, ori_var1_y)
			plt.plot(yxvals, slope * yxvals,linewidth=.8, color='#FFA500')
			otherstr = "y = " + str(round(slope[0], 3)) + "x + 0"

			# format the plot:
			#plt.axis('square')
			plt.title(f"The average total monomer counts for the {len(var0_x)}"
					  f" proteins in variant {var_num} plotted against the "
					  f"control variant")
			plt.legend(["monomer count data","New Gene",legstr,"y=x",otherstr])
			plt.xlabel("variant 0 (control)")
			plt.ylabel(f"variant {var_num} (containing the new gene)")
			plt.tight_layout()


	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		"""
		Call graph generating functions
		"""
		# Extract cell paths for each variant:
		all_variants = self.ap.get_variants()
		control_var = all_variants[0] # the first variant is the control
		experimental_vars = all_variants[1:] # remove the control variant

		# Determine the number of generations within the system:
		self.n_total_gens = self.ap.n_generation

		# Generate the plots for each variant:
		for variant_2 in experimental_vars:
			self.variant_pair = [control_var, variant_2]
			experimental_var = self.variant_pair[1]

			# define/initialize commonly used variables to be generated:
			self.all_monomer_ids = []
			self.original_monomer_ids = []
			self.new_gene_monomer_ids = []
			self.total_protein_counts = []
			protein_counts = self.generate_data(simDataFile)
			monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
										  enumerate(self.all_monomer_ids)}

			# Plots the total protein counts with no filter applied to the data:
			if plot_protein_count_comparisons == 1:
				# obtain the data to be plotted:
				PCs = self.total_protein_counts
				IDs = self.all_monomer_ids
				original_PCs = protein_counts
				original_IDs = self.original_monomer_ids

				# generate the plot:
				self.generate_plot(experimental_var, PCs, original_PCs)
				words = ('_total_monomer_count_comparison_with_Variant' +
						 str(experimental_var))
				exportFigure(plt, plotOutDir, plotOutFileName + '_' +
							 str(len(PCs[0])) + words, metadata)
				plt.close('all')

			if plot_protein_count_comparisons_log_scale == 1:
				# obtain the log data to be plotted:
				PC_LogData_idxs = self.get_idxs(IDs)
				PC_LogData = self.get_LogData(PC_LogData_idxs, PCs)
				ori_PC_LogData_idxs = self.get_idxs(original_IDs)
				ori_PC_LogData = (
					self.get_LogData(ori_PC_LogData_idxs, original_PCs))

				# generate the plot:
				self.generate_plot(experimental_var, PC_LogData, ori_PC_LogData)
				exportFigure(plt, plotOutDir, plotOutFileName + '_' +
							 str(len(PCs[0])) + words + '_LogScale', metadata)
				plt.close('all')


			# Manditory Data Filtration
			F_PCs, F_PC_ids, F_PC_idxs = (self.filter_data(protein_counts,
										monomer_idx_dict_PreFilter, FILTER))

			# Plot filtered data:
			if plot_comparisons_with_filtered_data == 1:
				# plot the filtered data:
				self.generate_plot(experimental_var, F_PCs, F_PCs, 1)
				exportFigure(plt, plotOutDir, plotOutFileName + '_' +
							 str(len(F_PCs[0])) + '_filtered_total_monomer_'
							 'count_comparison_with_Variant'+str(experimental_var)
							 + '_Filter_' + str(FILTER), metadata)
				plt.close('all')

				# plot the log scale version of the filtered data:
				F_PC_LogData = self.get_LogData(F_PC_idxs, F_PCs)
				# can do (F_PC_idxs, total_monomer_counts, F_PC_idxs) too
				self.generate_plot(experimental_var, F_PC_LogData, F_PC_LogData,
								   1)
				exportFigure(plt, plotOutDir, plotOutFileName + '_' +
							 str(len(F_PCs[0])) +'_filtered_total_monomer_count'
							 '_comparison_with_Variant' + str(experimental_var)
							 + '_Filter_' + str(FILTER) + '_LogScale', metadata)
				plt.close('all')

if __name__ == "__main__":
	Plot().cli()
