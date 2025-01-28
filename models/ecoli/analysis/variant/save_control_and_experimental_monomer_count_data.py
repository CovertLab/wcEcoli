"""
Anaylsis Script for saving the average total monomer counts for variants.
Note that this file only works for saving data for simulations that contain at
least two variants (one control and at least one experimental variant).
"""

import pickle
import os
import csv
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from wholecell.utils.protein_counts import get_simulated_validation_counts
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.io.tablereader import TableReader

""" USER INPUTS """

#Indicate which generation(s) to skip within each seed (sometimes this number
# should be greater than 0 because the first few generations may not be
# representative of the true dynamics occuring in the cell).
IGNORE_FIRST_N_GENS = 1 # use 14 for Sherlock runs with 24 gens

# Set this number to the desired filter threshold in which the data must pass
# through to be plotted. For example, with a filter of 1, only proteins with
# average total monomer counts above 1 in both the control and experimental
# variant will be plotted. Can be equal to 0 or greater.
FILTER = 0

""" END USER INPUTS"""

def save_file(out_dir, filename, columns, values):
	"""
	Saves data to a csv file
	Args:
		out_dir: directory to save to
		filename: name of file to be saved
		columns: columns within the table to be saved
		values: values to be saved

	Returns: A csv file with desired the data
	"""
	output_file = os.path.join(out_dir, filename)
	print(f'Saving data to {output_file}')
	with open(output_file, 'w') as f:
		writer = csv.writer(f)

		# Header for columns
		writer.writerow(columns)

		# Data rows
		for i in range(values.shape[0]):
			writer.writerow(values[i,:])

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
				self.ap.get_cells(variant=[variant],
								  generation=np.arange(IGNORE_FIRST_N_GENS,
													   self.n_total_gens),
								  only_successful=True))

			# Get the average total monomer counts in each cell generation:
			average_total_counts = (
				read_stacked_columns(all_cells,'MonomerCounts',
									 'monomerCounts',
									 fun=lambda x: np.mean(x[:], axis=0)))
			# Average together over all generations and seeds:
			avg_total_monomer_counts = np.mean(average_total_counts, axis=0)

			# Define the average total monomer counts for all proteins:
			self.total_protein_counts[var_idx] = avg_total_monomer_counts

			# Extract the protein counts for the original/native proteins:
			old_gene_idxs = [monomer_idx_dict.get(monomer_id)
							 for monomer_id in self.original_monomer_ids]
			avg_native_monomer_counts = avg_total_monomer_counts[old_gene_idxs]
			protein_counts[var_idx] = avg_native_monomer_counts

		# Return the protein counts for all proteins and the original proteins:
		self.total_protein_counts = np.array(self.total_protein_counts)
		protein_counts = np.array(protein_counts)

		return protein_counts

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
			pass  # if no filter, then return the data as is.
		else:
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

	def transpose_data(self, data):
		"""
		Transpose the data to obtain the desired format for saving
		Args:
			data: data to be transposed

		Returns: the transposed data
		"""
		data = [np.array(data)]; data = np.transpose(data)
		return data


	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Create save paths:
		outDir = plotOutDir[:-8] # remove part of the path
		save_pth = os.path.join(outDir, "variant_average_monomer_count_data")
		pth = os.path.join(save_pth, "unfiltered_data")
		log_pth = os.path.join(pth, "log_data")
		F_pth = os.path.join(save_pth, "filtered_data")
		F_log_pth = os.path.join(F_pth, "log_data")

		if not os.path.exists(save_pth):
			os.mkdir(save_pth)
			os.mkdir(pth)
			os.mkdir(F_pth)
			os.mkdir(log_pth)
			os.mkdir(F_log_pth)

		# Get data for all variants
		all_variants = self.ap.get_variants()
		control_var = all_variants[0] # control variant
		experimental_vars = all_variants[1:] # experimental variants
		for variant_2 in experimental_vars:
			self.variant_pair = [control_var, variant_2]
			experimental_var = self.variant_pair[1]

			# define/initialize commonly used variables
			self.n_total_gens = self.ap.n_generation
			self.all_monomer_ids = []
			self.original_monomer_ids = []
			self.new_gene_monomer_ids = []
			self.total_protein_counts = []
			protein_counts = self.generate_data(simDataFile)
			monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
										  enumerate(self.all_monomer_ids)}

			# Manditory Data Filtration
			F_PCs, F_PC_ids, F_PC_idxs = (
				self.filter_data(protein_counts, monomer_idx_dict_PreFilter,
								 FILTER))

			# Save unfiltered data:
			ids = self.transpose_data(self.all_monomer_ids)
			PCs_current = self.total_protein_counts
			values = np.concatenate((ids, PCs_current.T), axis=1)
			expstr = ("Variant" + str(experimental_var) +
					  " Average Monomer Counts")
			col_labels = ["Monomer IDs", "Control Average Monomer Counts", expstr]
			save_file(pth,f'AvgProteinCounts_startGen_'
						  f'{IGNORE_FIRST_N_GENS}_Variant{experimental_var}.csv',
					  col_labels, values)

			# Save unfiltered log data:
			log_PCs = self.get_LogData(self.all_monomer_ids,
									   self.total_protein_counts)
			log_PC_ids = self.transpose_data(self.all_monomer_ids)
			log_values = np.concatenate((log_PC_ids, log_PCs.T), axis=1)
			log_col_labels = \
				["Monomer ID", "Control Log10 Average Monomer Counts",
				 f"Variant{experimental_var} Log10 Average Monomer Counts"]
			save_file(log_pth,
					  f'LogAvgProteinCounts_startGen_'
						f'{IGNORE_FIRST_N_GENS}_Variant{experimental_var}.csv',
					  log_col_labels, log_values)


			# Save filtered data:
			Fcol_labels = ["Filtered Monomer IDs",
						   "Control Average Monomer Counts", expstr]
			F_PC_ids = self.transpose_data(F_PC_ids)
			Fvalues = np.concatenate((F_PC_ids, F_PCs.T), axis=1)
			save_file(F_pth,
				f'AvgProteinCounts_startGen_{IGNORE_FIRST_N_GENS}_'
					  f'filter_{FILTER}_Variant{experimental_var}.csv',
				Fcol_labels, Fvalues)

			# Save filtered log data:
			F_log_PCs = self.get_LogData(F_PC_idxs, F_PCs)
			F_log_values = (
				np.concatenate((F_PC_ids, F_log_PCs.T), axis=1))
			F_log_col_labels = \
				["Filtered Monomer ID", "Control Log10 Average Monomer Counts",
				 f"Variant{experimental_var} Log10 Average Monomer Counts"]
			save_file(F_log_pth, f'LogAvgProteinCounts_startGen_'
								 f'{IGNORE_FIRST_N_GENS}_filter_{FILTER}'
								 f'_Variant{experimental_var}.csv',
					  F_log_col_labels, F_log_values)

if __name__ == "__main__":
	Plot().cli()
