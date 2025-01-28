"""
Save protein count data from cohort simulations (with multiple seeds) as
csv files. This script is intended to be run with sherlock in order to easily
extract and save protein count data from cohort simulations.

Note: this script calculates the average total monomer counts, not the free
monomer counts.
"""

import pickle
import os
import csv
import numpy as np
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts

""" USER INPUTS """

# Indicate the number of generations to be ignored at the start of each seed:
IGNORE_FIRST_N_GENS = 2 # 2 for local, 14 for Sherlock (w/ 24 total gens)

# Set this value to 1 to create this comparison graph, 0 otherwise:
plot_protein_count_comparisons = 1

# Set this to 1 if you wish to create the scatter plot on a log scale as well:
plot_protein_count_comparisons_log_scale = 1

# Set this number to the desired filter threshold in which the data must pass
# through to be plotted. For example, with a filter of 1, only proteins with
# average total monomer counts above 1 in both the control and experimental
# variant will be plotted. Can be equal to 0 or greater.
FILTER = 0

# Set this number equal to 1 to also save the filtered data to a csv file:
save_filtered_data = 1

# Set this number equal to 1 to also save the data with corresponding
# validation data in the csv file:
save_with_validation_data = 1

""" END USER INPUTS """

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

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def generate_data(self, simDataFile):
		"""
        Extracts the average total monomer counts for each protein in the
        simulation
        Args:
            simDataFile: simulation data file

        Returns:
            total_protein_counts: average total monomer counts for all
            proteins in the simulation
            self.all_monomer_ids: list of all the monomer ids in the simulation
        """
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)

		# Extract monomer count data for each protein:
		self.all_monomer_ids = monomer_sim_data['id']
		all_cells = self.ap.get_cells(
			generation=np.arange(IGNORE_FIRST_N_GENS,self.n_total_gens),
			only_successful=True)

		# Calculate the average total monomer counts across all cells:
		total_protein_counts = (read_stacked_columns(
			all_cells, 'MonomerCounts',
			'monomerCounts')).mean(axis=0)

		# Make into a np.array:
		total_protein_counts = np.array(total_protein_counts)

		return total_protein_counts

	def get_validation_data(self, simDataFile, validationDataFile):
		# adapted from multigen and single/proteinCountsValidation.py
		"""
		Extract the protein counts for proteins that are present in either
		validation data set
		Args:
			simDataFile: Simulation data file
			validationDataFile: validation data file

		Returns: lists of proteins that are present in both the simulation and
		validation datasets, as well as the protein counts for each protein.
		"""
		# Get simulation data:
		sim_data = self.read_pickle_file(simDataFile)
		sim_monomer_ids = sim_data.process.translation.monomer_data["id"]

		# Get validation data:
		validation_data = self.read_pickle_file(validationDataFile)
		wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
		schmidt_ids = validation_data.protein.schmidt2015Data["monomerId"]
		wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		# Get the simulation cell directories:
		all_cells = self.ap.get_cells(generation=np.arange(IGNORE_FIRST_N_GENS,
														   self.n_total_gens),
									  only_successful=True)

		# Initialize lists to store data for overlapping protein counts:
		sim_schmidt_counts_multigen = []
		sim_wisniewski_counts_multigen = []

		# Extract data for proteins that are present in the validation data:
		for simDir in all_cells:
			simOutDir = os.path.join(simDir, "simOut")
			monomer_counts_reader = (
				TableReader(os.path.join(simOutDir, "MonomerCounts")))
			monomer_counts = monomer_counts_reader.readColumn("monomerCounts")

			# Obtain the protein counts for protiens that are present in both
			# the simluation and validation datasets:
			sim_schmidt_counts, val_schmidt_counts, schmidt_overlap_ids = (
				get_simulated_validation_counts(
				schmidt_counts, monomer_counts, schmidt_ids, sim_monomer_ids))
			sim_wisniewski_counts, val_wisniewski_counts, wisniewski_overlap_ids \
				= get_simulated_validation_counts(
				wisniewski_counts, monomer_counts, wisniewski_ids, sim_monomer_ids)

			# Append the protein counts for the current cell:
			sim_schmidt_counts_multigen.append(sim_schmidt_counts)
			sim_wisniewski_counts_multigen.append(sim_wisniewski_counts)

		# Average over all the cells:
		sim_schmidt_counts_multigen = (
			(np.array(sim_schmidt_counts_multigen)).mean(axis=0))
		sim_wisniewski_couts_multigen = (
			(np.array(sim_wisniewski_counts_multigen)).mean(axis=0))

		return (sim_schmidt_counts_multigen, val_schmidt_counts,
				schmidt_overlap_ids, sim_wisniewski_couts_multigen,
				val_wisniewski_counts, wisniewski_overlap_ids)

	def get_ids(self, monomer_idx_dict, protein_idxs):
		"""
        Obtain the protein ids for each protein based on its respective index
        Args:
            monomer_idx_dict: a dictionary that maps protein names to their idx
            protein_idxs: an array of indices for each protein, respectively

        Returns: the corresponding id for each respective protein index
        """
		inv_monomer_idx_dict = {idx: i for i, idx in monomer_idx_dict.items()}
		protein_ids = \
			[inv_monomer_idx_dict.get(monomer_id) for monomer_id in protein_idxs]

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

        Returns: the log10 values for the total average counts for each protein
        """
		avg_log_interest_proteins = np.zeros(len(protein_idxs))
		for idx in range(len(protein_idxs)):
			if len(index_vals) == 0:
				index = idx
			else:
				index = index_vals[idx]
			avg_log_interest_proteins[idx] = (
				np.log10(interest_protein_counts[index] + 1))
		return avg_log_interest_proteins

	def filter_data(self, nonfiltered_protein_counts,
					nonfiltered_monomer_idx_dict, filter_num=0):
		"""
        Filters the data to extract all proteins at least a specified number of
		protein counts in both variants.
        Args:
            nonfiltered_protein_counts: array of counts for all proteins
            nonfiltered_monomer_idx_dict: dictionary that maps ID to index
            filter_num: the minimum number of counts a protein must have to
            avoid being filtered out from the data (set to zero by default).

        Returns: arrays of counts, IDs, and indexes for all proteins that
        successfully pass the user defined filter (or the default set to 0)
        """
		# Extract all proteins with non-zero count values:
		nonzero_idxs = np.nonzero(nonfiltered_protein_counts)
		nonzero_PCs = nonfiltered_protein_counts[nonzero_idxs]
		nonzero_idxs = nonzero_idxs[0] # remove the extra array layer
		nonzero_ids = self.get_ids(nonfiltered_monomer_idx_dict,
									  nonzero_idxs)

		if filter_num == 0:
			# if no filter, then return the data as is.
			pass
		else:
			# reassign variable values if the filter is applied:
			filter_PC_idxs = np.nonzero(nonzero_PCs > filter_num)
			nonzero_PCs = nonzero_PCs[filter_PC_idxs]
			nonzero_ids = np.array(nonzero_ids)
			nonzero_ids = nonzero_ids[filter_PC_idxs]
			nonzero_idxs = filter_PC_idxs
		return nonzero_PCs, nonzero_ids, nonzero_idxs

	def transpose_and_reshape(self, data):
		"""
		Transpose and reshape the data to obtain the desired format for saving
		Args:
			data: data to be transposed and reshaped

		Returns: the transposed and reshaped data
		"""
		data = np.transpose(np.array(data)); data = data.reshape(-1, 1)
		return data


	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# Create the output directories for the data to be saved to:
		save_pth = os.path.join(variantDir,
								"cohort_average_monomer_count_data")
		pth = os.path.join(save_pth, "unfiltered_data")
		log_pth = os.path.join(pth, "log_data")

		if not os.path.exists(save_pth):
			os.mkdir(save_pth)
			os.mkdir(pth)
			os.mkdir(log_pth)

		# define/initialize commonly used variables:
		self.n_total_gens = self.ap.n_generation
		startGen = IGNORE_FIRST_N_GENS # note python numbering starts at 0
		self.all_monomer_ids = []
		total_protein_counts = self.generate_data(simDataFile)
		monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
									  enumerate(self.all_monomer_ids)}

		# Save the unfiltered data:
		ids = self.transpose_and_reshape(self.all_monomer_ids)
		current_counts = self.transpose_and_reshape(total_protein_counts)
		values = np.concatenate((ids, current_counts), axis=1)
		col_labels = ["Monomer ID", "Average Monomer Counts"]
		save_file(pth,f'AvgProteinCounts_startGen_{startGen}.csv',
				  col_labels, values)

		# Save the unfiltered data on a log scale:
		log_PCs = self.get_LogData(self.all_monomer_ids, total_protein_counts)
		log_PC_ids = self.transpose_and_reshape(self.all_monomer_ids)
		log_PCs = self.transpose_and_reshape(log_PCs)
		log_values = np.concatenate((log_PC_ids, log_PCs), axis=1)
		log_col_labels = ["Monomer ID", "Log10 Average Monomer Counts"]
		save_file(log_pth,f'LogAvgProteinCounts_startGen_{startGen}.csv',
				  log_col_labels, log_values)

		# Create and save filtered data
		if save_filtered_data == 1:
			# Generate the filtered data directories:
			F_pth = os.path.join(save_pth, "filtered_data")
			F_log_pth = os.path.join(F_pth, "log_data")
			if not os.path.exists(F_pth):
				os.mkdir(F_pth)
				os.mkdir(F_log_pth)

			# Conduct mandatory data filtration:
			F_PCs, F_PC_ids, F_PC_idxs = (
				self.filter_data(total_protein_counts,
								 monomer_idx_dict_PreFilter, FILTER))

			# Save filtered data:
			F_PC_ids = self.transpose_and_reshape(F_PC_ids)
			F_PCs = self.transpose_and_reshape(F_PCs)
			Fvalues = np.concatenate((F_PC_ids, F_PCs), axis=1)
			Fcol_labels = ["Filtered Monomer ID", "Average Monomer Counts"]
			save_file(F_pth,
					  f'Filtered_AvgProteinCounts_startGen_{startGen}_'
					  f'filter_{FILTER}.csv', Fcol_labels, Fvalues)

			# Save filtered log data:
			F_log_PCs = self.get_LogData(F_PC_idxs, F_PCs)
			F_log_PC_ids = self.transpose_and_reshape(F_PC_ids)
			F_log_PCs = self.transpose_and_reshape(F_log_PCs)
			F_log_values = np.concatenate(
				(F_log_PC_ids, F_log_PCs), axis=1)
			F_log_col_labels = \
				["Filtered Monomer ID", "Log10 Average Monomer Counts"]
			save_file(F_log_pth,
					f'Filtered_LogAvgProteinCounts_startGen_{startGen}'
					f'_filter_{FILTER}.csv', F_log_col_labels, F_log_values)

		# Save protein counts with available validation data:
		if save_with_validation_data == 1:
			# Generate the validation data directories:
			validation_pth = os.path.join(save_pth, "validation_data")
			validation_log_path = os.path.join(validation_pth, "log_data")
			if not os.path.exists(validation_pth):
				os.mkdir(validation_pth)
				os.mkdir(validation_log_path)

			# Obtain the validation data:
			(sim_schmidt_counts, val_schmidt_counts, s_ids, sim_wisniewski_counts,
			 val_wisniewski_counts, w_ids) = self.get_validation_data(
				 simDataFile, validationDataFile)

			# Save Wisniewski data:
			ids = self.transpose_and_reshape(w_ids)
			sim_wisniewski_counts = (
				self.transpose_and_reshape(sim_wisniewski_counts))
			val_wisniewski_counts = (
				self.transpose_and_reshape(val_wisniewski_counts))
			values = np.concatenate((ids, sim_wisniewski_counts,
									 val_wisniewski_counts), axis=1)
			col_labels = ["Monomer ID", "Simulated Counts",
						  "Wisniewski Validation Counts"]
			save_file(validation_pth,f'Wisniewski_Comparison_startGen_'
									 f'{startGen}.csv',
					  col_labels, values)

			# Save Schmidt data:
			ids = self.transpose_and_reshape(s_ids)
			sim_schmidt_counts = self.transpose_and_reshape(sim_schmidt_counts)
			val_schmidt_counts = self.transpose_and_reshape(val_schmidt_counts)
			values = np.concatenate((ids, sim_schmidt_counts,
									 val_schmidt_counts), axis=1)
			col_labels = \
				["Monomer ID", "Simulated Counts", "Schmidt Validation Counts"]
			save_file(validation_pth,f'Schmidt_Comparison_startGen_'
									 f'{startGen}.csv', col_labels, values)

			# Save log data for Wisniewski validation data:
			log_sim_wisniewski_counts = self.get_LogData(sim_wisniewski_counts,
														 sim_wisniewski_counts)
			log_val_wisniewski_counts = self.get_LogData(val_wisniewski_counts,
														 val_wisniewski_counts)
			log_ids = self.transpose_and_reshape(w_ids)
			log_sim_wisniewski_counts = (
				self.transpose_and_reshape(log_sim_wisniewski_counts))
			log_val_wisniewski_counts = (
				self.transpose_and_reshape(log_val_wisniewski_counts))
			log_values = (
				np.concatenate((log_ids, log_sim_wisniewski_counts,
								log_val_wisniewski_counts), axis=1))
			log_col_labels = ["Monomer ID", "Log10 Simulated Counts",
							  "Log10 Wisniewski Validation Counts"]
			save_file(validation_log_path,
					  f'Log10_Wisniewski_Comparison_startGen_'
					  f'{startGen}.csv', log_col_labels, log_values)

			# Save log data for Schmidt validation data:
			log_sim_schmidt_counts = self.get_LogData(sim_schmidt_counts,
													  sim_schmidt_counts)
			log_val_schmidt_counts = self.get_LogData(val_schmidt_counts,
													  val_schmidt_counts)
			log_ids = self.transpose_and_reshape(s_ids)
			log_sim_schmidt_counts = (
				self.transpose_and_reshape(log_sim_schmidt_counts))
			log_val_schmidt_counts = (
				self.transpose_and_reshape(log_val_schmidt_counts))
			log_values = np.concatenate((log_ids, log_sim_schmidt_counts,
										 log_val_schmidt_counts), axis=1)
			log_col_labels = ["Monomer ID", "Log10 Simulated Counts",
							  "Log10 Schmidt Validation Counts"]
			save_file(validation_log_path,
				f'Log10_Schmidt_Comparison_startGen_'
				f'{startGen}.csv', log_col_labels, log_values)


if __name__ == '__main__':
	Plot().cli()
