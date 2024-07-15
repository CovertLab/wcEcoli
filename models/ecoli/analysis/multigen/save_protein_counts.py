"""
Anaylsis Script for saving all protein counts from a given simulation. Note that
this file will average over seeds and generations for each .
"""

import pickle
import os
import csv
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.utils.protein_counts import get_simulated_validation_counts
from wholecell.io.tablereader import TableReader
#TODO: consider whether to maybe combine this script with the new gene data saving script, and just have a clause that says "if new genes are present, save this way instead"
"""
Indicate which generation the data should start being collected from (sometimes 
this number should be greater than 0 because the first few generations may not 
be representative of accurate cellular dynamics as the simulation takes a little 
time to initialize).
"""

IGNORE_FIRST_N_GENS = 1

"""
Enter a "filter number" below to filter out any proteins that have counts less
than a given filter number. This is useful for removing proteins that have very
low counts in a variant(s) that could likely be due to sub-generational gene expression. 
The default is set to 0, but can be changed to any number greater than 0.
"""
# Number to be set as the minimum threshold PC value proteins must have in both
# variants in order to be used in the plots to follow (set to 0 for default):
#todo: edit this to have a "produce_filtered_csv" question and then if that is set to one, you can specifiy a filter number (cant be the same bc filter num can be zero)
filter_num = 0

save_filtered_data = 1

#todo write an explaination here
save_with_validation_data = 1

def save_file(out_dir, filename, columns, values):
	output_file = os.path.join(out_dir, filename)
	print(f'Saving data to {output_file}')
	with open(output_file, 'w') as f:
		writer = csv.writer(f)

		# Header for columns
		writer.writerow(columns)

		# Data rows
		for i in range(values.shape[0]):
			writer.writerow(values[i,:])

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def generate_data(self, simDataFile, validationDataFile):
		"""
		Generates csv files of protein count data from simulations
		Args:
			simDataFile: simulation data file
		Returns:
			#TODO: update this description
			protein_counts: protein count (PC) data for all proteins (originally
			 present on the E.coli chromosome) in the simulation for each variant
			  (the PC for each protein is averaged over all the generations)
			self.total_protein_counts: the original PCs and new gene (NG) PCs
			in one variable
			self.new_gene_monomer_ids: protein ids for new genes inserted into
			the E.coli genome
			self.original_gene_ids: protein ids for the original proteins on
			the E.coli genome
			self.all_monomer_ids: list of all the monomer ids (NG protein ids
			and orginal proteins' gene ids)
		"""
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)
		validation_data = self.read_pickle_file(validationDataFile)
		self.all_monomer_ids = monomer_sim_data['id']
		all_cells = self.ap.get_cells(generation=np.arange(IGNORE_FIRST_N_GENS, self.n_total_gens), only_successful=True)
		# Get the protein counts for each gene/protein
		gen_avg_monomer_counts = (read_stacked_columns(all_cells, 'MonomerCounts', 'monomerCounts', fun=lambda x: np.mean(x[:], axis=0)))

		# todo: need to figure out how to get it to average across all the seeds too. (or figure out if it is automatically doing that)
		# attempt to average over seeds:
		(gen_avg_monomer_counts_seeds,) =  read_stacked_bulk_molecules(
			all_cells, self.all_monomer_ids, ignore_exception=True)
		# average gen_avg_monomer_counts_seeds along the 0th axis
		avg_monomer_counts = np.mean(gen_avg_monomer_counts_seeds, axis=0)

		# todo: make sure to check these numbers when debugging:
		total_avg_gene_counts = np.mean((gen_avg_monomer_counts), axis=0)
		total_protein_counts = np.array(total_avg_gene_counts)

		return total_protein_counts

	def get_validation_data(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		#adapted from multigen and single/proteinCountsValidation.py
		sim_data = self.read_pickle_file(simDataFile)
		validation_data = self.read_pickle_file(validationDataFile)

		monomer_ids = sim_data.process.translation.monomer_data["id"]
		sim_monomer_ids = sim_data.process.translation.monomer_data["id"]
		wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
		schmidt_ids = validation_data.protein.schmidt2015Data["monomerId"]
		wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		# Get all cells
		allDir = self.ap.get_cells(generation=np.arange(IGNORE_FIRST_N_GENS, self.n_total_gens), only_successful=True)
		sim_schmidt_counts_multigen = []
		sim_wisniewski_counts_multigen = []

		for simDir in allDir:

			simOutDir = os.path.join(simDir, "simOut")

			monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			monomer_counts = monomer_counts_reader.readColumn("monomerCounts")

			# Obtain the protein counts for protiens that are present in both the simluation and validation datasets:
			sim_schmidt_counts, val_schmidt_counts, schmidt_overlap_ids = get_simulated_validation_counts(
				schmidt_counts, monomer_counts, schmidt_ids, monomer_ids)
			sim_wisniewski_counts, val_wisniewski_counts, wisniewski_overlap_ids = get_simulated_validation_counts(
				wisniewski_counts, monomer_counts, wisniewski_ids, sim_monomer_ids)

			sim_schmidt_counts_multigen.append(sim_schmidt_counts)
			sim_wisniewski_counts_multigen.append(sim_wisniewski_counts)

		sim_schmidt_counts_multigen = (np.array(sim_schmidt_counts_multigen)).mean(axis=0)
		sim_wisniewski_couts_multigen = (np.array(sim_wisniewski_counts_multigen)).mean(axis=0)

		return sim_schmidt_counts_multigen, val_schmidt_counts, schmidt_overlap_ids, sim_wisniewski_couts_multigen, val_wisniewski_counts, wisniewski_overlap_ids


	def get_ids(self, monomer_idx_dict, protein_idxs):
		"""
		Obtain the protein ids for each protein based on its respective index
		Args:
			monomer_idx_dict: a dictionary that maps protein names to their idx
			protein_idxs: an array of indices for proteins of interest
		Returns: the corresponding id for each respective protein index
		"""
		inv_monomer_idx_dict = {idx: i for i, idx in monomer_idx_dict.items()}
		protein_ids = [inv_monomer_idx_dict.get(monomer_id) for monomer_id in protein_idxs]
		return protein_ids

	def get_LogData(self, protein_idxs, interest_protein_counts, index_vals=[]):
		"""
		todo: update full description
		Covert normal protein count data to their log10 values
		Args:
			protein_idxs: an array of the indices for proteins of interest (this should be smaller than the next entry if the data has been filtered)
			interest_protein_counts: the full data structure of all original proteins and their respective counts
			(usually size variants by # of proteins), either filtered or unfiltered
			index_vals: if the protein idxs are not in sequential order (usually happens
			after filtering the data)
		Returns: an data structure of the log version of interest PCs
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
		Filter the data to extract all proteins with 0 (or below a specified filter_num)
		protein counts (PCs)
		Args:
			nonfiltered_protein_counts: array of PCs for all proteins
			nonfiltered_monomer_idx_dict: dictionary that maps id to index
			filter_num: the minimum number of PCs a protein must have to avoid being
			discarded from the new filtered data (set to zero by default).
		Returns: an array of PCs, ids, and indexes for the remaining proteins
		that successfully pass the user defined filter (or the default set to 0)
		"""
		# Extract all proteins with non-zero protein counts in both variants:
		nonzero_PC_idxs = np.nonzero(nonfiltered_protein_counts)
		nonzero_PCs = nonfiltered_protein_counts[nonzero_PC_idxs]
		nonzero_PC_idxs = nonzero_PC_idxs[0]
		#todo why are there nonzero_PC_ids and nonzero_ids? do they make the same output?
		nonzero_PC_ids = self.get_ids(nonfiltered_monomer_idx_dict,
									   nonzero_PC_idxs)

		if filter_num == 0:
			pass
		else:
			filter_PC_idxs = np.nonzero(nonzero_PCs > filter_num)
			nonzero_PCs = nonzero_PCs[filter_PC_idxs]
			nonzero_PC_ids = np.array(nonzero_PC_ids)
			nonzero_PC_ids = nonzero_PC_ids[filter_PC_idxs]
			nonzero_PC_idxs = filter_PC_idxs

		return nonzero_PCs, nonzero_PC_ids, nonzero_PC_idxs

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		"""
		Call data generating functions for saving data to csv files
		"""

		# Create saving paths
		#todo Figure out how to make this be placed in just the "out" folder!
		# todo: note that right now I think I am not saving the data and averaging over seeds, just saving the data per seed (and averaging over all those gens. ask for a second opinion on whether this is what I want)

		outDir = plotOutDir[:-8]
		save_pth = os.path.join(outDir, "multigen_saved_protein_count_data")
		pth = os.path.join(save_pth, "unfiltered_data")
		log_pth = os.path.join(pth, "log_data")

		if not os.path.exists(save_pth):
			os.mkdir(save_pth)
			os.mkdir(pth)
			os.mkdir(log_pth)

		# define/initialize commonly used variables
		self.n_total_gens = self.ap.n_generation
		seeds = self.ap.get_seeds()
		self.all_monomer_ids = []
		total_protein_counts = self.generate_data(simDataFile, validationDataFile)
		monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
									  enumerate(self.all_monomer_ids)}

		# Save unfiltered data
		startGen = IGNORE_FIRST_N_GENS + 1 # account for the fact that python numbering starts at 0
		col_labels = ["Monomer ID", "Average Protein Count"]
		ids = np.transpose(np.array(self.all_monomer_ids)); ids = ids.reshape(-1, 1)
		PCs_current = np.transpose(np.array(total_protein_counts))
		PCs_current = PCs_current.reshape(-1, 1)
		values = np.concatenate((ids, PCs_current), axis=1)
		save_file(
			pth,
			f'AvgProteinCounts_startGen_'
			f'{startGen}.csv',
			col_labels, values)

		# Save log data
		log_PCs = self.get_LogData(self.all_monomer_ids, total_protein_counts)
		log_col_labels = ["Monomer ID", "Log10 Average Protein Count"]
		log_PC_ids = np.transpose(np.array(self.all_monomer_ids)); log_PC_ids = log_PC_ids.reshape(-1, 1)
		log_PCs = np.transpose(np.array(log_PCs)); log_PCs = log_PCs.reshape(-1, 1)
		log_values = np.concatenate((log_PC_ids, log_PCs), axis=1)
		save_file(
			log_pth,
			f'LogAvgProteinCounts_startGen_{startGen}.csv', log_col_labels, log_values)

		# Create and save filtered data
		if save_filtered_data == 1:
			F_pth = os.path.join(save_pth, "filtered_data")
			F_log_pth = os.path.join(F_pth, "log_data")
			if not os.path.exists(F_pth):
				os.mkdir(F_pth)
				os.mkdir(F_log_pth)
			# Manditory Data Filtration
			F_PCs, F_PC_ids, F_PC_idxs = (
				self.filter_data(total_protein_counts, monomer_idx_dict_PreFilter,
								filter_num))
			Fcol_labels = ["Filtered Monomer ID", "Average Protein Count"]
			F_PC_ids = np.transpose(np.array(F_PC_ids)); F_PC_ids = F_PC_ids.reshape(-1, 1)
			F_PCs = np.transpose(np.array(F_PCs)); F_PCs = F_PCs.reshape(-1, 1)
			Fvalues = np.concatenate((F_PC_ids, F_PCs), axis=1)
			save_file(
				F_pth,
				f'Filtered_AvgProteinCounts_startGen_{startGen}.csv', Fcol_labels, Fvalues)

			# Save log data
			# todo: double check that entry #2 in the following line should be total_protein_counts
			F_log_PCs = self.get_LogData(F_PC_idxs, total_protein_counts, F_PC_idxs)
			F_log_col_labels = ["Filtered Monomer ID", "Log10 Average Protein Count"]
			F_log_PC_ids = np.transpose(np.array(F_PC_ids)); F_log_PC_ids = F_log_PC_ids.reshape(-1, 1)
			F_log_PCs = np.transpose(np.array(F_log_PCs)); F_log_PCs = F_log_PCs.reshape(-1, 1)
			F_log_values = np.concatenate((F_log_PC_ids, F_log_PCs), axis=1)
			save_file(
				F_log_pth,
				f'Filtered_LogAvgProteinCounts_startGen_{startGen}.csv', F_log_col_labels, F_log_values)

		# Save protein counts with available validation data:
		if save_with_validation_data == 1:
			validation_pth = os.path.join(save_pth, "validation_data")
			validation_log_path = os.path.join(validation_pth, "log_data")
			if not os.path.exists(validation_pth):
				os.mkdir(validation_pth)
				os.mkdir(validation_log_path)
			# Get validation data
			sim_schmidt_counts, val_schmidt_counts, s_ids, sim_wisniewski_counts, val_wisniewski_counts, w_ids = self.get_validation_data(
				inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata)

			# Save Wisniewski data:
			col_labels = ["Monomer ID", "Simulated Wisniewski Counts", "Validation Wisniewski Counts"]
			ids = np.transpose(np.array(w_ids));
			ids = ids.reshape(-1, 1)
			sim_wisniewski_counts = np.transpose(np.array(sim_wisniewski_counts));
			sim_wisniewski_counts = sim_wisniewski_counts.reshape(-1, 1)
			val_wisniewski_counts = np.transpose(np.array(val_wisniewski_counts));
			val_wisniewski_counts = val_wisniewski_counts.reshape(-1, 1)
			values = np.concatenate((ids, sim_wisniewski_counts, val_wisniewski_counts), axis=1)
			save_file(
				validation_pth,
				f'Wisniewski_Comparison_startGen_{IGNORE_FIRST_N_GENS + 1}.csv', col_labels, values)

			# Save Schmidt data:
			col_labels = ["Monomer ID", "Simulated Schmidt Counts", "Validation Schmidt Counts"]
			ids = np.transpose(np.array(s_ids));
			ids = ids.reshape(-1, 1)
			sim_schmidt_counts = np.transpose(np.array(sim_schmidt_counts));
			sim_schmidt_counts = sim_schmidt_counts.reshape(-1, 1)
			val_schmidt_counts = np.transpose(np.array(val_schmidt_counts));
			val_schmidt_counts = val_schmidt_counts.reshape(-1, 1)
			values = np.concatenate((ids, sim_schmidt_counts, val_schmidt_counts), axis=1)
			save_file(
				validation_pth,
				f'Schmidt_Comparison_startGen_{IGNORE_FIRST_N_GENS + 1}.csv', col_labels, values)

			# Save log data for Wisniewski validation data:
			log_sim_wisniewski_counts = self.get_LogData(sim_wisniewski_counts, sim_wisniewski_counts)
			log_val_wisniewski_counts = self.get_LogData(val_wisniewski_counts, val_wisniewski_counts)
			log_col_labels = ["Monomer ID", "Log10 Simulated Wisniewski Counts", "Log10 Validation Wisniewski Counts"]
			log_ids = np.transpose(np.array(w_ids));
			log_ids = log_ids.reshape(-1, 1)
			log_sim_wisniewski_counts = np.transpose(np.array(log_sim_wisniewski_counts));
			log_sim_wisniewski_counts = log_sim_wisniewski_counts.reshape(-1, 1)
			log_val_wisniewski_counts = np.transpose(np.array(log_val_wisniewski_counts));
			log_val_wisniewski_counts = log_val_wisniewski_counts.reshape(-1, 1)
			log_values = np.concatenate((log_ids, log_sim_wisniewski_counts, log_val_wisniewski_counts), axis=1)
			save_file(
				validation_log_path,
				f'Log10_Wisniewski_Comparison_startGen_{IGNORE_FIRST_N_GENS + 1}.csv', log_col_labels, log_values)

			# Save log data for Schmidt validation data:
			log_sim_schmidt_counts = self.get_LogData(sim_schmidt_counts, sim_schmidt_counts)
			log_val_schmidt_counts = self.get_LogData(val_schmidt_counts, val_schmidt_counts)
			log_col_labels = ["Monomer ID", "Log10 Simulated Schmidt Counts", "Log10 Validation Schmidt Counts"]
			log_ids = np.transpose(np.array(s_ids));
			log_ids = log_ids.reshape(-1, 1)
			log_sim_schmidt_counts = np.transpose(np.array(log_sim_schmidt_counts));
			log_sim_schmidt_counts = log_sim_schmidt_counts.reshape(-1, 1)
			log_val_schmidt_counts = np.transpose(np.array(log_val_schmidt_counts));
			log_val_schmidt_counts = log_val_schmidt_counts.reshape(-1, 1)
			log_values = np.concatenate((log_ids, log_sim_schmidt_counts, log_val_schmidt_counts), axis=1)
			save_file(
				validation_log_path,
				f'Log10_Schmidt_Comparison_startGen_{IGNORE_FIRST_N_GENS + 1}.csv', log_col_labels, log_values)

if __name__ == "__main__":
	Plot().cli()

