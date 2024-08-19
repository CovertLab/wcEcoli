"""
Anaylsis Script for saving all protein counts impacted by new genes. Note that
this file only works for saving data for simulations that contain at least two
variants (one control and one experimental variant). If you are looking to save
individual variants, using the #todo write the cohort script name here
"""

import pickle
import os
import csv
from typing import cast, List, Tuple
from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.io.tablereader import TableReader

"""
Indicate which generation the data should start being collected from (sometimes 
this number should be greater than 0 because the first few generations may not 
be representative of the true dynamics occuring in the cell).
"""

IGNORE_FIRST_N_GENS = 14

"""
Enter a "filter number" below to filter out any proteins that have counts less
than a given filter number. This is useful for removing proteins that have very
low counts in a variant(s) that could likely be due to sub-generational gene expression. 
The default is set to 0, but can be changed to any number greater than 0.
"""
# Set to True if you want to save the filtered data:
save_filtered_data = True

# Number to be set as the minimum threshold PC value proteins must have in both
# variants in order to be used in the plots to follow (set to 0 for default):
filter_num = 0

"""
Indicate whether the proteins with corresponding validation data should also 
be saved here. True will save the validation data, False will not.
"""
save_validation_data = True

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

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def generate_data(self, simDataFile):
		"""
		Generates data that will be used for saving protein count data
		Args:
			simDataFile: simulation data file
		Returns:
			protein_counts: protein count (PC) data for all proteins (originally
			present on the E.coli chromosome) in the simulation for each variant
			(the PC for each protein is averaged over all seeds & generations)

			self.total_protein_counts: the original E. coli chromosome proteins'
			PCs and new genes' (NG) PCs

			self.new_gene_monomer_ids: protein ids for new genes inserted into
			the E.coli genome

			self.original_gene_ids: protein ids for the original proteins on
			the E.coli genome

			self.all_monomer_ids: list of all the monomer ids (NG protein ids
			and orginal proteins' gene ids)
		"""
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = (
			sim_data.process.transcription.cistron_data.struct_array)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)
		new_gene_mRNA_ids = mRNA_sim_data[
			mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		self.new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in new_gene_mRNA_ids]

		self.all_monomer_ids = monomer_sim_data['id']
		self.original_monomer_ids = np.delete(self.all_monomer_ids, np.where(
			self.all_monomer_ids == self.new_gene_monomer_ids))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(self.all_monomer_ids)}
		protein_counts = np.zeros((len(self.variant_pair),
								   len(self.original_monomer_ids)))
		self.total_protein_counts = np.zeros((len(self.variant_pair),
											  len(self.all_monomer_ids)))
		for var_idx in range(len(self.variant_pair)):
			variant = self.variant_pair[var_idx]
			all_cells = self.ap.get_cells(
				variant=[variant],generation=np.arange(IGNORE_FIRST_N_GENS,
				self.n_total_gens), only_successful=True)

			# Get the protein counts for each gene/protein
			gen_avg_monomer_counts = (
				read_stacked_columns(all_cells,'MonomerCounts',
									 'monomerCounts',
									 fun=lambda x: np.mean(x[:], axis=0)))
			total_avg_gene_counts = np.mean(gen_avg_monomer_counts, axis=0)
			self.total_protein_counts[var_idx] = total_avg_gene_counts
			old_gene_idxs = [monomer_idx_dict.get(monomer_id)
							 for monomer_id in self.original_monomer_ids]
			avg_gene_monomer_counts = total_avg_gene_counts[old_gene_idxs]
			protein_counts[var_idx] = avg_gene_monomer_counts

		protein_counts = np.array(protein_counts)
		self.total_protein_counts = np.array(self.total_protein_counts)

		return

	def get_simulated_validation_counts(self,
			validation_counts, monomer_counts, validation_ids, simulation_ids):
		# type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray) -> Tuple[np.ndarray, np.ndarray]
		"""
		ADAPTED FROM wholecell/utils/protein_counts.py

		Get simulated counts and validation counts of monomers that exist in both
		the simulation and the validation data

		Arguments:
			validation_counts: Monomer counts from validation data.
			monomer_counts: Simulated monomer counts (from translation
				process).
			validation_ids: Monomer IDs from validation data. IDs
				must appear in same order as in validation_counts.
			simulation_ids: IDs of monomers in the same order as
				monomer_counts.

		Returns:
			The simulated counts of the monomers that appear in the
			validation data, and the validation counts of the monomers in the same
			order.
		"""

		avg_sim_counts = monomer_counts

		sim_ids_lst = cast(List[str], simulation_ids.tolist())
		val_ids_lst = cast(List[str], validation_ids.tolist())
		overlapping_ids_set = set(sim_ids_lst) & set(val_ids_lst)

		sim_id_to_index_map = {
			sim_id: i for i, sim_id in enumerate(sim_ids_lst)
			if sim_id in overlapping_ids_set}
		val_id_to_index_map = {
			val_id: i for i, val_id in enumerate(val_ids_lst)
			if val_id in overlapping_ids_set}

		overlapping_ids_list = list(overlapping_ids_set)
		sim_filtered_idx = np.array([
			sim_id_to_index_map[monomer_id] for monomer_id in overlapping_ids_list
		])
		val_filtered_idx = np.array([
			val_id_to_index_map[monomer_id] for monomer_id in overlapping_ids_list
		])
		hello = 2

		return avg_sim_counts[:, sim_filtered_idx], validation_counts[
			val_filtered_idx], overlapping_ids_list

	def get_validation_data(self, simDataFile, validationDataFile):
		# todo: add docstring
		# adapted from multigen and single/proteinCountsValidation.py
		sim_data = self.read_pickle_file(simDataFile)
		validation_data = self.read_pickle_file(validationDataFile)

		# Aquire the monomer ids and counts for simulation and validation data:
		monomer_ids = sim_data.process.translation.monomer_data["id"]
		wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
		schmidt_ids = validation_data.protein.schmidt2015Data["monomerId"]
		wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]
		monomer_counts = self.total_protein_counts

		# Obtain the protein counts for protiens that are present in both the simluation and validation datasets:
		sim_schmidt_counts, val_schmidt_counts, schmidt_overlap_ids = self.get_simulated_validation_counts(
			schmidt_counts, monomer_counts, schmidt_ids, monomer_ids)
		sim_wisniewski_counts, val_wisniewski_counts, wisniewski_overlap_ids = self.get_simulated_validation_counts(
			wisniewski_counts, monomer_counts, wisniewski_ids, monomer_ids)


		return sim_schmidt_counts, val_schmidt_counts, schmidt_overlap_ids, sim_wisniewski_counts, val_wisniewski_counts, wisniewski_overlap_ids

	def get_ids(self, monomer_idx_dict, protein_idxs):
		"""
		Obtain the protein ids for proteins based on their indexes
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
		# todo edit this docstring
		Covert normal protein count data to their log values
		Args:
			protein_idxs: an array of the indices for proteins of interest

			interest_protein_counts: the full data structure  of protein counts
			(usually size variants by # of proteins), either filtered or unfiltered

			index_vals: if the protein idxs are not in sequential order

		Returns: an data structure of the log version of interest PCs
		"""
		hi = 5
		avg_log_interest_proteins = np.zeros((
			len(self.variant_pair), len(protein_idxs)))
		for variant in range(len(interest_protein_counts)):
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
		Filter the data to extract all proteins with 0 PCs in at least variant
		Args:
			nonfiltered_protein_counts: array of PCs for all proteins # todo: should this always include the new gene as well????
			nonfiltered_monomer_idx_dict: dictionary that maps id to index
			filter_num: the minimum number of PCs a protein must have in both
			variants to avoid being discarded from the new filtered data
			(by default this is set to zero).
		Returns: an array of PCs, ids, and indexes for the remaining proteins
		that successfully pass the user defined filter (or the default set to 0)
		"""
		# Extract all proteins with non-zero protein counts in both variants:
		nonzero_p_counts_var0_idxs = np.nonzero(nonfiltered_protein_counts[0])
		nonzero_p_counts_var1_idxs = np.nonzero(nonfiltered_protein_counts[1])
		shared_nonzero_PCs_idxs = np.intersect1d(nonzero_p_counts_var0_idxs,
												 nonzero_p_counts_var1_idxs)
		nonzero_PCs = nonfiltered_protein_counts[:, shared_nonzero_PCs_idxs]
		#nonzero_PCs_ids = self.get_ids(nonfiltered_monomer_idx_dict,
									   #shared_nonzero_PCs_idxs)
		# todo: do I need this line:
		nonzero_PCs_ids = self.all_monomer_ids[shared_nonzero_PCs_idxs]
		nonzero_idxs = shared_nonzero_PCs_idxs

		if filter_num == 0:
			pass
		else:
			var0_filter_PCs_idxs = np.nonzero(nonzero_PCs[0] > filter_num)
			var1_filter_PCs_idxs = np.nonzero(nonzero_PCs[1] > filter_num)
			shared_filtered_PC_idxs = np.intersect1d(var0_filter_PCs_idxs,
													 var1_filter_PCs_idxs)
			nonzero_PCs = nonzero_PCs[:, shared_filtered_PC_idxs]
			nonzero_PCs_ids = np.array(nonzero_PCs_ids)
			nonzero_PCs_ids = nonzero_PCs_ids[shared_filtered_PC_idxs]
			nonzero_idxs = shared_filtered_PC_idxs
		return nonzero_PCs, nonzero_PCs_ids, nonzero_idxs

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		"""
		# todo: add docstring
		Call data generating functions for saving
		"""

		# Create saving paths
		save_pth = os.path.join(inputDir, "saved_data")
		pth = os.path.join(save_pth, "unfiltered_data")
		log_pth = os.path.join(pth, "log_data")

		F_pth = os.path.join(save_pth, "filtered_data")

		if not os.path.exists(save_pth):
			os.mkdir(save_pth)
			os.mkdir(pth)
			os.mkdir(log_pth)

		# Get data for all variants
		all_variants = self.ap.get_variants()
		control_var = all_variants[0]
		experimental_vars = all_variants[1:]
		for variant_2 in experimental_vars:
			self.variant_pair = [control_var, variant_2]
			experimental_var = self.variant_pair[1]

			# define/initialize commonly used variables
			self.n_total_gens = self.ap.n_generation
			self.all_monomer_ids = []
			self.original_monomer_ids = []
			self.new_gene_monomer_ids = []
			self.total_protein_counts = []

			# Generate data
			self.generate_data(simDataFile)
			monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
										  enumerate(self.all_monomer_ids)}

			# Save unfiltered data
			expstr = "var_" + str(experimental_var) + "_avg_PCs"
			col_labels = ["all_monomer_ids", "var_0_avg_PCs", expstr]
			ids = [np.array(self.all_monomer_ids)]
			ids = np.transpose(ids)
			PCs_current = self.total_protein_counts
			values = np.concatenate((ids, PCs_current.T), axis=1)
			save_file(
				pth,
				f'AvgProteinCounts_Variant_{experimental_var}_startGen_'
				f'{IGNORE_FIRST_N_GENS}.csv',
				col_labels, values)

			# Save unfiltered log data
			log_C_PCs = (self.get_LogData(self.all_monomer_ids, self.total_protein_counts))[0]
			log_E_PCs = (self.get_LogData(self.all_monomer_ids, self.total_protein_counts))[1]
			log_col_labels = ["Monomer ID", "Control Variant Log10 Average Protein Count", "Experimental Variant Log10 Average Protein Count"]
			log_PC_ids = self.all_monomer_ids
			log_values = np.vstack((log_PC_ids, log_C_PCs, log_E_PCs))
			log_values = np.transpose(log_values)

			save_file(
				log_pth,
				f'LogAvgProteinCounts_Variant_{experimental_var}_startGen_'
				f'{IGNORE_FIRST_N_GENS}.csv',
				log_col_labels, log_values)

			# Save filtered data
			if save_filtered_data == True:

				# Manditory Data Filtration
				# todo: should protein_counts be replaced by self.total_protein_counts? bc monomer_idx_dict_PreFilter is based on self.all_monomer_ids!!!
				F_PCs, F_PC_ids, F_PC_idxs = \
					(self.filter_data(self.total_protein_counts, monomer_idx_dict_PreFilter,
									 filter_num))


				F_pth = os.path.join(save_pth, "filtered_data")
				F_log_pth = os.path.join(F_pth, "log_data")
				if not os.path.exists(F_pth):
					os.mkdir(F_pth)
					os.mkdir(F_log_pth)

				# save filtered data:
				expstr = "Experimental Variant" + str(experimental_var) + " Log10 Average Protein Count"
				Fcol_labels = ["Filtered Monomer ID", "Control Varaint Average Protein Count", expstr]
				F_PC_ids = [np.array(F_PC_ids)];
				F_PC_ids = np.transpose(F_PC_ids)
				Fvalues = np.concatenate((F_PC_ids, F_PCs.T), axis=1)
				save_file(
					F_pth,
					f'AvgProteinCounts_Variant_{experimental_var}_'
					f'filter_{filter_num}_startGen_{IGNORE_FIRST_N_GENS}.csv', Fcol_labels, Fvalues)

				# save filtered log data:
				F_log_C_PCs = (self.get_LogData(F_PC_idxs, self.total_protein_counts, F_PC_idxs))[0]
				F_log_E_PCs = (self.get_LogData(F_PC_idxs, self.total_protein_counts, F_PC_idxs))[1]
				F_log_col_labels = ["Filtered Monomer ID", "Control Varaint Log10 Average Protein Count", "Experimental Variant Log10 Average Protein Count"]
				F_log_PC_ids = np.squeeze(F_PC_ids)
				F_log_values = np.vstack((F_log_PC_ids, F_log_C_PCs, F_log_E_PCs))
				F_log_values = np.transpose(F_log_values)
				save_file(
					F_log_pth,
					f'LogAvgProteinCounts_Variant_{experimental_var}_'
					f'filter_{filter_num}_startGen_{IGNORE_FIRST_N_GENS}.csv', F_log_col_labels,
					F_log_values)

			# Save validation data
			if save_validation_data:
				validation_pth = os.path.join(save_pth, "validation_data")
				validation_log_path = os.path.join(validation_pth, "log_data")
				schmidt_pth = os.path.join(validation_pth, "schmidt_data")
				wisniewski_pth = os.path.join(validation_pth, "wisniewski_data")
				schmidt_log_pth = os.path.join(validation_log_path, "schmidt_data")
				wisniewski_log_pth = os.path.join(validation_log_path, "wisniewski_data")
				if not os.path.exists(validation_pth):
					os.mkdir(validation_pth)
					os.mkdir(validation_log_path)
					os.mkdir(schmidt_pth)
					os.mkdir(wisniewski_pth)
					os.mkdir(schmidt_log_pth)
					os.mkdir(wisniewski_log_pth)

				# Get validation data for the control and experimental variants:
				(sim_schmidt_counts, val_schmidt_counts, s_ids, sim_wisniewski_counts,
				 val_wisniewski_counts, w_ids) = self.get_validation_data(
					simDataFile, validationDataFile)

				# Save Schmidt data:
				expstr = "Experimental Variant " + str(experimental_var) + " Schmidt Counts"
				col_labels = ["Monomer ID", "Validation Schmidt Counts",
							  "Control Variant Schmidt Counts", expstr]

				values = np.vstack((s_ids, val_schmidt_counts, sim_schmidt_counts))
				values = np.transpose(values)
				save_file(
					schmidt_pth,
					f'SchmidtValidation_Variant_{experimental_var}'
					f'_startGen_{IGNORE_FIRST_N_GENS}.csv',
					col_labels, values)

				# Save log data for Schmidt validation data:
				# todo: should I be making these into functions?t
				log_sim_schmidt_counts = self.get_LogData(s_ids,
															 sim_schmidt_counts)
				log_val_schmidt_counts = self.get_LogData(s_ids,
															 [val_schmidt_counts,
															  val_schmidt_counts])[0]

				expstr = "Log10 Experimental Variant " + str(
					experimental_var) + " Schmidt Counts"
				log_col_labels = ["Monomer ID", "Log10 Validation Schmidt Counts",
								  "Log10 Control Variant Schmidt Counts",
								  expstr]

				log_s_values = np.vstack(
					(s_ids, log_val_schmidt_counts, log_sim_schmidt_counts))
				log_s_values = np.transpose(log_s_values)
				save_file(
					schmidt_log_pth,
					f'Log10_SchmidtValidation_Variant_{experimental_var}'
					f'_startGen_{IGNORE_FIRST_N_GENS}.csv',
					log_col_labels, log_s_values)


				# Save Wisniewski data:
				expstr = "Experimental Variant " + str(experimental_var) + " Wisniewski Counts"
				col_labels = ["Monomer ID", "Validation Wisniewski Counts",
							  "Control Variant Wisniewski Counts", expstr]

				values = np.vstack((w_ids, val_wisniewski_counts, sim_wisniewski_counts))
				values = np.transpose(values)
				save_file(
					wisniewski_pth,
					f'WisniewskiValidation_Variant_{experimental_var}'
					f'_startGen_{IGNORE_FIRST_N_GENS}', col_labels,
					values)

				# Save log data for Wisniewski validation data:
				# todo: should I be making these into functions?
				log_sim_wisniewski_counts = self.get_LogData(w_ids,
															 sim_wisniewski_counts)
				log_val_wisniewski_counts = self.get_LogData(w_ids,
															 [val_wisniewski_counts,val_wisniewski_counts])[0]

				expstr = "Log10 Experimental Variant " + str(experimental_var) + " Wisniewski Counts"
				log_col_labels = ["Monomer ID", "Log10 Validation Wisniewski Counts",
								  "Log10 Control Variant Wisniewski Counts",
								  expstr]

				log_w_values = np.vstack((w_ids, log_val_wisniewski_counts, log_sim_wisniewski_counts))
				log_w_values = np.transpose(log_w_values)
				save_file(
					wisniewski_log_pth,
					f'Log10_WisniewskiValidation_Variant_{experimental_var}'
					f'_startGen_{IGNORE_FIRST_N_GENS}.csv',
					log_col_labels, log_w_values)


if __name__ == "__main__":
	Plot().cli()
