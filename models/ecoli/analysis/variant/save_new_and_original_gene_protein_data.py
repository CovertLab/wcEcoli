"""
Anaylsis Script for saving all protein counts impacted by new genes. Note that
this file only works for saving data for simulations that contain at least two
variants (one control and one experimental variant).
"""

import pickle
import os
import csv
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

IGNORE_FIRST_N_GENS = 1

"""
Enter a "filter number" below to filter out any proteins that have counts less
than a given filter number. This is useful for removing proteins that have very
low counts in a variant(s) that could likely be due to sub-generational gene expression. 
The default is set to 0, but can be changed to any number greater than 0.
"""
# Number to be set as the minimum threshold PC value proteins must have in both
# variants in order to be used in the plots to follow (set to 0 for default):
filter_num = 1

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
		Generates data variables that will be used for subsequent graphs
		Args:
			simDataFile: simulation data file
		Returns:
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
			all_cells = self.ap.get_cells(variant=[variant],
										  generation=np.arange(IGNORE_FIRST_N_GENS,
															   self.n_total_gens),
										  only_successful=True)
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
		return protein_counts

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
		Covert normal protein count data to their log values
		Args:
			protein_idxs: an array of the indices for proteins of interest
			interest_protein_counts: the full data structure  of protein counts
			(usually size variants by # of proteins), either filtered or unfiltered
			index_vals: if the protein idxs are not in sequential order
		Returns: an data structure of the log version of interest PCs
		"""
		avg_log_interest_proteins = np.zeros((
			len(self.variant_pair), len(protein_idxs)))
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
		Filter the data to extract all proteins with 0 PCs in at least variant
		Args:
			nonfiltered_protein_counts: array of PCs for all proteins
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
		nonzero_PCs_ids = self.get_ids(nonfiltered_monomer_idx_dict,
									   shared_nonzero_PCs_idxs)
		nonzero_ids = self.all_monomer_ids[shared_nonzero_PCs_idxs]

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
			nonzero_ids = shared_filtered_PC_idxs
		return nonzero_PCs, nonzero_PCs_ids, nonzero_ids

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		"""
		Call data generating functions for saving
		"""

		# Create saving paths
		outDir = plotOutDir[:-8]
		save_pth = os.path.join(outDir, "saved_data")
		pth = os.path.join(save_pth, "unfiltered_data")
		F_pth = os.path.join(save_pth, "filtered_data")

		if not os.path.exists(save_pth):
			os.mkdir(save_pth)
			os.mkdir(pth)
			os.mkdir(F_pth)

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
			protein_counts = self.generate_data(simDataFile)
			monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
										  enumerate(self.all_monomer_ids)}

			# Manditory Data Filtration
			F_PCs, F_PC_ids, F_PC_idxs = (
				self.filter_data(protein_counts, monomer_idx_dict_PreFilter,
								 filter_num))

			# Save unfiltered data
			expstr = "var_" + str(experimental_var) + "_avg_PCs"
			col_labels = ["all_monomer_ids","var_0_avg_PCs", expstr]
			ids = [np.array(self.all_monomer_ids)]
			ids = np.transpose(ids)
			PCs_current = self.total_protein_counts
			values = np.concatenate((ids, PCs_current.T), axis=1)
			save_file(
				pth,
				f'wcm_full_monomers_var_{experimental_var}_startGen_'
				f'{IGNORE_FIRST_N_GENS}.csv',
				col_labels, values)

			# Save filtered data
			Fcol_labels = ["filtered_monomer_ids", "var_0_avg_PCs", expstr]
			F_PC_ids = [np.array(F_PC_ids)]; F_PC_ids = np.transpose(F_PC_ids)
			F_PCs_current = F_PCs
			Fvalues = np.concatenate((F_PC_ids, F_PCs.T), axis=1)
			save_file(
				F_pth,
				f'wcm_filter_{filter_num}_monomers_var_{experimental_var}'
				f'_startGen_{IGNORE_FIRST_N_GENS}.csv',
				Fcol_labels, Fvalues)

if __name__ == "__main__":
	Plot().cli()
