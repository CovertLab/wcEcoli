"""
Anaylsis Script for analyzing all protein counts impacted by new genes.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
											   read_bulk_molecule_counts, read_stacked_bulk_molecules,
											   read_stacked_columns, stacked_cell_threshold_mask)
from wholecell.io.tablereader import TableReader

"""
1 to plot early (before MIN_LATE_CELL_INDEX), and late generations in
addition to all generations
"""
exclude_early_gens = 1

FONT_SIZE=9
MAX_VARIANT = 10 # do not include any variant >= this index

"""
generations before this may not be representative of dynamics 
due to how they are initialized
"""
IGNORE_FIRST_N_GENS = 4


"""
Indicate here the number of proteins that should be plotted for comparison in 
comparison plots 1 and 3  (ideally a value between 5 and 100 should plot well,
but technically any value between 2 and 4308 should work for graphs 1 and 3). 
Setting one of the nubmers corresponding to graphs 1-4 equal to zero will prevent
that graph from being generated. The user may need to manually edit the graphing 
specification code if the number given is not close to the range the graph was 
designed for originally (see more about this in the descriptions for each graph
within the code).  
"""

"""
The first four graph types do not need any filter applied in order to be generated
(and will include any proteins with zero PCs in the resultant graph).
"""

'''
Graph 1: Random protein count (PC) comparisons
'''
# Number of random proteins to compare protein counts (PCs) of between variants:
randnum = 11
'''
Graph 2a: PC comparisons between two variants
Graph 2b: Graph the data from Graph 2a using a log scale
'''
# If the NG is to be included on the plots, set this variable equal to 1:
include_NG_G2 = 1
# a: Set this value to 1 to create this comparison graph, 0 otherwise:
var_PC_comparison = 1
# b: Set this value to 1 to create Graph 2a on a log scale as well:
var_PC_comparison_LogScale = 1
'''
Graph 3: Largest absolute PC differences 
'''
# Number of proteins with the largest PC increase between variants
# to observe (this one will include the new gene in the graph):
num = 12
'''
Graph 4: Visualize proteins that appear uniquely in only one variant
'''
# Set this value to 1 to generate this graph, 0 otherwise:
unique_PC_appearances = 1

"""
The following graphs require a manditory filter to be applied to the data first
in order to be generated. This filter by default is set to 0 (where in any 
proteins with a PC value of 0 in at least one variant will be filtered out). 
A number greater than zero can also be specified below in "filter_num" to 
filter out proteins that have PCs below a desired threshold. 

Note: New genes (that appear in the experimental variant but the the control)
will be completely filtered out by the manditory filtration as it will have 0
PCs in the control variable, so it will not show up in any of these graphs.
"""
# Number to be set as the minimum threshold PC value proteins must have in both
# variants in order to be used in the plots to follow (set to 0 by default):
filter_num = 1

'''
Graph 5a: Filtered PC comparisons between two variants
Graph 5b: Graph the data from Graph 5a on a log scale
'''
# a: Set this value to 1 to create this comparison graph, 0 otherwise:
var_PC_comparison_wF = 1
# b: Set this value to 1 to create Graph 5a on a log scale as well:
var_PC_comparison_wF_LogScale = 1
''' 
Graph 6a: Plot proteins with the smallest difference in PCs between variants
Graph 6b: Same as Graph 6a but the change in PCs between vars is included 
'''
# Number of proteins with the smallest PC change between variants:
min_num = 13
# set to 1 to create a graph that displays the % differnce on the graph:
show_PC_diff_6b = 1
'''
Graph 7a: Plot proteins with the greatest difference in PCs
Graph 7b: Same as Graph 7a but the change in PCs between vars is included 
'''
# Number of proteins with the greatest change in PCs between variants:
max_num = 14
# set this to 1 to create a graph that displays the % differnce on the graph:
show_PC_diff_7b = 1
'''
Graph 8a: Proteins with the max fold increase between variants
Graph 8b: Same as Graph 8a but NOT plotted on a log scale 
'''
# Number of proteins to observe with the greatest fold increase in PCs:
max_fold_num = 15
# set this to 1 to create the same graph plotted without a log scale:
max_fold_num_woLogScale = 0
''' 
Graph 9a: Proteins with the max fold decrease between variants
Graph 9b: Same as Graph 9a but NOT plotted on a log scale 
'''
# Number of proteins to observe with the largest fold decrease PCs:
min_fold_num = 6
# set this to 1 to create the same graph plotted without a log scale:
min_fold_num_woLogScale = 1
'''
Graph 10a: visualize the max increases and decreases in PCs side by side
Graph 10b: Graph 10a on a log scale
'''
# Number of proteins to compare max increases and decreases in PCs:
shared_diff_num = 9
# set this to 1 to plot the same comparison but with a log scale:
shared_diff_LogScale = 0
'''
Graph 11: visualize the max fold increases and decreases in PCs side by side
'''
# Number of proteins to observe the max fold increases and decreases in PCs:
sharednum = 9

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def generate_data(self, simDataFile):
		"""
		Generates data variables that will be used for subsequen graphs
		Args:
			simDataFile: simulation data file
		Returns:
			protein_counts: protein count (PC) data for all proteins (originally present on the E.coli chromosome) in the simulation for each variant (the PC for each protein is averaged over all the generations)
			self.total_protein_counts: the original PCs and new gene (NG) PCs in one variable
			self.new_gene_monomer_ids: protein ids for new genes inserted into the E.coli genome
			self.original_gene_ids: protein ids for the original proteins on the E.coli genome
			self.all_monomer_ids: list of all the monomer ids (NG protein ids and orginal proteins' gene ids)
		"""
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[
			mRNA_sim_data['is_new_gene']]['id'].tolist()
		mRNA_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										monomer_sim_data['id']))
		self.new_gene_monomer_ids = [mRNA_monomer_id_dict.get(mRNA_id)
								for mRNA_id in new_gene_mRNA_ids]
		if len(new_gene_mRNA_ids) == 0:
			print("This plot is intended to be run on simulations where the"
				  " new gene option was enabled, but no new gene mRNAs were "
				  "found.")
			return
		if len(self.new_gene_monomer_ids) == 0:
			print("This plot is intended to be run on simulations where the "
				  "new gene option was enabled, but no new gene proteins "
				  "were "
				  "found.")
			return
		assert len(self.new_gene_monomer_ids) == len(new_gene_mRNA_ids),\
			'number of new gene monomers and mRNAs should be equal'

		self.all_monomer_ids = monomer_sim_data['id']
		self.original_monomer_ids = np.delete(self.all_monomer_ids, np.where(
			self.all_monomer_ids == self.new_gene_monomer_ids))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(self.all_monomer_ids)}
		protein_counts = np.zeros((len(self.variants), len(self.original_monomer_ids)))
		self.total_protein_counts = np.zeros((len(self.variants), len(self.all_monomer_ids)))
		for variant in self.variants:
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
			self.total_protein_counts[variant] = total_avg_gene_counts
			old_gene_idxs = [monomer_idx_dict.get(monomer_id)
							 for monomer_id in self.original_monomer_ids]
			avg_gene_monomer_counts = total_avg_gene_counts[old_gene_idxs]
			protein_counts[variant] = avg_gene_monomer_counts

		protein_counts = np.array(protein_counts)
		self.total_protein_counts = np.array(self.total_protein_counts)
		return protein_counts

	def get_idxs(self, monomer_names):
		"""
		Obtain the indexes for specific proteins using their monomer ids
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
			len(self.variants), len(protein_idxs)))
		for variant in self.variants:
			for idx in range(len(protein_idxs)):
				if len(index_vals) == 0:
					index = idx
				else:
					index = index_vals[idx]
				avg_log_interest_proteins[variant][idx] = \
					np.log10(interest_protein_counts[variant][index] + 1)
		return avg_log_interest_proteins

	def find_unique_proteins(self, nonfiltered_protein_counts, nonfiltered_monomer_idx_dict):
		"""
		Obtain the proteins that have nonzero PCs in one variant only (these
		are most likely the proteins that experience subgenerational gene expression
		in the simulation).
		Args:
			nonfiltered_protein_counts: PC array of all original proteins
			nonfiltered_monomer_idx_dict: dictionary mapping protein names to
			their respective index in the array
		Returns: the PCs for proteins with unique appearances in only one variant,
		as well as their id and index within the array of all proteins.
		"""
		# Extract and view the proteins that have nonzero counts for a variant:
		p_counts = np.array(nonfiltered_protein_counts)
		# Determine indexes for proteins with protein count values of zero:
		variant0_zeros_idxs = np.where(p_counts[0] == 0)[0]
		variant1_zeros_idxs = np.where(p_counts[1] == 0)[0]
		shared_zeros_idxs = np.intersect1d(variant0_zeros_idxs,
										   variant1_zeros_idxs)
		# Obtain indexes of proteins w/ 0 counts that are unique for variants:
		var0_unique_0_PC_idxs = [idx for idx in variant0_zeros_idxs
								 if idx not in shared_zeros_idxs]
		var1_unique_0_PC_idxs = [idx for idx in variant1_zeros_idxs
								 if idx not in shared_zeros_idxs]
		# Obtain the protein counts of ALL the unique-zero protein counts:
		nonshared_0_PC_idxs = np.append(var0_unique_0_PC_idxs,
										var1_unique_0_PC_idxs, axis=0)
		nonshared_0_PC_ids = self.get_ids(nonfiltered_monomer_idx_dict, nonshared_0_PC_idxs)
		nonshared_0_PCs = p_counts[:, nonshared_0_PC_idxs]
		return nonshared_0_PCs, nonshared_0_PC_ids, nonshared_0_PC_idxs

	def filter_data(self, nonfiltered_protein_counts, nonfiltered_monomer_idx_dict, filter_num=0):
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
		nonzero_PCs_ids = self.get_ids(nonfiltered_monomer_idx_dict, shared_nonzero_PCs_idxs)
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

	def get_max_fold_increase(self, max_fold_num, filtered_PCs, filtered_ids, max_fold_num_woLogScale=0):
		"""
		Obtains the proteins with the max fold increase in PCs from the control
		variant (variant 0) and an experiemental variant (variant 1)
		Args:
			max_fold_num: # of proteins with max fold increases to observe
			filtered_PCs: array of filtered PCs (fails with non-filtered PCs)
			filtered_ids: corresponding ids for the filtered PCs
			max_fold_num_woLogScale: by default this output provides the PC data
			on a log scale, but setting this to 1 returns it without log scale
		Returns: the PCs, ids, and % difference in PCs from the experimental
		variant to the control variant (var 1 to var 0) for proteins with the
		greatest fold change between varients
		"""
		PC_max_folds = abs(filtered_PCs[1] / filtered_PCs[0])
		PC_max_fold_percents = (100 * (filtered_PCs[1] - filtered_PCs[0])
								/ filtered_PCs[0])
		psortdiff = np.argsort(PC_max_folds)
		PC_max_fold_idxs = psortdiff[-max_fold_num:]
		max_fold_changes = []
		max_fold_PCs_ids = []
		for idx in PC_max_fold_idxs:
			max_fold_changes.append(PC_max_fold_percents[idx])
			max_fold_PCs_ids.append(filtered_ids[idx])

		max_fold_PC_Log_values = self.get_LogData(PC_max_fold_idxs, filtered_PCs, PC_max_fold_idxs)
		max_fold_PC_values = max_fold_PC_Log_values

		if max_fold_num_woLogScale == 1:
			max_fold_PC_values = np.zeros((len(self.variants),
										   len(PC_max_fold_idxs)))
			for variant in self.variants:
				for idx in range(len(PC_max_fold_idxs)):
					index = PC_max_fold_idxs[idx]
					max_fold_PC_values[variant][idx] = filtered_PCs[variant][index]
		return max_fold_PC_values, max_fold_PCs_ids, max_fold_changes

	def get_max_fold_decrease(self, min_fold_num, filtered_PCs, filtered_ids, min_fold_num_woLogScale=0):
		"""
		Determines the proteins with the greatest fold decrease in PCs from
		the experimental variant to the control variant
		Args:
			min_fold_num: # of proteins with max fold decreases to observe
			filtered_PCs: array of filtered PCs (fails with non-filtered PCs)
			filtered_ids: corresponding ids for the filtered PCs
			min_fold_num_no_LogScale: by default this output provides the PC data
			on a log scale, but setting this to 1 returns it without log scale
		Returns: the PCs, ids, and % difference in PCs from the experimental
		variant to the control variant (var 1 to var 0) for proteins with the
		greatest fold decrease between variants
		"""
		PC_min_folds = abs(filtered_PCs[0] / filtered_PCs[1])
		PC_min_fold_percents = (100 * (filtered_PCs[1] - filtered_PCs[0])
								/ filtered_PCs[0])
		min_psortdiff = np.argsort(PC_min_folds)
		PC_min_fold_idxs = min_psortdiff[-min_fold_num:]
		min_fold_changes = []
		min_fold_PCs_ids = []
		for idx in PC_min_fold_idxs:
			min_fold_changes.append(PC_min_fold_percents[idx])
			min_fold_PCs_ids.append(filtered_ids[idx])

		min_fold_PC_Log_values = self.get_LogData(PC_min_fold_idxs, filtered_PCs, PC_min_fold_idxs)
		min_fold_PC_values = min_fold_PC_Log_values

		if min_fold_num_woLogScale == 1:
			min_fold_PC_values = np.zeros((len(self.variants),
										   len(PC_min_fold_idxs)))
			for variant in self.variants:
				for idx in range(len(PC_min_fold_idxs)):
					index = PC_min_fold_idxs[idx]
					min_fold_PC_values[variant][idx] = filtered_PCs[variant][index]
		return min_fold_PC_values, min_fold_PCs_ids, min_fold_changes

	def gen_G1(self, randnum):
		"""
		Plots the PCs for a randomly selected set of proteins
		Args:
			randnum: # of proteins to be randomly selected and plotted
			NOTE: the graph specifications (figsize=(10, 6)) are set to optimally
			plot for 10 randomly selected proteins (so this may need to be altered
			if larger numbers are used)
		Returns: a plot comparing the  the random proteins PCs' for each variant
		"""
		rand_monomers = np.random.choice(self.original_monomer_ids, size=randnum, replace=False)
		monomer_names = np.append(rand_monomers, self.new_gene_monomer_ids)
		protein_idxs = self.get_idxs(monomer_names)
		rand_PCs = self.total_protein_counts[:, protein_idxs]
		log_rand_PCs = self.get_LogData(protein_idxs, rand_PCs)

		# Make plot for the random protein selections:
		plt.figure(figsize=(10, 6))  # 8.5, 11 is good for 100 genes!

		plt.barh(monomer_names, log_rand_PCs[0], 0.1, align='edge',
				 label='no New Gene')
		plt.barh(monomer_names, log_rand_PCs[1], -0.1, align='edge', label='New Gene')

		plt.xlabel("log(Average Protein Count)", fontweight='bold')
		plt.ylabel("Protein ID", fontweight='bold')
		plt.legend()
		plt.title("Protein count comparisons between variants "
				  "\n for " + str(randnum)
				  + " randomly selected proteins ")
		plt.tight_layout()

	def gen_G2(self, protein_counts):
		"""
		Generates a plot of the control variant's PC data plotted against the
		experimental variants' PC data
		Args:
			protein_counts: the PCs for the desired proteins to be plotted
		Returns: an x-y style comparison plot
		"""
		var0_x = protein_counts[0]
		var1_y = protein_counts[1]

		plt.figure(figsize=(10, 10))
		plt.scatter(var0_x, var1_y, 1)
		m, b = np.polyfit(var0_x, var1_y, 1)
		plt.plot(var0_x, m * var0_x + b, linewidth=.5, color='#bcbd22')
		legstr = "linear fit: y = " + str(round(m, 2)) + "x + " + str(round(b, 2))
		max_pt = np.argsort(protein_counts)
		max_pt_idx = max_pt[0][-1]
		max_pt_x = var0_x[max_pt_idx]
		max_pt_y = var1_y[max_pt_idx]
		slope = max_pt_y / max_pt_x
		plt.plot(var0_x, slope * var0_x, linewidth=.5, color='#FFA500')
		otherstr = "y = " + str(round(slope, 2)) + "x"
		plt.legend(["PC data", legstr, otherstr])
		plt.xlabel("variant 0 (no New Gene)")
		plt.ylabel("variant 1 (New Gene)")
		plt.title(f"The {len(var0_x)} proteins plotted against each other")
		plt.tight_layout()

	def gen_G3(self, num, protein_counts, monomer_idx_dict):
		"""
		Find proteins with the largest absolute difference in PCs between variants
		Args:
			num: # of proteins to compare
			protein_counts: array of the all the protein count data
			monomer_idx_dict: protein id to index dictonary
		Returns: a plot of the proteins with the largest absolute difference
		in protein counts
		"""
		diff = abs(protein_counts[1] - protein_counts[0])
		sortdiff = np.argsort(diff)
		maxdiff = sortdiff[-num:]
		protein_ids = self.get_ids(monomer_idx_dict, maxdiff)
		protein_ids = np.append(protein_ids, self.new_gene_monomer_ids)
		protein_idxs = [monomer_idx_dict.get(monomer_id)
								 for monomer_id in protein_ids]
		interest_protein_counts = self.total_protein_counts[:, protein_idxs]
		PC_LogData = self.get_LogData(protein_idxs, interest_protein_counts)
		var0 = PC_LogData[0]
		var1 = PC_LogData[1]
		# Make plot for genes with the largest protein count differences:
		plt.figure(figsize=(10, 6))
		plt.barh(protein_ids, var0, 0.1, align='edge',
				 label='no New Gene')
		plt.barh(protein_ids, var1, -0.1, align='edge',
				 label='New Gene')
		plt.xlabel("log(Average Protein Counts)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {num} proteins with the greatest difference in protein "
				  f"counts")
		plt.tight_layout()

	def gen_G4(self, nonshared_0_PCs, nonshared_0_PC_ids, nonshared_0_PC_idxs):
		"""
		Plots the proteins that have unique PC appearences (PCs in only one
		variant) that are likely due to subgeneral gene expression instances
		Args:
			nonfiltered_protein_counts: nonfiltered_protein_counts can contain
			new genes if desired, it may just throw off the scaling used to compare
			nonshared_0_PC_ids: ids of proteins with unique appearances
			nonshared_0_PC_idxs: idxs of proteins with unique appearances
		Returns: a graph of the proteins with 0 protein counts in one variant
		plotted from most PC counts (in the nonzero variant) to least
		"""
		# Order the protein count amounts by least to greatest:
		Ns_0_diff_counts = abs(nonshared_0_PCs[1] - nonshared_0_PCs[0])
		max_change_order = np.argsort(Ns_0_diff_counts)
		NS_0_PC_ids = []
		for idx in max_change_order:
			NS_0_PC_ids.append(nonshared_0_PC_ids[idx])
		avg_log_NS_0_PCs = self.get_LogData(nonshared_0_PC_idxs, nonshared_0_PCs, max_change_order)

		# Plot the results:
		plt.figure(figsize=(50, 60))
		plt.barh(NS_0_PC_ids, avg_log_NS_0_PCs[0], 0.1, align='edge', label='no New Gene')
		plt.barh(NS_0_PC_ids, avg_log_NS_0_PCs[1], -0.1, align='edge', label='New Gene')
		plt.xlabel("log(Average Protein Count)", fontweight='bold')
		plt.ylabel("Protein ID", fontweight='bold')
		plt.legend()
		plt.title(f"The {len(nonshared_0_PC_ids)} proteins with zero protein"
				  f"counts in one variant")
		plt.tight_layout()

	def gen_G5(self, filtered_PCs):
		# TODO: should i delete this since its so simple
		"""
		Generates the same graph as the gen_G2 function but with filtered PCs
		"""
		self.gen_G2(filtered_PCs)

	def gen_G6(self, min_num, filtered_PCs, filtered_ids, include_PC_change=0):
		"""
		Obtain proteins with the smallest absolute difference in protein counts
		Args:
			min_num: # of proteins to compare
			filtered_PCs: array of filtered PCs
			filtered_ids: array of filltered protein ids
			include_PC_change (optional): include the change in PCs from the
			experemental variant to the control variant (set  to 1 to include)
		Returns: A plot of the proteins with the smallest absolute change in PCs
		for those that have non-zero PCs in both variants
		"""
		min_diffs = abs(filtered_PCs[1] - filtered_PCs[0])
		min_PC_diff_changes = (filtered_PCs[1] - filtered_PCs[0])
		sortdiff = np.argsort(min_diffs)
		min_diff_idxs = sortdiff[:min_num]
		min_changes = []
		min_PC_diffs_ids = []
		#TODO: make this into a funciton going forward since I use the structure quite often?
		for idx in min_diff_idxs:
			min_changes.append(min_PC_diff_changes[idx])
			min_PC_diffs_ids.append(filtered_ids[idx])
		avg_log_min_PC_diffs = self.get_LogData(min_diff_idxs, filtered_PCs, min_diff_idxs)
		min_diff_var0 = avg_log_min_PC_diffs[0]
		min_diff_var1 = avg_log_min_PC_diffs[1]

		if include_PC_change == 0:
			plt.figure(figsize=(10, 6))
			plt.barh(min_PC_diffs_ids, min_diff_var0,
					 0.1, align='edge', label='no New Gene')
			plt.barh(min_PC_diffs_ids, min_diff_var1,
					 height=-0.1, align='edge', label='New Gene')
			plt.xlabel("log(Average Protein Count)", fontweight='bold')
			plt.ylabel("Protein ID", fontweight='bold')
			plt.legend()
			plt.title(f"The {min_num} proteins with the smallest change "
					  f" \nin in protein counts when a new gene is added ")
			plt.tight_layout()
		else:
			# plot with percent difference on the side
			fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
			# create labels for the % changes:
			min_change_str = []
			for m in range(min_num):
				min_value = ' ' + str(round(min_changes[m], 3))
				min_change_str.append(min_value)

			ax.barh(min_PC_diffs_ids, min_diff_var0,
					0.1, align='edge', label='no New Gene')
			ax.barh(min_PC_diffs_ids, min_diff_var1,
					height=-0.1, align='edge', label='New Gene')
			ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
			ax.set_ylabel("Protein ID", fontweight='bold')
			ax2 = ax.twinx()
			ax2.set_ylim(ax.get_ylim())
			min_change_str = np.array(min_change_str)
			ax2.set_yticks(np.arange(len(min_PC_diffs_ids)), labels=min_change_str)
			ax2.set_ylabel(r'$\Delta$ in Protein Count (from var1 to var0)', fontweight='bold')
			ax.legend()
			ax.set_title(f"The {min_num} proteins with the smallest absolute difference"
						 f" \nin protein counts when a new gene is added ")
			plt.tight_layout()

	def gen_G7(self, max_num, filtered_PCs, filtered_ids, include_PC_change=0):
		"""
		Obtains filtered proteins with the largest absolute difference in PCs
		Args:
			max_num: # of proteins to compare
			filtered_PCs: filtered protein count data
			filtered_ids: filtered proteins' ids
			include_PC_change (optional): plot the difference in PCs between
			the experiemental variant and control variant (set to 1 to turn on)
		Returns: a plot of the proteins with the largest absolute difference in
		PCs between variants
		"""
		max_diffs = abs(filtered_PCs[1] - filtered_PCs[0])
		max_PC_diff_changes = (filtered_PCs[1] - filtered_PCs[0])
		sortdiff = np.argsort(max_diffs)
		max_diff_idxs = sortdiff[-max_num:]
		max_changes = []
		max_PC_diffs_ids = []
		for idx in max_diff_idxs:
			max_changes.append(max_PC_diff_changes[idx])
			max_PC_diffs_ids.append(filtered_ids[idx])
		avg_log_max_PC_diffs = self.get_LogData(max_diff_idxs, filtered_PCs, max_diff_idxs)

		max_diff_var0 = avg_log_max_PC_diffs[0]
		max_diff_var1 = avg_log_max_PC_diffs[1]

		if include_PC_change == 0:
			plt.figure(figsize=(10, 6))
			plt.barh(max_PC_diffs_ids, max_diff_var0,
					 0.1, align='edge', label='no New Gene')
			plt.barh(max_PC_diffs_ids, max_diff_var1,
					 height=-0.1, align='edge', label='New Gene')

			plt.xlabel("log(Average Protein Count)")
			plt.ylabel("Protein ID", fontweight='bold')
			plt.legend()
			plt.title(f"The {max_num} proteins with the greatest change in protein"
					  f" \ncounts when a new gene is added")
			plt.tight_layout()
		else:
			# plot with percent difference on the side
			fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
			# create labels for the % changes:
			max_change_str = []
			for m in range(max_num):
				max_value = ' ' + str(round(max_changes[m], 2))
				max_change_str.append(max_value)

			ax.barh(max_PC_diffs_ids, max_diff_var0,
					0.1, align='edge', label='no New Gene')
			ax.barh(max_PC_diffs_ids, max_diff_var1,
					height=-0.1, align='edge', label='New Gene')
			ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
			ax.set_ylabel("Protein ID", fontweight='bold')
			ax2 = ax.twinx()
			ax2.set_ylim(ax.get_ylim())
			max_change_str = np.array(max_change_str)
			ax2.set_yticks(np.arange(len(max_PC_diffs_ids)), labels=max_change_str)
			ax2.set_ylabel(r'$\Delta$ in Protein Count (from var1 to var0)', fontweight='bold')
			ax.legend()
			ax.set_title(f"The {max_num} proteins with the largest absolute change"
						 f" \nin protein counts when a new gene is added "
						 )
			plt.tight_layout()

	def gen_G8(self, max_fold_num, filtered_PCs, filtered_ids, max_fold_num_woLogScale=0):
		"""
		Find proteins with the max fold increase in PCs from the control variant
		to the experimental variant
		Args:
			max_fold_num: # of proteins to compare on the plot
			filtered_PCs: filtered PC data
			filtered_ids: ids for filtered proteins
			max_fold_num_woLogScale (optional): by default this graph is plotted
			on a log scale; set this value equal to 1 to plot with out log scale
		Returns: a plot of the proteins with the largest fold increase in PCs
		from the control variant to the experiemental variant
		"""
		log_max_fold_PCs, max_fold_PCs_ids, max_fold_changes = self.get_max_fold_increase(max_fold_num, filtered_PCs, filtered_ids)
		max_fold_PC_values = log_max_fold_PCs

		if max_fold_num_woLogScale == 1:
			max_fold_PCs, max_fold_PCs_ids, max_fold_changes = self.get_max_fold_increase(max_fold_num, filtered_PCs, filtered_ids, 1)
			max_fold_PC_values = max_fold_PCs

		# Plot the results:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
		max_var0 = max_fold_PC_values[0]
		max_var1 = max_fold_PC_values[1]
		max_fold_val = [max(max_var0),
						max(max_var1)]
		max_fold_val = max(max_fold_val)
		max_fold_val = round(max_fold_val) + 1

		# create labels for the % changes:
		max_change_str = []
		for m in range(max_fold_num):
			max_value = '+' + str(round(max_fold_changes[m], 2)) + '%'
			max_change_str.append(max_value)

		ax.barh(max_fold_PCs_ids, max_var0,
				0.1, align='edge', label='no New Gene')
		ax.barh(max_fold_PCs_ids, max_var1,
				height=-0.1, align='edge', label='New Gene')
		ax.set_xlim([0, max_fold_val])
		for m, n in enumerate(max_var1):
			ax.text(n + .01, m - .15, max_change_str[m],
					color='black', fontweight='bold', size=10)
		ax.set_ylabel("Protein ID", fontweight='bold')
		ax.legend()
		if max_fold_num_woLogScale == 1:
			ax.set_xlabel("Average Protein Count", fontweight='bold')
		else:
			ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
		ax.set_title(f"The {max_fold_num} proteins with the greatest fold"
					 f" \n increase in protein count between variants"
					 )
		plt.tight_layout()

	def gen_G9(self, min_fold_num, filtered_PCs, filtered_ids, min_fold_num_woLogScale=0):
		"""
		Obtain proteins with the greatest fold decrease in protein counts from
		the control variant to the experimental variant
		Args:
			min_fold_num: # of proteins to compare on the plot
			filtered_PCs: array of the filtered proteins' PCs
			filtered_ids: filtered proteins' ids
			min_fold_num_woLogScale (optional): by default this graph is plotted
			on a log scale; set this value equal to 1 to plot with out log scale
		Returns: a plot of the proteins with the greatest fold increases in PCs
		from the control variant to the experiemental variant
		"""
		self.get_max_fold_decrease(min_fold_num, filtered_PCs, filtered_ids)
		log_min_fold_PCs, min_fold_PCs_ids, min_fold_changes = self.get_max_fold_decrease(min_fold_num, filtered_PCs, filtered_ids)
		min_fold_PC_values = log_min_fold_PCs

		if min_fold_num_woLogScale == 1:
			min_fold_PCs, min_fold_PCs_ids, min_fold_changes = self.get_max_fold_decrease(min_fold_num, filtered_PCs, filtered_ids, 1)
			min_fold_PC_values = min_fold_PCs

		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
		min_var0 = min_fold_PC_values[0]
		min_var1 = min_fold_PC_values[1]
		min_fold_val = [max(min_var0),
						max(min_var1)]
		min_fold_val = max(min_fold_val)
		min_fold_val = round(min_fold_val) + 1

		min_change_str = []
		for m in range(min_fold_num):
			min_value = ' ' + str(round(min_fold_changes[m], 2)) + '%'
			min_change_str.append(min_value)

		ax.barh(min_fold_PCs_ids, min_var0,
				0.1, align='edge', label='no New Gene')
		ax.barh(min_fold_PCs_ids, min_var1,
				height=-0.1, align='edge', label='New Gene')
		ax.set_xlim([0, min_fold_val])
		for m, n in enumerate(min_var1):
			ax.text(n + .01, m - .3, min_change_str[m],
					color='black', fontweight='bold', size=10)
		ax.set_ylabel("Protein ID", fontweight='bold')
		ax.legend()
		if min_fold_num_woLogScale == 1:
			ax.set_xlabel("Average Protein Count", fontweight='bold')
		else:
			ax.set_xlabel("log(Average Protein Count)", fontweight='bold')

		ax.set_title(f"The {min_fold_num} proteins with the largest fold "
						f" \ndecrease in protein counts between variants"
						)
		plt.tight_layout()

	def gen_G10(self, shared_diff_num, filtered_PCs, filtered_ids, PlotwLogScale=0):
		"""
		Obtain and plot the proteins with the greatest increases and
		decreases side by side so the magintudes can be compared to eachother
		Args:
			shared_diff_num: # of proteins to compare on the plot
			filtered_PCs: array of the filtered proteins' PCs
			filtered_ids: filtered proteins' ids
			PlotwLogScale (optional): if dsires to plot this on a log scale,
			set this equal to 1.
		Returns: a plot of greatest differences in PC increases and decreases
		"""
		max_diffs = filtered_PCs[1] - filtered_PCs[0]
		sortdiff = np.argsort(max_diffs)
		max_diff_idxs = sortdiff[-shared_diff_num:]
		max_changes = []
		max_PC_diffs_ids = []
		for idx in max_diff_idxs:
			max_changes.append(max_diffs[idx])
			max_PC_diffs_ids.append(filtered_ids[idx])
		avg_max_PC_diffs = np.zeros((len(self.variants), len(max_diff_idxs)))
		for variant in self.variants:
			for idx in range(len(max_diff_idxs)):
				index = max_diff_idxs[idx]
				avg_max_PC_diffs[variant][idx] = filtered_PCs[variant][index]

		min_diffs = filtered_PCs[0] - filtered_PCs[1]
		sortdiff = np.argsort(min_diffs)
		min_diff_idxs = sortdiff[-shared_diff_num:]
		min_changes = []
		min_PC_diffs_ids = []
		for idx in min_diff_idxs:
			min_changes.append(min_diffs[idx])
			min_PC_diffs_ids.append(filtered_ids[idx])
		avg_min_PC_diffs = np.zeros((len(self.variants), len(min_diff_idxs)))
		for variant in self.variants:
			for idx in range(len(min_diff_idxs)):
				index = min_diff_idxs[idx]
				avg_min_PC_diffs[variant][idx] = filtered_PCs[variant][index]

		if PlotwLogScale == 1:
			avg_log_max_PC_diffs = self.get_LogData(max_diff_idxs, filtered_PCs, max_diff_idxs)
			avg_max_PC_diffs = avg_log_max_PC_diffs
			avg_log_min_PC_diffs = self.get_LogData(min_diff_idxs, filtered_PCs, min_diff_idxs)
			avg_min_PC_diffs = avg_log_min_PC_diffs

		# Plot the max increases and decreases  together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
		# designate variables
		max_var0 = avg_max_PC_diffs[0]
		max_var1 = avg_max_PC_diffs[1]
		min_var0 = avg_min_PC_diffs[0]
		min_var1 = avg_min_PC_diffs[1]
		min_max_val = [max(min_var0), max(min_var1)]
		max_max_val = [max(max_var0), max(max_var1)]
		min_max_val = max(min_max_val)
		max_max_val = max(max_max_val)

		if PlotwLogScale == 1:
			max_val = [min_max_val, max_max_val]
			max_val = max(max_val)
			max_val = round(max_val) + 1
		else:
			min_max_val = round(min_max_val) + 10000
			max_max_val = round(max_max_val) + 1000

		max_change_str = []
		min_change_str = []
		for m in range(shared_diff_num):
			max_value = ' +' + str(round(max_changes[m], 2))
			min_value = ' -'+ str(round(min_changes[m], 2))
			max_change_str.append(max_value)
			min_change_str.append(min_value)
		# maximum change plot specs:
		ax[0].barh(max_PC_diffs_ids, max_var1,
				   align='center', label='New Gene', color='#17becf')
		ax[0].barh(max_PC_diffs_ids, max_var0,
				   align='center', label='no New Gene', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False,
						  labelleft=True, labelright=False)
		if PlotwLogScale == 1:
			ax[0].set_xlim([0, max_val])
			for m, n in enumerate(max_var1):
				ax[0].text(n + .05, m - .1, max_change_str[m],
						   color='black', fontweight='bold', size=6)
			ax[0].set_xlabel("Log(Average Protein Counts)", fontweight='bold')
		else:
			ax[0].set_xlim([0, max_max_val])
			for m, n in enumerate(max_var1):
				ax[0].text(n + 5, m - .1, max_change_str[m],
						   color='black', fontweight='bold', size=6)
			ax[0].set_xlabel("Average Protein Counts", fontweight='bold')
		ax[0].set_ylabel("Protein ID", fontweight='bold')
		ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[0].set_title(f"The {shared_diff_num} proteins with the greatest "
						f"increase \n in protein count with the addition "
						f"of a new gene", size=10)
		# minimum change plot specs:
		ax[1].barh(min_PC_diffs_ids, min_var0,
				   align='center', label='no New Gene', color='#bcbd22')
		ax[1].barh(min_PC_diffs_ids, min_var1,
				   align='center', label='New Gene', color='#17becf')
		ax[1].tick_params(bottom=True, top=False, left=False, right=True)
		ax[1].tick_params(axis='both', which='major', labelsize=6)
		ax[1].tick_params(labelbottom=True, labeltop=False,
						  labelleft=False, labelright=True)
		if PlotwLogScale == 1:
			ax[1].set_xlim([0, max_val])
			for m, n in enumerate(min_var1):
				ax[1].text(n - .9, m - .1, min_change_str[m],
						   color='black', fontweight='bold', size=6)
			ax[1].set_xlabel("Log(Average Protein Counts)", fontweight='bold')
		else:
			ax[1].set_xlim([0, min_max_val])
			for m, n in enumerate(min_var1):
				ax[1].text(n + 5, m - .1, min_change_str[m],
						   color='black', fontweight='bold', size=6)
			ax[1].set_xlabel("Average Protein Counts", fontweight='bold')
		ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[1].set_title(f"The {shared_diff_num} proteins with the greatest "
						f"decrease \n in protein count with the addition "
						f"of a new gene", size=10)
		fig.suptitle("Maximum Increases and Decreases in Protein Counts "
					 "\nwith the Addition of a New Gene")
		plt.tight_layout()

	def gen_G11(self, sharednum, filtered_PCs, filtered_ids):
		"""
		Obtain and plot the proteins with the greatest fold increases and
		decreases side by side so the magintudes can be compared to eachother
		Args:
			sharednum: # of proteins to compare on the plot
			filtered_PCs: array of the filtered proteins' PCs
			filtered_ids: filtered proteins' ids
		Returns: a plot of the proteins with the greatest fold increases and
		decreases compared side by side with each other
		"""
		max_fold_PC_values, max_fold_PCs_ids, max_fold_changes = self.get_max_fold_increase(sharednum, filtered_PCs, filtered_ids)
		min_fold_PC_values, min_fold_PCs_ids, min_fold_changes = self.get_max_fold_decrease(sharednum, filtered_PCs, filtered_ids)
		# Plot the max and min fold changes together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
		# designate variables
		max_var0 = max_fold_PC_values[0]
		max_var1 = max_fold_PC_values[1]
		min_var0 = min_fold_PC_values[0]
		min_var1 = min_fold_PC_values[1]

		# create labels for the % changes:
		max_change_str = []
		min_change_str = []
		for m in range(sharednum):
			max_val = '+' + str(round(max_fold_changes[m], 2)) + '%'
			min_val = ' ' + str(round(min_fold_changes[m], 2)) + '%'
			max_change_str.append(max_val)
			min_change_str.append(min_val)

		# maximum fold change plot specs:
		ax[0].barh(max_fold_PCs_ids, max_var1, align='center',
				   label='New Gene', color='#17becf')
		ax[0].barh(max_fold_PCs_ids, max_var0, align='center',
				   label='no New Gene', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False,
						  labelleft=True, labelright=False)
		ax[0].set_xlim([0, 3])
		for m, n in enumerate(max_var1):
			ax[0].text(n + .05, m - .1, max_change_str[m], color='black',
					   fontweight='bold', size=6)
		ax[0].set_xlabel("log(Average Protein Counts)", fontweight='bold')
		ax[0].set_ylabel("Protein ID", fontweight='bold')
		# ax[0].legend(loc='best')
		ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[0].set_title(f"The {sharednum} proteins with the greatest fold "
						f"increase \n in protein counts with the"
						f" addition of a New Gene", size=10)
		# minimum fold change plot specs:
		ax[1].barh(min_fold_PCs_ids, min_var0, align='center',
				   label='no New Gene', color='#bcbd22')
		ax[1].barh(min_fold_PCs_ids, min_var1, align='center',
				   label='New Gene', color='#17becf')
		ax[1].tick_params(bottom=True, top=False, left=False, right=True)
		ax[1].tick_params(axis='both', which='major', labelsize=6)
		ax[1].tick_params(labelbottom=True, labeltop=False,
						  labelleft=False, labelright=True)
		ax[1].set_xlim([0, 3])
		for m, n in enumerate(min_var1):
			ax[1].text(n + .05, m - .1, min_change_str[m], color='black',
					   fontweight='bold', size=6)
		ax[1].set_xlabel("log(Average Protein Counts)", fontweight='bold')
		# ax[1].legend(loc='best')
		ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[1].set_title(f"The {sharednum} proteins with the greatest fold "
						f"decrease \n in protein counts with the addition of"
						f" a New Gene", size=10)
		fig.suptitle("Maximum and Minimum Fold Change Comparisons for"
					 " Protein Counts "
					 "\nwith the Addition of a New Gene")
		plt.tight_layout()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		"""
		Call graph generating functions
		"""
		# define/initialize commonly used variables
		self.variants = self.ap.get_variants()
		self.n_total_gens = self.ap.n_generation
		self.all_monomer_ids = []
		self.original_monomer_ids = []
		self.new_gene_monomer_ids = []
		self.total_protein_counts = []
		protein_counts = self.generate_data(simDataFile)
		monomer_idx_dict_PreFilter = {monomer: i for i, monomer in
									  enumerate(self.all_monomer_ids)}
		# Plot 1
		if randnum > 0:
			self.gen_G1(randnum)
			exportFigure(plt, plotOutDir, plotOutFileName + '_PCs_for_'
						 + str(randnum) + '_random_proteins_wNG_noFilter', metadata)
			plt.close('all')

		# Plots 2a and 2b
		if include_NG_G2 == 0:
			PCs = protein_counts
			IDs = self.original_monomer_ids
			words = '_original_PC_comparisons_woNG'
		else:
			PCs = self.total_protein_counts
			IDs = self.all_monomer_ids
			words = '_original_PC_comparisons_wNG'

		if var_PC_comparison == 1:
			self.gen_G2(PCs)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(PCs[0])) +
						 words,
						 metadata)
			plt.close('all')

		if var_PC_comparison_LogScale == 1:
			PC_LogData_idxs = self.get_idxs(IDs)
			PC_LogData = self.get_LogData(PC_LogData_idxs, PCs)
			self.gen_G2(PC_LogData)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(PCs[0])) +
						 words +'_LogScale',
						 metadata)
			plt.close('all')

		# Plot 3:
		if num > 0:
			self.gen_G3(num, protein_counts, monomer_idx_dict_PreFilter)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(num) +
						 '_largest_absolute_PC_diffs_wNG_noFilter',
						 metadata)
			plt.close('all')

		# Plot 4:
		if unique_PC_appearances == 1:
			nonshared_0_PCs, nonshared_0_PC_ids, nonshared_0_PC_idxs = self.find_unique_proteins(protein_counts, monomer_idx_dict_PreFilter)
			self.gen_G4(nonshared_0_PCs, nonshared_0_PC_ids, nonshared_0_PC_idxs)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' +
						 str(len(nonshared_0_PC_ids)) +
						 '_unique_PC_appearances_noNG_noFilter', metadata)
			plt.close('all')

		# Manditory Data Filtration
		F_PCs, F_PC_ids, F_PC_idxs = self.filter_data(protein_counts, monomer_idx_dict_PreFilter, filter_num)
		F_monomer_idx_dict = {monomer: i for i, monomer in enumerate(F_PC_ids)}

		# Plot 5:
		if var_PC_comparison_wF == 1:
			self.gen_G5(F_PCs)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(F_PCs[0])) +
						 '_PC_comparisons_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')

		if var_PC_comparison_wF_LogScale == 1:
			F_PC_LogData = self.get_LogData(F_PC_idxs, F_PCs)
			self.gen_G5(F_PC_LogData)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(F_PC_LogData[0])) +
						 '_PC_comparisons_LogScale_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')

		# Plot 6
		if min_num > 0:
			self.gen_G6(min_num, F_PCs, F_PC_ids)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(min_num) +
						 '_min_PC_diffs_woD_Filter_'
						 + str(filter_num), metadata)
			plt.close('all')
		if show_PC_diff_6b == 1:
			self.gen_G6(min_num, F_PCs, F_PC_ids, show_PC_diff_6b)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(min_num) +
						 '_min_PC_diffs_wD_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')

		# Plot 7:
		if max_num > 0:
			self.gen_G7(max_num, F_PCs, F_PC_ids)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(max_num) +
						 '_max_PC_diffs_woD_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')
		if show_PC_diff_7b == 1:
			self.gen_G7(max_num, F_PCs, F_PC_ids, show_PC_diff_7b)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(max_num) +
						 '_max_PC_diffs_wD_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')

		# Plot 8
		if max_fold_num > 0:
			self.gen_G8(max_fold_num, F_PCs, F_PC_ids)
			exportFigure(plt, plotOutDir,
						 plotOutFileName + '_' + str(max_fold_num) +
						 '_max_PC_fold_increases_LogScale_wPD_Filter_' +
						 str(filter_num),
						 metadata)
			plt.close('all')
			if max_fold_num_woLogScale == 1:
				self.gen_G8(max_fold_num, F_PCs, F_PC_ids, 1)
				exportFigure(plt, plotOutDir,
						 	plotOutFileName + '_' + str(max_fold_num) +
						 	'_max_PC_fold_increases_wPD_Filter_' +
						 	str(filter_num),
						 	metadata)
				plt.close('all')

		# Plot 9
		if min_fold_num > 0:
			self.gen_G9(min_fold_num, F_PCs, F_PC_ids)
			exportFigure(plt, plotOutDir,
						 plotOutFileName + '_' + str(min_fold_num) +
						 '_max_PC_fold_decreases_LogScale_wPD_Filter_' +
						 str(filter_num),
						 metadata)
			plt.close('all')
			if min_fold_num_woLogScale == 1:
				self.gen_G9(min_fold_num, F_PCs, F_PC_ids, 1)
				exportFigure(plt, plotOutDir,
						 	plotOutFileName + '_' + str(min_fold_num) +
						 	'_max_PC_fold_decreases_wPD_Filter_' +
						 	str(filter_num),
						 	metadata)
				plt.close('all')

		# Plot 10
		if shared_diff_num > 0:
			self.gen_G10(shared_diff_num, F_PCs, F_PC_ids)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' +
						 str(shared_diff_num) +
						 '_max_PC_diff_comparisons_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')
			if shared_diff_LogScale == 1:
				self.gen_G10(shared_diff_num, F_PCs, F_PC_ids, 1)
				exportFigure(plt, plotOutDir, plotOutFileName + '_' +
							 str(shared_diff_num) +
							 '_max_PC_diff_comparisons_LogScale_Filter_' + str(filter_num),
							 metadata)
				plt.close('all')

		# Plot 11
		if sharednum > 0:
			self.gen_G11(sharednum, F_PCs, F_PC_ids)
			exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(sharednum) +
						 '_max_PC_fold_comparisons_Filter_' + str(filter_num),
						 metadata)
			plt.close('all')

if __name__ == "__main__":
	Plot().cli()
