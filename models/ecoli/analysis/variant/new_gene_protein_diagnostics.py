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

# 1 to exclude cells that took full MAX_CELL_LENGTH, 0 otherwise
exclude_timeout_cells = 1

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

MAX_CELL_LENGTH = 180
if (exclude_timeout_cells==0):
	MAX_CELL_LENGTH += 1000000

MAX_YLIM_PLOT = MAX_CELL_LENGTH

"""
designate any numbers for how many 
proteins should be plotted in the plots here:
"""
# Number of random proteins to observe the protein counts (PCs) between variants:
randnum = 10
# Number of proteins with the largest PC increase between variants
# to observe (this one will include GFP and other variants with zero PCs
# in a variant:
num = 10
# Number of proteins with the smallest PC change between variants:
min_num = 20
# Number of proteins with the greatest change in PCs between variants
# to observe (not including GFP and proteins with zero counts in a variant):
max_num = 20
# Number of proteins to visualize the max increases and max decreases in
# protein counts between variants side by side:
shared_diff_num = 10
# Number of proteins to observe with the greatest fold increase in protein
# counts when GFP is added:
max_fold_num = 12
# Number of proteins to observe with the largest fold decrease in
# protein counts (not including GFP):
min_fold_num = 200
# Number of proteins to observe the max fold increase and decrease in protein
# counts side by side:
sharednum = 10
# change this number to be the minimum PC value to plot with:
filter_num = 1


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[
			mRNA_sim_data['is_new_gene']]['id'].tolist()
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
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids),\
			'number of new gene monomers and mRNAs should be equal'


		# Randomly select proteins to observe with the new protein(s):

		all_monomer_ids = monomer_sim_data['id']
		original_monomer_ids = np.delete(all_monomer_ids, np.where(
			all_monomer_ids == new_gene_monomer_ids))
		rand_monomers = np.random.choice(original_monomer_ids,
										 size = randnum, replace = False)
		monomer_names = np.append(rand_monomers, new_gene_monomer_ids)

		# Find gene indexes:
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(all_monomer_ids)}
		protein_idxs = [monomer_idx_dict.get(monomer_id)
									for monomer_id in monomer_names]
		protein_idxs = np.array(protein_idxs)

		variants = self.ap.get_variants()
		n_total_gens = self.ap.n_generation

		# Data extraction
		each_gen_avg_monomer_counts = [{} for id_ in monomer_names]
		all_avg_monomer_counts = [{} for id_ in monomer_names]
		each_gen_avg_log_monomer_counts = [{} for id_ in monomer_names]
		all_avg_log_monomer_counts = [{} for id_ in monomer_names]

		for variant in variants:
			if variant >= MAX_VARIANT:
				continue

			all_cells = self.ap.get_cells(
				variant=[variant],
				generation=np.arange(IGNORE_FIRST_N_GENS, n_total_gens),
				only_successful=True)

			if len(all_cells) == 0:
				continue

			# Load Data
			gen_avg_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				fun=lambda x: np.mean(x[:, protein_idxs], axis=0))
			avg_monmer_counts = np.mean(gen_avg_monomer_counts, axis=0)

			for m in range(len(protein_idxs)):
				each_gen_avg_monomer_counts[m][variant] \
					= gen_avg_monomer_counts[:,m]
				each_gen_avg_log_monomer_counts[m][variant] = \
					np.log10(gen_avg_monomer_counts[:,m] + 1)
				all_avg_monomer_counts[m][variant] = avg_monmer_counts[m]
				all_avg_log_monomer_counts[m][variant] = \
					np.log10(avg_monmer_counts[m] + 1)

		# Make plot for the random protein selections:
		plt.figure(figsize = (10, 6)) # 8.5, 11 is good for 100 genes!

		logdata = all_avg_log_monomer_counts
		ld=np.zeros((len(variants), len(protein_idxs)))
		for variant in variants:
			vardata=np.zeros((len(protein_idxs)))
			for m in range(len(protein_idxs)):
				data = logdata[m][variant]
				vardata[m] = data
			ld[variant] = vardata

		var0 = ld[0]
		var1 = ld[1]

		plt.barh(monomer_names, var0, 0.1, align = 'edge',
				 label='no GFP')
		plt.barh(monomer_names, var1, -0.1, align = 'edge', label='GFP')

		plt.xlabel("log(Average Protein Count)", fontweight='bold')
		plt.ylabel("Protein ID", fontweight='bold')
		plt.legend()
		plt.title("Protein count comparisons between variants "
				  "\n for "+str(randnum)
				  +" randomly selected proteins ")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_PCs_for_'
					 +str(randnum) +'_random_proteins_wGFP_noFilter', metadata)
		plt.close('all')



		# Observe the largest protein count size differences between variants:

		protein_counts = np.zeros((len(variants), len(original_monomer_ids)))
		total_protein_counts = np.zeros((len(variants), len(all_monomer_ids)))
		for variant in variants:
			all_cells = self.ap.get_cells(variant=[variant],
										  generation=np.arange(IGNORE_FIRST_N_GENS,
															   n_total_gens),
										  only_successful=True)

			# Get the protein counts for each gene/protein
			gen_avg_monomer_counts = read_stacked_columns(all_cells,
				'MonomerCounts', 'monomerCounts',
				fun=lambda x: np.mean(x[:], axis=0))
			total_avg_gene_counts = np.mean(gen_avg_monomer_counts, axis=0)
			total_protein_counts[variant] = total_avg_gene_counts
			old_gene_idxs = [monomer_idx_dict.get(monomer_id)
							 for monomer_id in original_monomer_ids]
			avg_gene_monomer_counts = total_avg_gene_counts[old_gene_idxs]
			protein_counts[variant] = avg_gene_monomer_counts

		# compare biggest changes between variants:
		protein_counts = np.array(protein_counts)
		total_protein_counts = np.array(total_protein_counts)

		# Retreive the indicies and ids of the proteins with the maximum count
		# differences between variants (including the newly added genes):
		diff = abs(protein_counts[1] - protein_counts[0])
		sortdiff = np.argsort(diff)
		maxdiff = sortdiff[-num:]
		inv_monomer_idx_dict = {idx: i for i, idx in monomer_idx_dict.items()}
		interest_proteins = [
			inv_monomer_idx_dict.get(monomer_id) for monomer_id in maxdiff]
		interest_proteins = np.append(interest_proteins, new_gene_monomer_ids)
		interest_protein_idxs = [monomer_idx_dict.get(monomer_id)
									for monomer_id in interest_proteins]
		interest_protein_counts = total_protein_counts[:, interest_protein_idxs]

		# get the log values of the proteins of interest:
		avg_log_interest_proteins = np.zeros((
			len(variants), len(interest_proteins)))
		for variant in variants:
			for m in range(len(interest_protein_idxs)):
				avg_log_interest_proteins[variant][m] = \
					np.log10(interest_protein_counts[variant][m] + 1)

		# Make plot for genes with the largest protein count differences:
		plt.figure(figsize = (10, 6))

		var0 = avg_log_interest_proteins[0]
		var1 = avg_log_interest_proteins[1]

		plt.barh(interest_proteins, var0, 0.1, align='edge',
				 label='no GFP')
		plt.barh(interest_proteins, var1, -0.1, align='edge',
				 label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {num} proteins with the greatest difference in protein "
				  f"counts")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(num) +
					 '_largest_absolute_PC_diffs_wGFP_noFilter',
					 metadata)
		plt.close('all')

		# Observe the differences between var 0 and var 1 counts:
		var0_x = protein_counts[0]
		var1_y = protein_counts[1]

		plt.figure(figsize= (10,10))
		plt.scatter(var0_x, var1_y, 1)
		m, b = np.polyfit(var0_x, var1_y, 1)
		plt.plot(var0_x, m*var0_x +b, linewidth=.5, color='#bcbd22')
		legstr = "linear fit: y = " + str(round(m,2)) + "x + " + str(round(b,2))
		plt.legend(["PC data", legstr])
		plt.xlabel("variant 0 (no GFP)")
		plt.ylabel("variant 1 (GFP)")
		plt.title(f"The {len(var0_x)} proteins plotted against each other")
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(var0_x)) +
					 '_original_PC_comparisons',
					 metadata)
		plt.close('all')

		# Observe the differences between var 0 and var 1 counts on a log scale:
		avg_log_interest_proteins = np.zeros((
			len(variants), len(protein_counts[0])))
		for variant in variants:
			for m in range(len(protein_counts[0])):
				avg_log_interest_proteins[variant][m] = \
					np.log10(protein_counts[variant][m] + 1)

		var0_x = avg_log_interest_proteins[0]
		var1_y = avg_log_interest_proteins[1]

		plt.figure(figsize=(10, 10))
		plt.scatter(var0_x, var1_y, 1)
		m, b = np.polyfit(var0_x, var1_y, 1)
		plt.plot(var0_x, m * var0_x + b, linewidth=.5, color='#bcbd22')
		legstr = "linear fit: y = " + str(round(m, 2)) + "x + " + str(round(b, 2))
		plt.legend(["PC data", legstr])
		plt.xlabel("log(variant 0 (no GFP))")
		plt.ylabel("log(variant 1 (GFP))")
		plt.title(f"The {len(var0_x)} proteins plotted against each other")
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(var0_x)) +
					 '_original_PC_comparisons_LogScale',
					 metadata)
		plt.close('all')


		# FILTER OUT ZEROS!

		# Extract and view the proteins that have nonzero counts for a variant:
		p_counts = np.array(protein_counts) # not including newly added proteins

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
										   var1_unique_0_PC_idxs, axis = 0)
		nonshared_0_PC_ids = [inv_monomer_idx_dict.get(monomer_id)
								 for monomer_id in nonshared_0_PC_idxs]
		nonshared_0_PCs = p_counts[:, nonshared_0_PC_idxs]


		# Order the protein count amounts by least to greatest:
		Ns_0_diff_counts = abs(nonshared_0_PCs[1] - nonshared_0_PCs[0])
		max_change_order = np.argsort(Ns_0_diff_counts)
		NS_0_PC_ids = []
		for idx in max_change_order:
			NS_0_PC_ids.append(nonshared_0_PC_ids[idx])
		avg_log_NS_0_PCs = np.zeros((len(variants), len(nonshared_0_PC_idxs)))
		for variant in variants:
			for idx in range(len(nonshared_0_PC_idxs)):
				index = max_change_order[idx]
				avg_log_NS_0_PCs[variant][idx] = \
					np.log10(nonshared_0_PCs[variant][index] + 1)


		# Plot the results:
		plt.figure(figsize=(50, 60))
		var0 = avg_log_NS_0_PCs[0]
		var1 = avg_log_NS_0_PCs[1]

		plt.barh(NS_0_PC_ids, var0, 0.1, align='edge', label='no GFP')
		plt.barh(NS_0_PC_ids, var1, -0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)", fontweight='bold')
		plt.ylabel("Protein ID", fontweight='bold')
		plt.legend()
		plt.title(f"The {len(nonshared_0_PC_ids)} proteins with zero protein"
				  f"counts in one variant")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' +
					 str(len(nonshared_0_PC_ids)) +
					 '_unique_PC_appearances_noGFP_noFilter', metadata)
		plt.close('all')


		# Filter data to obtain all non-zero :
		# Extract all proteins with non-zero protein counts in both variants:
		nonzero_p_counts_var0_idxs = np.nonzero(p_counts[0])
		nonzero_p_counts_var1_idxs = np.nonzero(p_counts[1])
		shared_nonzero_PCs_idxs = np.intersect1d(nonzero_p_counts_var0_idxs,
												 nonzero_p_counts_var1_idxs)
		nonzero_PCs = p_counts[:, shared_nonzero_PCs_idxs]
		nonzero_PCs_ids = [inv_monomer_idx_dict.get(monomer_id)
						   for monomer_id in shared_nonzero_PCs_idxs]
		nonzero_ids = all_monomer_ids[shared_nonzero_PCs_idxs]


		# filter again to a chosen minimum value:
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

		# Filtered x vs y comparisons of protein counts:

		# Observe the differences between var 0 and var 1 counts:
		var0_x = nonzero_PCs[0]
		var1_y = nonzero_PCs[1]

		plt.figure(figsize=(10, 10))
		plt.scatter(var0_x, var1_y, 1)
		m, b = np.polyfit(var0_x, var1_y, 1)
		plt.plot(var0_x, m * var0_x + b, linewidth=.5, color='#bcbd22')
		legstr = "linear fit: y = " + str(round(m, 2)) + "x + " + str(round(b, 2))
		plt.legend(["PC data", legstr])
		plt.xlabel("variant 0 (no GFP)")
		plt.ylabel("variant 1 (GFP)")
		plt.title(f"{len(var0_x)} proteins plotted against each other")
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(var0_x)) +
					 '_PC_comparisons_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')

		# Observe the differences between var 0 and var 1 counts on a log scale:
		avg_log_nonzero_PCs = np.zeros((len(variants), len(nonzero_PCs[0])))
		for variant in variants:
			for idx in range(len(nonzero_PCs[0])):
				avg_log_nonzero_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][idx] + 1)

		var0_x = avg_log_nonzero_PCs[0]
		var1_y = avg_log_nonzero_PCs[1]

		plt.figure(figsize=(10, 10))
		plt.scatter(var0_x, var1_y, 1)
		m, b = np.polyfit(var0_x, var1_y, 1)
		plt.plot(var0_x, m * var0_x + b, linewidth=.5, color='#bcbd22')
		legstr = "linear fit: y = " + str(round(m, 2)) + "x + " + str(round(b, 2))
		plt.legend(["PC data", legstr])
		plt.xlabel("log(variant 0 (no GFP))")
		plt.ylabel("log(variant 1 (GFP))")
		plt.title(f"{len(var0_x)} proteins plotted against each other")
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(len(var0_x)) +
					 '_PC_comparisons_LogScale_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')

		# Observe which protiens have the smallest decrease in
		# protein counts with the addition of GFP:

		min_diffs = abs(nonzero_PCs[1] - nonzero_PCs[0])
		min_PC_diff_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		sortdiff = np.argsort(min_diffs)
		min_diff_idxs = sortdiff[:min_num]
		min_changes = []
		min_PC_diffs_ids = []
		for idx in min_diff_idxs:
			min_changes.append(min_PC_diff_percents[idx])
			min_PC_diffs_ids.append(nonzero_PCs_ids[idx])
		avg_log_min_PC_diffs = np.zeros((len(variants), len(min_diff_idxs)))
		for variant in variants:
			for idx in range(len(min_diff_idxs)):
				index = min_diff_idxs[idx]
				avg_log_min_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Make plot for genes with the smallest protein count differences:
		plt.figure(figsize=(10, 6))

		# create labels for the % changes:
		min_change_str = []
		for m in range(min_num):
			min_value = ' ' + str(round(min_changes[m], 3)) + '%'
			min_change_str.append(min_value)

		min_diff_var0 = avg_log_min_PC_diffs[0]
		min_diff_var1 = avg_log_min_PC_diffs[1]

		plt.barh(min_PC_diffs_ids, min_diff_var0,
				0.1, align='edge', label='no GFP')
		plt.barh(min_PC_diffs_ids, min_diff_var1,
				height=-0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)", fontweight='bold')
		plt.ylabel("Protein ID", fontweight='bold')
		plt.legend()
		plt.title(f"The {min_num} proteins with the smallest change "
					 f" \nin in protein counts when GFP is added "
					 )

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(min_num) +
					 '_min_PC_diffs_woGFP_woPD_Filter_'
					 + str(filter_num), metadata)
		plt.close('all')

		# plot with percent difference on the side
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

		ax.barh(min_PC_diffs_ids, min_diff_var0,
				0.1, align='edge', label='no GFP')
		ax.barh(min_PC_diffs_ids, min_diff_var1,
				height=-0.1, align='edge', label='GFP')


		ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
		ax.set_ylabel("Protein ID", fontweight='bold')

		ax2 = ax.twinx()
		ax2.set_ylim(ax.get_ylim())
		min_change_str = np.array(min_change_str)
		ax2.set_yticks(np.arange(len(min_PC_diffs_ids)), labels=min_change_str)
		ax2.set_ylabel('% difference', fontweight='bold')

		ax.legend()
		ax.set_title(f"The {min_num} proteins with the smallest change "
					 f" \nin in protein counts when GFP is added "
					 )

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(min_num) +
					 '_min_PC_diffs_woGFP_wPD_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')


		# Observe which protiens have the greatest change in protein counts
		# with the addition of GFP (but not including GFP or other proteins
		# with zero PCs for a variant):

		max_diffs = abs(nonzero_PCs[1] - nonzero_PCs[0])
		max_PC_diff_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		sortdiff = np.argsort(max_diffs)
		max_diff_idxs = sortdiff[-max_num:]
		max_changes = []
		max_PC_diffs_ids = []
		for idx in max_diff_idxs:
			max_changes.append(max_PC_diff_percents[idx])
			max_PC_diffs_ids.append(nonzero_PCs_ids[idx])
		avg_log_max_PC_diffs = np.zeros((len(variants), len(max_diff_idxs)))
		for variant in variants:
			for idx in range(len(max_diff_idxs)):
				index = max_diff_idxs[idx]
				avg_log_max_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Make plot for genes with the largest protein count differences:
		plt.figure(figsize=(10, 6))
		max_diff_var0 = avg_log_max_PC_diffs[0]
		max_diff_var1 = avg_log_max_PC_diffs[1]

		plt.barh(max_PC_diffs_ids, max_diff_var0,
				 0.1, align='edge', label='no GFP')
		plt.barh(max_PC_diffs_ids, max_diff_var1,
				 height=-0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID", fontweight='bold')
		plt.legend()
		plt.title(f"The {max_num} proteins with the greatest change in protein"
				  f" \ncounts when GFP is added")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(max_num) +
					 '_max_PC_diffs_woGFP_woPD_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')

		# plot with percent difference on the side

		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

		# create labels for the % changes:
		max_change_str = []
		for m in range(max_num):
			max_value = ' ' + str(round(max_changes[m], 2)) + '%'
			max_change_str.append(max_value)

		ax.barh(max_PC_diffs_ids, max_diff_var0,
				0.1, align='edge', label='no GFP')
		ax.barh(max_PC_diffs_ids, max_diff_var1,
				height=-0.1, align='edge', label='GFP')

		ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
		ax.set_ylabel("Protein ID", fontweight='bold')

		ax2 = ax.twinx()
		ax2.set_ylim(ax.get_ylim())
		min_change_str = np.array(max_change_str)
		ax2.set_yticks(np.arange(len(max_PC_diffs_ids)), labels=max_change_str)
		ax2.set_ylabel('% difference', fontweight='bold')

		ax.legend()
		ax.set_title(f"The {max_num} proteins with the largest change "
					 f" \nin in protein counts when GFP is added "
					 )

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(min_num) +
					 '_max_PC_diffs_woGFP_wPD_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')



		# Observe the max and min changes in protein counts side by side:

		max_PC_diffs = (nonzero_PCs[1] - nonzero_PCs[0])
		max_PC_diff_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		sortdiff = np.argsort(max_PC_diffs)
		max_diff_idxs = sortdiff[-shared_diff_num:]
		max_PC_diff_changes = []
		max_PC_diffs_ids = []
		numbersmax = []
		for idx in max_diff_idxs:
			max_PC_diff_changes.append(max_PC_diff_percents[idx])
			max_PC_diffs_ids.append(nonzero_PCs_ids[idx])
			numbersmax.append(max_PC_diffs[idx])
		avg_log_max_PC_diffs = np.zeros((len(variants), len(max_diff_idxs)))
		avg_max_PC_diffs = np.zeros((len(variants), len(max_diff_idxs)))
		for variant in variants:
			for idx in range(len(max_diff_idxs)):
				index = max_diff_idxs[idx]
				avg_max_PC_diffs[variant][idx] = nonzero_PCs[variant][index]
				avg_log_max_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		min_PC_diffs = (nonzero_PCs[0] - nonzero_PCs[1])
		min_PC_diff_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		sortdiff = np.argsort(min_PC_diffs)
		min_diff_idxs = sortdiff[-shared_diff_num:]
		min_PC_diff_changes = []
		min_PC_diff_ids = []
		numbersmin = []
		for idx in min_diff_idxs:
			min_PC_diff_changes.append(min_PC_diff_percents[idx])
			min_PC_diff_ids.append(nonzero_PCs_ids[idx])
			numbersmin.append(min_PC_diffs[idx])
		avg_log_min_PC_diffs = np.zeros((len(variants), len(min_diff_idxs)))
		avg_min_PC_diffs = np.zeros((len(variants), len(min_diff_idxs)))
		for variant in variants:
			for idx in range(len(min_diff_idxs)):
				index = min_diff_idxs[idx]
				avg_min_PC_diffs[variant][idx] = nonzero_PCs[variant][index]
				avg_log_min_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the max increases and decreases  together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))

		# Name the results:
		max_var0 = avg_max_PC_diffs[0]
		max_var1 = avg_max_PC_diffs[1]
		min_var0 = avg_min_PC_diffs[0]
		min_var1 = avg_min_PC_diffs[1]

		min_max_val = [max(avg_min_PC_diffs[0]), max(avg_min_PC_diffs[1])]
		max_max_val = [max(avg_max_PC_diffs[0]), max(avg_max_PC_diffs[1])]
		min_max_val = max(min_max_val)
		max_max_val = max(max_max_val)
		min_max_val = round(min_max_val) + 10000
		max_max_val = round(max_max_val) + 1000

		# create labels for the % changes:
		max_change_str = []
		min_change_str = []
		for m in range(shared_diff_num):
			max_value = '+' + str(round(max_PC_diff_changes[m], 2)) + '%'
			min_value = str(round(min_PC_diff_changes[m], 2)) + '%'
			max_change_str.append(max_value)
			min_change_str.append(min_value)

		max_change_str1 = []
		min_change_str1 = []
		for m in range(shared_diff_num):
			max_value = '+ ' + str(round(numbersmax[m], 2))
			min_value = '- ' + str(round(numbersmin[m], 2))
			max_change_str1.append(max_value)
			min_change_str1.append(min_value)

		# maximum fold change plot specs:
		ax[0].barh(max_PC_diffs_ids, max_var1,
				   align='center', label='GFP', color='#17becf')
		ax[0].barh(max_PC_diffs_ids, max_var0,
				   align='center', label='no GFP', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False,
						  labelleft=True, labelright=False)
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
						f"of GFP",size=10)

		# minimum fold change plot specs:
		ax[1].barh(min_PC_diff_ids, min_var0,
				   align='center', label='no GFP', color='#bcbd22')
		ax[1].barh(min_PC_diff_ids, min_var1,
				   align='center', label='GFP', color='#17becf')
		ax[1].tick_params(bottom=True, top=False, left=False, right=True)
		ax[1].tick_params(axis='both', which='major', labelsize=6)
		ax[1].tick_params(labelbottom=True, labeltop=False,
						  labelleft=False, labelright=True)
		ax[1].set_xlim([0, min_max_val])
		for m, n in enumerate(min_var1):
			ax[1].text(n + 5, m - .1, min_change_str[m],
					   color='black', fontweight='bold', size=6)
		ax[1].set_xlabel("Average Protein Counts", fontweight='bold')
		ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[1].set_title(f"The {shared_diff_num} proteins with the greatest "
						f"decrease \n in protein count with the addition "
						f"of GFP", size=10)

		fig.suptitle("Maximum Increases and Decreases in Protein Counts "
					 "\nwith the Addition of GFP")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' +
					 str(shared_diff_num) +
					 '_max_PC_diff_comparisons_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')



		# Plot the maximum increases and decreases on a log scale together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))

		# Name the results:
		max_var0 = avg_log_max_PC_diffs[0]
		max_var1 = avg_log_max_PC_diffs[1]
		min_var0 = avg_log_min_PC_diffs[0]
		min_var1 = avg_log_min_PC_diffs[1]

		max_val = [max(avg_log_min_PC_diffs[0]), max(avg_log_min_PC_diffs[1]),
				   max(avg_log_max_PC_diffs[0]), max(avg_log_max_PC_diffs[1])]
		max_val = max(max_val)
		max_val = round(max_val) + 1

		# create labels for the % changes:
		max_change_str = []
		min_change_str = []
		for m in range(shared_diff_num):
			max_value = '+' + str(round(max_PC_diff_changes[m], 2)) + '%'
			min_value = str(round(min_PC_diff_changes[m], 2)) + '%'
			max_change_str.append(max_value)
			min_change_str.append(min_value)

		# maximum fold change plot specs:
		ax[0].barh(max_PC_diffs_ids, max_var1,
				   align='center', label='GFP', color='#17becf')
		ax[0].barh(max_PC_diffs_ids, max_var0,
				   align='center', label='no GFP', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False,
						  labelleft=True, labelright=False)
		ax[0].set_xlim([0, max_val])
		for m, n in enumerate(max_var1):
			ax[0].text(n + .05, m - .1, max_change_str[m],
					   color='black', fontweight='bold', size=6)
		ax[0].set_xlabel("log(Average Protein Counts)")
		ax[0].set_ylabel("Protein ID")
		ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[0].set_title(f"The {shared_diff_num} proteins with the greatest "
						f"increase \n in protein count with the addition "
						f"of GFP", size=10)

		# minimum fold change plot specs:
		ax[1].barh(min_PC_diff_ids, min_var0,
				   align='center', label='no GFP', color='#bcbd22')
		ax[1].barh(min_PC_diff_ids, min_var1,
				   align='center', label='GFP', color='#17becf')
		ax[1].tick_params(bottom=True, top=False, left=False, right=True)
		ax[1].tick_params(axis='both', which='major', labelsize=6)
		ax[1].tick_params(labelbottom=True, labeltop=False,
						  labelleft=False, labelright=True)
		ax[1].set_xlim([0, max_val])
		for m, n in enumerate(min_var1):
			ax[1].text(n - .7, m - .1, min_change_str[m],
					   color='black', fontweight='bold', size=6)
		ax[1].set_xlabel("log(Average Protein Counts)")
		ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[1].set_title(f"The {shared_diff_num} proteins with the greatest "
						f"decrease \n in protein count with the addition "
						f"of GFP", size=10)

		fig.suptitle("Maximum Increases and Decreases in Protein Counts "
					 "\nwith the Addition of GFP")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' +
					 str(shared_diff_num) +
					 '_max_PC_diff_comparisons_logscale_Filter_' +
					 str(filter_num),
					 metadata)
		plt.close('all')



		# Observe which proteins have the greatest fold increase in
		# protein counts:

		# Obtain the proteins with the max fold change in protein counts
		# between variants:
		PC_max_folds = abs(nonzero_PCs[1] / nonzero_PCs[0])
		PC_max_fold_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		psortdiff = np.argsort(PC_max_folds)
		PC_max_fold_idxs = psortdiff[-max_fold_num:]
		max_fold_changes = []
		NZ_max_fold_PCs_ids = []
		for idx in PC_max_fold_idxs:
			max_fold_changes.append(PC_max_fold_percents[idx])
			NZ_max_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_max_fold_PCs = np.zeros((len(variants),
											len(PC_max_fold_idxs)))
		avg_NZ_max_fold_PCs = np.zeros((len(variants),
											len(PC_max_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_max_fold_idxs)):
				index = PC_max_fold_idxs[idx]
				avg_NZ_max_fold_PCs[variant][idx] = nonzero_PCs[variant][index]
				avg_log_NZ_max_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the results:
		#plt.figure(figsize=(8.5, 20))
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

		max_var0 = avg_log_NZ_max_fold_PCs[0]
		max_var1 = avg_log_NZ_max_fold_PCs[1]

		max_fold_val = [max(avg_log_NZ_max_fold_PCs[0]),
						max(avg_log_NZ_max_fold_PCs[1])]
		max_fold_val = max(max_fold_val)
		max_fold_val = round(max_fold_val) + 1

		# create labels for the % changes:
		max_change_str = []
		for m in range(max_fold_num):
			max_value = '+' + str(round(max_fold_changes[m], 2)) + '%'
			max_change_str.append(max_value)

		ax.barh(NZ_max_fold_PCs_ids, max_var0,
				 0.1, align='edge', label='no GFP')
		ax.barh(NZ_max_fold_PCs_ids, max_var1,
				 height=-0.1, align='edge', label='GFP')
		ax.set_xlim([0, max_fold_val])
		for m, n in enumerate(max_var1):
			ax.text(n + .01, m - .15, max_change_str[m],
					   color='black', fontweight='bold', size=10)
		ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
		ax.set_ylabel("Protein ID", fontweight='bold')
		ax.legend()
		ax.set_title(f"The {max_fold_num} proteins with the greatest fold"
				  f" \n increase in protein count between variants"
				  )

		plt.tight_layout()
		exportFigure(plt, plotOutDir,
					 plotOutFileName + '_' + str(max_fold_num) +
					 '_max_PC_fold_increases_woGFP_wPD_Filter_' +
					 str(filter_num),
					 metadata)
		plt.close('all')



		# Observe which proteins have the largest fold decrease in
		# protein counts:

		PC_min_folds = abs(nonzero_PCs[0] / nonzero_PCs[1])
		PC_min_fold_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		min_psortdiff = np.argsort(PC_min_folds)
		PC_min_fold_idxs = min_psortdiff[-min_fold_num:]
		min_changes = []
		NZ_min_fold_PCs_ids = []
		for idx in PC_min_fold_idxs:
			min_changes.append(PC_min_fold_percents[idx])
			NZ_min_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_min_fold_PCs = np.zeros((len(variants),
											len(PC_min_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_min_fold_idxs)):
				index = PC_min_fold_idxs[idx]
				avg_log_NZ_min_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the results:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 50))

		min_var0 = avg_log_NZ_min_fold_PCs[0]
		min_var1 = avg_log_NZ_min_fold_PCs[1]

		min_fold_val = [max(avg_log_NZ_min_fold_PCs[0]),
						max(avg_log_NZ_min_fold_PCs[1])]
		min_fold_val = max(min_fold_val)
		min_fold_val = round(min_fold_val) + 1

		# create labels for the % changes:
		min_change_str = []
		for m in range(min_fold_num):
			min_value = ' ' + str(round(min_changes[m], 2)) + '%'
			min_change_str.append(min_value)

		ax.barh(NZ_min_fold_PCs_ids, min_var0,
				 0.1, align='edge', label='no GFP')
		ax.barh(NZ_min_fold_PCs_ids, min_var1,
				 height=-0.1, align='edge', label='GFP')
		ax.set_xlim([0, min_fold_val])
		for m, n in enumerate(min_var1):
			ax.text(n + .01, m - .3, min_change_str[m],
					color='black', fontweight='bold', size=10)
		ax.set_xlabel("log(Average Protein Count)", fontweight='bold')
		ax.set_ylabel("Protein ID", fontweight='bold')
		ax.legend()
		ax.set_title(f"The {min_fold_num} proteins with the largest fold "
				  f" \ndecrease in protein counts between variants"
				  )

		plt.tight_layout()
		exportFigure(plt, plotOutDir,
					 plotOutFileName + '_' + str(min_fold_num) +
					 '_max_PC_fold_decreases_woGFP_wPD_Filter_' +
					 str(filter_num),
					 metadata)
		plt.close('all')


		# Observe the max and min fold changes in protein counts side by side:

		# Obtain the proteins with the max fold change in protein counts
		# between variants:
		PC_max_folds = abs(nonzero_PCs[1] / nonzero_PCs[0])
		PC_max_fold_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		psortdiff = np.argsort(PC_max_folds)
		PC_max_fold_idxs = psortdiff[-sharednum:]
		max_changes = []
		NZ_max_fold_PCs_ids = []
		for idx in PC_max_fold_idxs:
			max_changes.append(PC_max_fold_percents[idx])
			NZ_max_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_max_fold_PCs = np.zeros((len(variants),
											len(PC_max_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_max_fold_idxs)):
				index = PC_max_fold_idxs[idx]
				avg_log_NZ_max_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Obtain the proteins with the min fold change in protein counts
		# between variants:
		PC_min_folds = abs(nonzero_PCs[0] / nonzero_PCs[1])
		PC_min_fold_percents = (100 * (nonzero_PCs[1] - nonzero_PCs[0])
								/ nonzero_PCs[0])
		min_psortdiff = np.argsort(PC_min_folds)
		PC_min_fold_idxs = min_psortdiff[-sharednum:]
		min_changes = []
		NZ_min_fold_PCs_ids = []
		for idx in PC_min_fold_idxs:
			min_changes.append(PC_min_fold_percents[idx])
			NZ_min_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_min_fold_PCs = np.zeros((len(variants),
											len(PC_min_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_min_fold_idxs)):
				index = PC_min_fold_idxs[idx]
				avg_log_NZ_min_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the max and min fold changes together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))

		# Name the results:
		max_var0 = avg_log_NZ_max_fold_PCs[0]
		max_var1 = avg_log_NZ_max_fold_PCs[1]
		min_var0 = avg_log_NZ_min_fold_PCs[0]
		min_var1 = avg_log_NZ_min_fold_PCs[1]

		# create labels for the % changes:
		max_change_str = []
		min_change_str = []
		for m in range(sharednum):
			max_val = '+' + str(round(max_changes[m],2)) + '%'
			min_val = ' ' + str(round(min_changes[m], 2)) + '%'
			max_change_str.append(max_val)
			min_change_str.append(min_val)

		# maximum fold change plot specs:
		ax[0].barh(NZ_max_fold_PCs_ids, max_var1, align='center',
				   label='GFP', color='#17becf')
		ax[0].barh(NZ_max_fold_PCs_ids, max_var0, align='center',
				   label='no GFP', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False,
						  labelleft=True, labelright=False)
		ax[0].set_xlim([0, 3])
		for m, n in enumerate(max_var1):
			ax[0].text(n + .05, m - .1, max_change_str[m], color='black',
					   fontweight='bold', size=6)
		ax[0].set_xlabel("log(Average Protein Counts)", fontweight='bold')
		ax[0].set_ylabel("Protein ID", fontweight='bold')
		#ax[0].legend(loc='best')
		ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[0].set_title(f"The {sharednum} proteins with the greatest fold "
						f"increase \n in protein counts with the"
						f" addition of GFP", size=10)

		# minimum fold change plot specs:
		ax[1].barh(NZ_min_fold_PCs_ids, min_var0, align='center',
				   label='no GFP', color='#bcbd22')
		ax[1].barh(NZ_min_fold_PCs_ids, min_var1, align='center',
				   label='GFP', color='#17becf')
		ax[1].tick_params(bottom=True, top=False, left=False, right=True)
		ax[1].tick_params(axis='both', which='major', labelsize=6)
		ax[1].tick_params(labelbottom=True, labeltop=False,
						  labelleft=False, labelright=True)
		ax[1].set_xlim([0, 3])
		for m, n in enumerate(min_var1):
			ax[1].text(n + .05, m - .1, min_change_str[m], color='black',
					   fontweight='bold', size=6)
		ax[1].set_xlabel("log(Average Protein Counts)", fontweight='bold')
		#ax[1].legend(loc='best')
		ax[1].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=2,
					 fontsize='small')
		ax[1].set_title(f"The {sharednum} proteins with the greatest fold "
						f"decrease \n in protein counts with the addition of"
						f" GFP", size=10)

		fig.suptitle("Maximum and Minimum Fold Change Comparisons for"
					 " Protein Counts "
					 "\nwith the Addition of GFP")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(sharednum) +
					 '_max_PC_fold_comparisons_Filter_' + str(filter_num),
					 metadata)
		plt.close('all')



if __name__ == "__main__":
	Plot().cli()
