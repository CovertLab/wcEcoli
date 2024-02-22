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

		# Randomly select genes to observe with the new gene(s):
		# edit this number as desired:
		randnum = 100

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

			for i in range(len(protein_idxs)):
				each_gen_avg_monomer_counts[i][variant] \
					= gen_avg_monomer_counts[:,i]
				each_gen_avg_log_monomer_counts[i][variant] = \
					np.log10(gen_avg_monomer_counts[:,i] + 1)
				all_avg_monomer_counts[i][variant] = avg_monmer_counts[i]
				all_avg_log_monomer_counts[i][variant] = \
					np.log10(avg_monmer_counts[i] + 1)

		# Make plot for the random protein selections:
		plt.figure(figsize = (8.5, 20))

		logdata = all_avg_log_monomer_counts
		ld=np.zeros((len(variants), len(protein_idxs)))

		for variant in variants:
			vardata=np.zeros((len(protein_idxs)))
			for i in range(len(protein_idxs)):
				data = logdata[i][variant]
				vardata[i] = data
			ld[variant] = vardata

		var0 = ld[0]
		var1 = ld[1]

		plt.barh(monomer_names, var0, 0.1, align = 'edge', label='no GFP')
		plt.barh(monomer_names, var1, -0.1, align = 'edge', label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title("Protein Count Comparisons "
				  "\n for "+str(randnum)
				  +" randomly selected proteins ")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_Protein_Counts_for_'
					 +str(randnum) +'_random_proteins', metadata)
		plt.close('all')


		# Observe the protein count size differences between the two variants:
		# TODO: change this number as desired:
		num = 22

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
		# differences between variants:
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
			for i in range(len(interest_protein_idxs)):
				avg_log_interest_proteins[variant][i] = \
					np.log10(interest_protein_counts[variant][i] + 1)

		# Make plot for genes with the largest protein count differences:
		plt.figure(figsize = (8.5, 20))
		var0 = avg_log_interest_proteins[0]
		var1 = avg_log_interest_proteins[1]

		plt.barh(interest_proteins, var0, 0.1, align='edge', label='no GFP')
		plt.barh(interest_proteins, var1, -0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {num} proteins with the greatest change in protein "
				  f"\ncount between two varaints when a new gene is added")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName+'_max_' +
					 str(num) +
					 '_protein_count_differences_with_GFP',
					 metadata)
		plt.close('all')


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
		plt.figure(figsize=(8.5, 70))
		var0 = avg_log_NS_0_PCs[0]
		var1 = avg_log_NS_0_PCs[1]

		plt.barh(NS_0_PC_ids, var0, 0.1, align='edge', label='no GFP')
		plt.barh(NS_0_PC_ids, var1, -0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {len(nonshared_0_PC_ids)} proteins counts for proteins"
				  f" that appear in only one variant"
				  )

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' +
					 str(len(nonshared_0_PC_ids)) +
					 '_unique_protein_appearances_PCs', metadata)
		plt.close('all')



		# Extract all proteins with non-zero protein counts in both variants:
		nonzero_p_counts_var0_idxs = np.nonzero(p_counts[0])
		nonzero_p_counts_var1_idxs = np.nonzero(p_counts[1])
		shared_nonzero_PCs_idxs = np.intersect1d(nonzero_p_counts_var0_idxs,
												 nonzero_p_counts_var1_idxs)
		nonzero_PCs = p_counts[:, shared_nonzero_PCs_idxs]
		nonzero_PCs_ids = [inv_monomer_idx_dict.get(monomer_id)
						   for monomer_id in shared_nonzero_PCs_idxs]


		# Observe which protiens have the greatest increase in
		# protein counts with the addition of GFP:
		# TODO: Change this number as desired:
		min_num = 21

		min_diffs = abs(nonzero_PCs[1] - nonzero_PCs[0])
		sortdiff = np.argsort(min_diffs)
		min_diff_idxs = sortdiff[:min_num]
		min_changes = []
		min_PC_diffs_ids = []
		for idx in min_diff_idxs:
			min_changes.append(min_diffs[idx])
			min_PC_diffs_ids.append(nonzero_PCs_ids[idx])
		avg_log_min_PC_diffs = np.zeros((len(variants), len(min_diff_idxs)))
		for variant in variants:
			for idx in range(len(min_diff_idxs)):
				index = min_diff_idxs[idx]
				avg_log_min_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Make plot for genes with the largest protein count differences:
		plt.figure(figsize=(8.5, 20))
		min_diff_var0 = avg_log_min_PC_diffs[0]
		min_diff_var1 = avg_log_min_PC_diffs[1]

		plt.barh(min_PC_diffs_ids, min_diff_var0,
				 0.1, align='edge', label='no GFP')
		plt.barh(min_PC_diffs_ids, min_diff_var1,
				 height=-0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {min_num} proteins with the smallest change in protein"
				  f" \ncount between two varaints when a new gene is added")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(min_num) +
					 '_min_protein_count_change', metadata)
		plt.close('all')

		# Observe which protiens have the greatest increase in protein counts
		# with the addition of GFP (but not including GFP):
		# TODO: Change this number as desired:
		max_num = 20

		max_diffs = abs(nonzero_PCs[1] - nonzero_PCs[0])
		sortdiff = np.argsort(max_diffs)
		max_diff_idxs = sortdiff[-max_num:]
		max_changes = []
		max_PC_diffs_ids = []
		for idx in max_diff_idxs:
			max_changes.append(max_diffs[idx])
			max_PC_diffs_ids.append(nonzero_PCs_ids[idx])
		avg_log_max_PC_diffs = np.zeros((len(variants), len(max_diff_idxs)))
		for variant in variants:
			for idx in range(len(max_diff_idxs)):
				index = max_diff_idxs[idx]
				avg_log_max_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Make plot for genes with the largest protein count differences:
		plt.figure(figsize=(8.5, 20))
		max_diff_var0 = avg_log_max_PC_diffs[0]
		max_diff_var1 = avg_log_max_PC_diffs[1]

		plt.barh(max_PC_diffs_ids, max_diff_var0,
				 0.1, align='edge', label='no GFP')
		plt.barh(max_PC_diffs_ids, max_diff_var1,
				 height=-0.1, align='edge', label='GFP')

		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {max_num} proteins with the greatest increase in protein"
				  f" \ncounts between two varaints when GFP is added")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(max_num) +
					 '_max_protein_count_change', metadata)
		plt.close('all')


		# Observe the max and min fold changes in protein counts side by side:
		# TODO: Change the number below to as desired:
		shared_diff_num = 11

		min_PC_diffs = abs(nonzero_PCs[1] - nonzero_PCs[0])
		sortdiff = np.argsort(min_PC_diffs)
		min_diff_idxs = sortdiff[:shared_diff_num]
		min_PC_diff_changes = []
		min_PC_diff_ids = []
		for idx in min_diff_idxs:
			min_PC_diff_changes.append(min_PC_diffs[idx])
			min_PC_diff_ids.append(nonzero_PCs_ids[idx])
		avg_log_min_PC_diffs = np.zeros((len(variants), len(min_diff_idxs)))
		for variant in variants:
			for idx in range(len(min_diff_idxs)):
				index = min_diff_idxs[idx]
				avg_log_min_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		max_PC_diffs = abs(nonzero_PCs[1] - nonzero_PCs[0])
		sortdiff = np.argsort(max_PC_diffs)
		max_diff_idxs = sortdiff[-shared_diff_num:]
		max_PC_diff_changes = []
		max_PC_diffs_ids = []
		for idx in max_diff_idxs:
			max_PC_diff_changes.append(max_PC_diffs[idx])
			max_PC_diffs_ids.append(nonzero_PCs_ids[idx])
		avg_log_max_PC_diffs = np.zeros((len(variants), len(max_diff_idxs)))
		for variant in variants:
			for idx in range(len(max_diff_idxs)):
				index = max_diff_idxs[idx]
				avg_log_max_PC_diffs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the max and min fold changes together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

		# Name the results:
		max_var0 = avg_log_max_PC_diffs[0]
		max_var1 = avg_log_max_PC_diffs[1]
		min_var0 = avg_log_min_PC_diffs[0]
		min_var1 = avg_log_min_PC_diffs[1]

		# create labels for the % changes:
		max_change_str = []
		min_change_str = []
		for m in range(shared_diff_num):
			max_val = '+' + str(round(max_PC_diff_changes[m], 2)) + '%'
			min_val = '+' + str(round(min_PC_diff_changes[m], 2)) + '%'
			max_change_str.append(max_val)
			min_change_str.append(min_val)

		# maximum fold change plot specs:
		ax[0].barh(max_PC_diffs_ids, max_var1,
				   align='center', label='GFP', color='#17becf')
		ax[0].barh(max_PC_diffs_ids, max_var0,
				   align='center', label='no GFP', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False,
						  labelleft=True, labelright=False)
		ax[0].set_xlim([0, 3])
		for i, v in enumerate(max_var1):
			ax[0].text(v + .05, i - .1, max_change_str[i],
					   color='black', fontweight='bold', size=6)
		ax[0].set_xlabel("log(Average Protein Counts)")
		ax[0].set_ylabel("Protein ID")
		ax[0].legend(loc='best')
		ax[0].set_title(f"The {shared_diff_num} proteins with the greatest "
						f"increase \n in protein counts with the addition "
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
		ax[1].set_xlim([0, 3])
		for i, v in enumerate(min_var1):
			ax[1].text(v + .05, i - .1, min_change_str[i],
					   color='black', fontweight='bold', size=6)
		ax[1].set_xlabel("log(Average Protein Counts)")
		ax[1].legend(loc='best')
		ax[1].set_title(f"The {shared_diff_num} proteins with the smallest "
						f"increase \n in protein counts with the addition "
						f"of GFP", size=10)

		fig.suptitle("Maximum and Minimum Protein Count Change Comparisons "
					 "\nwith the Addition of GFP")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' +
					 str(shared_diff_num) +
					 '_max_and_min_PC_difference_comparisons',
					 metadata)
		plt.close('all')


		# Observe which proteins have the greatest fold increase in
		# protein counts:
		# Change this number as desired:
		max_fold_num = 12

		# Obtain the proteins with the max fold change in protein counts
		# between variants:
		PC_max_folds = abs(nonzero_PCs[1] / nonzero_PCs[0])
		psortdiff = np.argsort(PC_max_folds)
		PC_max_fold_idxs = psortdiff[-max_fold_num:]
		max_changes = []
		NZ_max_fold_PCs_ids = []
		for idx in PC_max_fold_idxs:
			max_changes.append(PC_max_folds[idx])
			NZ_max_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_max_fold_PCs = np.zeros((len(variants),
											len(PC_max_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_max_fold_idxs)):
				index = PC_max_fold_idxs[idx]
				avg_log_NZ_max_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the results:
		plt.figure(figsize=(8.5, 20))
		max_var0 = avg_log_NZ_max_fold_PCs[0]
		max_var1 = avg_log_NZ_max_fold_PCs[1]

		plt.barh(NZ_max_fold_PCs_ids, max_var0,
				 0.1, align='edge', label='no GFP')
		plt.barh(NZ_max_fold_PCs_ids, max_var1,
				 height=-0.1, align='edge', label='GFP')
		#ax.set_yticklabels(max_changes)
		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {max_fold_num} proteins with the greatest fold"
				  f" \n increase in protein counts between variants"
				  )

		plt.tight_layout()
		exportFigure(plt, plotOutDir,
					 plotOutFileName + '_' + str(max_fold_num) +
					 '_max_fold_protein_count_changes',
					 metadata)
		plt.close('all')

		# Observe which proteins have the largest fold decrease in protein counts:
		# TODO: Change this number as desired:
		min_fold_num = 13

		PC_min_folds = abs(nonzero_PCs[0] / nonzero_PCs[1])
		min_psortdiff = np.argsort(PC_min_folds)
		PC_min_fold_idxs = min_psortdiff[-min_fold_num:]
		min_changes = []
		NZ_min_fold_PCs_ids = []
		for idx in PC_min_fold_idxs:
			min_changes.append(PC_min_folds[idx])
			NZ_min_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_min_fold_PCs = np.zeros((len(variants),
											len(PC_min_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_min_fold_idxs)):
				index = PC_min_fold_idxs[idx]
				avg_log_NZ_min_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the results:
		plt.figure(figsize=(8.5, 20))
		min_var0 = avg_log_NZ_min_fold_PCs[0]
		min_var1 = avg_log_NZ_min_fold_PCs[1]
		plt.barh(NZ_min_fold_PCs_ids, min_var0,
				 0.1, align='edge', label='no GFP')
		plt.barh(NZ_min_fold_PCs_ids, min_var1,
				 height=-0.1, align='edge', label='GFP')
		plt.xlabel("log(Average Protein Count)")
		plt.ylabel("Protein ID")
		plt.legend()
		plt.title(f"The {min_fold_num} proteins with the largest fold "
				  f" \ndecrease in protein counts between variants"
				  )

		plt.tight_layout()
		exportFigure(plt, plotOutDir,
					 plotOutFileName + '_' + str(min_fold_num) +
					 '_min_fold_protein_count_changes',
					 metadata)
		plt.close('all')


		# Observe the max and min fold changes in protein counts side by side:
		# TODO: Change the number below to as desired:
		sharednum = 100

		# Obtain the proteins with the max fold change in protein counts
		# between variants:
		PC_max_folds = abs(nonzero_PCs[1] / nonzero_PCs[0])
		psortdiff = np.argsort(PC_max_folds)
		PC_max_fold_idxs = psortdiff[-sharednum:]
		max_changes = []
		NZ_max_fold_PCs_ids = []
		for idx in PC_max_fold_idxs:
			max_changes.append(PC_max_folds[idx])
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
		min_psortdiff = np.argsort(PC_min_folds)
		PC_min_fold_idxs = min_psortdiff[-sharednum:]
		min_changes = []
		NZ_min_fold_PCs_ids = []
		for idx in PC_min_fold_idxs:
			min_changes.append(PC_min_folds[idx])
			NZ_min_fold_PCs_ids.append(nonzero_PCs_ids[idx])
		avg_log_NZ_min_fold_PCs = np.zeros((len(variants), len(PC_min_fold_idxs)))
		for variant in variants:
			for idx in range(len(PC_min_fold_idxs)):
				index = PC_min_fold_idxs[idx]
				avg_log_NZ_min_fold_PCs[variant][idx] = \
					np.log10(nonzero_PCs[variant][index] + 1)

		# Plot the max and min fold changes together:
		fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 50))

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
			min_val = '-' + str(round(min_changes[m], 2)) + '%'
			max_change_str.append(max_val)
			min_change_str.append(min_val)

		# maximum fold change plot specs:
		ax[0].barh(NZ_max_fold_PCs_ids, max_var1, align='center', label='GFP', color='#17becf')
		ax[0].barh(NZ_max_fold_PCs_ids, max_var0, align='center', label='no GFP', color='#bcbd22')
		ax[0].tick_params(axis='both', which='major', labelsize=6)
		ax[0].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
		ax[0].set_xlim([0, 3])
		for i, v in enumerate(max_var1):
			ax[0].text(v + .05, i - .1, max_change_str[i], color='black', fontweight='bold', size=6)
		ax[0].set_xlabel("log(Average Protein Counts)")
		ax[0].set_ylabel("Protein ID")
		ax[0].legend(loc='best')
		ax[0].set_title(f"The {sharednum} proteins with the greatest fold "
						f"increase \n in protein counts with the addition of GFP", size=10)

		# minimum fold change plot specs:
		ax[1].barh(NZ_min_fold_PCs_ids, min_var0, align='center', label='no GFP', color='#bcbd22')
		ax[1].barh(NZ_min_fold_PCs_ids, min_var1, align='center', label='GFP', color='#17becf')
		ax[1].tick_params(bottom=True, top=False, left=False, right=True)
		ax[1].tick_params(axis='both', which='major', labelsize=6)
		ax[1].tick_params(labelbottom=True, labeltop=False, labelleft=False, labelright=True)
		ax[1].set_xlim([0, 3])
		for i, v in enumerate(min_var1):
			ax[1].text(v + .05, i - .1, min_change_str[i], color='black', fontweight='bold', size=6)
		ax[1].set_xlabel("log(Average Protein Counts)")
		ax[1].legend(loc='best')
		ax[1].set_title(f"The {sharednum} proteins with the greatest fold "
						f"decrease \n in protein counts with the addition of GFP", size=10)

		fig.suptitle("Maximum and Minimum Fold Change Comparisons for Protein Counts "
					 "\nwith the Addition of GFP")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + '_' + str(sharednum) +
					 '_max_and_min_fold_comparisons', metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
