"""
Plot mRNA and protein counts for genes across multiple generations
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Replace with the proteins you would like to visualize here:
interest_proteins = np.array([
	 #'ACRD-MONOMER[i]',
	 #'CYNX-MONOMER[i]',
	 #'B0270-MONOMER[i]',
	 #'G7634-MONOMER[i]',
	#'EG11854-MONOMER[c]',
	#'G6606-MONOMER[c]',
	#'MONOMER0-2678[c]',
	#'EG10037-MONOMER[c]',
	#'NG-GFP-MONOMER[c]',
	'TRYPSYN-APROTEIN[c]',
	'TRYPSYN-BPROTEIN[c]',
	#"ANTHRANSYNCOMPI-MONOMER[c]",
	#'PD00519[c]',
])

# if you know how many variants you have and want to specify which you would
# like to plot specifically, you can do so here. Otherwise, leave this as 0.
experimental_vars = 0
# example of specific variants (used in sherlock):
#experimental_vars = [16, 17, 18, 19, 20]

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def seed_plot(self, control_var, experimental_vars):
		# retrieve the seed data
		seeds = self.ap.get_seeds()
		# Plotting
		colors = [["turquoise", "yellowgreen", "mediumpurple", "deeppink", "darkturquoise"],
				  ["deepskyblue","lightcoral", "gold", "darkorange", "indigo"], ["darkred",
				  "darkgreen", "darkblue", "darkviolet", "cornflowerblue"],]
		ccolors = ["#FF796C", "slateblue", "darkviolet", "plum", "sandybrown",]
		# Protein Counts
		fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(11, 16))
		LS = ['-', ':', '-.', '--']
		for seed in range(len(seeds)):
			cseed_dir = self.ap.get_cells(variant=[control_var], seed=[seed],
										 only_successful=True)
			ctime = read_stacked_columns(cseed_dir, 'Main', 'time',
										ignore_exception=True)
			(cip_monomer_counts,) = read_stacked_bulk_molecules(
				cseed_dir, self.cistron_monomer_ids, ignore_exception=True)
			if len(self.cistron_monomer_ids) == 1:
				name = 'seed ' + str(seed) + ' cntrl var '
				ax[0].plot(ctime / 60., cip_monomer_counts,
						 linestyle=LS[seed],
						 label=name, color="#FF796C", linewidth=.5)
				var_num = 0
				for variant in experimental_vars:
					seed_dir = self.ap.get_cells(variant=[variant], seed=[seed],
												 only_successful=True)
					# Load data for the seed
					time = read_stacked_columns(seed_dir, 'Main', 'time',
												ignore_exception=True)
					(ip_monomer_counts,) = read_stacked_bulk_molecules(
						seed_dir, self.cistron_monomer_ids, ignore_exception=True)
					name = 'seed ' + str(seed) + ' exp  var ' + str(variant)
					ax[0].plot(time / 60., ip_monomer_counts, label=name,
							 color=colors[0][var_num], linestyle=LS[seed], linewidth=.5)
					var_num = var_num + 1
					print(var_num)
				# plot specs
				ax[0].set_title(f"Protein Counts for {self.cistron_monomer_ids[0]} in "
						  f"\nthe control variant"
						  f" and {len(experimental_vars)} experimental variant(s)")

		# finish protein count plot
		ax[0].set_xlabel("Time (min)")
		ax[0].set_ylabel("Protein Counts")
		ax[0].legend(loc='lower center', bbox_to_anchor=(0.5, -.2), ncols=3,
					 fontsize='small')

		# mRNA Counts
		for seed in range(len(seeds)):
			cseed_dir = self.ap.get_cells(variant=[control_var], seed=[seed],
										  only_successful=True)
			ctime = read_stacked_columns(cseed_dir, 'Main', 'time',
										 ignore_exception=True)
			cip_mRNA_counts = read_stacked_columns(
				cseed_dir, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True)[:, self.new_gene_mRNA_indexes]
			if len(self.cistron_ids) == 1:
				name = 'seed ' + str(seed) + ' cntrl var '
				ax[1].plot(ctime / 60., cip_mRNA_counts, linestyle=LS[seed],
						 label=name, color="#FF796C", linewidth=.5)
				var_num = 0
				for variant in experimental_vars:
					seed_dir = self.ap.get_cells(variant=[variant], seed=[seed],
												 only_successful=True)
					# Load data for the seed
					time = read_stacked_columns(seed_dir, 'Main', 'time',
												ignore_exception=True)
					ip_mRNA_counts = read_stacked_columns(
						seed_dir, 'RNACounts', 'mRNA_cistron_counts',
						ignore_exception=True)[:, self.new_gene_mRNA_indexes]
					name = 'seed ' + str(seed) + ' exp  var ' + str(variant)
					ax[1].plot(time / 60., ip_mRNA_counts, label=name,
							 color=colors[0][var_num], linestyle=LS[seed], linewidth=.5)
					var_num = var_num + 1
				# plot specs
				ax[1].set_title(f"mRNA Counts for {self.cistron_ids[0]} in the \n "
						  f"control variant"
						  f" and {len(experimental_vars)} experimental variant(s)")

		ax[1].set_xlabel("Time (min)")
		ax[1].set_ylabel("Cistron Counts")
		ax[1].legend(loc='upper center', bbox_to_anchor=(.5, 1.2), ncols=3,
					 fontsize='small')

		plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
		plt.tight_layout()

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		# extract shared information
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_cistron_sim_data = (
			sim_data.process.transcription.cistron_data.struct_array)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)

		# extract info about the protein(s) from the monomer data:
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			ip_monomer_data = monomer_sim_data[monomer_idx]
			ip_cistron_id = ip_monomer_data['cistron_id']
			self.cistron_ids = [ip_cistron_id]
			cistron_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
											   monomer_sim_data['id']))
			self.cistron_monomer_ids = [cistron_monomer_id_dict.get(mRNA_id)
										for mRNA_id in self.cistron_ids]

			# get the cell paths for the experimental variant
			cell_paths = self.ap.get_cells()
			sim_dir = cell_paths[0]
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Extract mRNA indexes for each gene of interest
			mRNA_counts_reader = TableReader(os.path.join(simOutDir,
														  'RNACounts'))
			mRNA_idx_dict = {rna: i for i, rna in enumerate(
				mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
			self.new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
									 self.cistron_ids]

			# get the data for all variants:
			all_variants = self.ap.get_variants()
			# specifiy the control and experimental variants:
			control_var = all_variants[0]
			# if the experimental vars list is empty, plot all other variants
			if experimental_vars == 0:
				experimental_variants = all_variants[1:]

			# plot the data:
			self.seed_plot(control_var, experimental_variants)
			plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + str(protein)
						 + '_all_seeds', metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
