"""
Plot mRNA and protein counts for genes across multiple seeds and generations

NOTE: it is tough to do this plot with out splitting up the seeds because the
doubling times are different across different seeds (so they will not line up
nicely)
"""

# todo: add descriptions for each function, add averaging functionality, figure out what seedOutDir must be to handle seeds that start at any seed, check that plot out names are ok with the code style

import pickle
import os
from matplotlib import pyplot as plt
import numpy as np
import io
from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH
from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Replace with the proteins you would like to visualize here:
interest_proteins = np.array([
	# 'ACRD-MONOMER[i]',
	# 'CYNX-MONOMER[i]',
	# 'B0270-MONOMER[i]',
	# 'G7634-MONOMER[i]',
	#'EG11854-MONOMER[c]',
	#'G6606-MONOMER[c]',
	#'EG10542-MONOMER[c]',
	'PD03938[c]', # metR
	'RPOS-MONOMER[c]', # rpoS
])

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# function to extract data from the simOut directory
		def extract_data(simOutDir, cistron_ids, cell_paths, cistron_monomer_ids):
			# Extract mRNA indexes for each gene of interest
			mRNA_counts_reader = TableReader(os.path.join(simOutDir,
														  'RNACounts'))
			mRNA_idx_dict = {rna: i for i, rna in enumerate(
				mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
			gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
									 cistron_ids]

			# Load data
			time = read_stacked_columns(cell_paths, 'Main',
										'time', ignore_exception=True)
			(ip_monomer_counts,) = read_stacked_bulk_molecules(
				cell_paths, cistron_monomer_ids, ignore_exception=True)
			ip_mRNA_counts = read_stacked_columns(
				cell_paths, 'RNACounts', 'mRNA_cistron_counts',
				ignore_exception=True)[:, gene_mRNA_indexes]

			return time, ip_monomer_counts, ip_mRNA_counts


		# function to plot the counts per seed:
		def plot_counts_per_seed(seed, cell_paths, cistron_monomer_ids, time,
								 ip_monomer_counts, ip_mRNA_counts,):

			# Plotting
			plt.figure(figsize=(8.5, 11))
			# Get doubling times from cells with this variant index
			dt = read_stacked_columns(
				cell_paths, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			dts = np.zeros(len(dt))
			for i in range(len(dt)):
				if i == 0:
					gt = dt[i]
					dts[i] = gt
				else:
					gt = dt[i] + dts[i - 1]
					dts[i] = gt

			# Protein Counts
			plt.subplot(2, 1, 1)
			for x in dts:
				plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)

			gene_list = []
			if len(cistron_monomer_ids) == 1:
				common_name = get_common_names(cistron_monomer_ids[0])
				name = cistron_monomer_ids[0] + ' (' + common_name + ')'
				gene_list.append(common_name)
				plt.plot(time / 60., ip_monomer_counts,
						 label=name)
			else:
				for m in range(len(cistron_monomer_ids)):
					common_name = get_common_names(cistron_monomer_ids[m])
					name = cistron_monomer_ids[m] + ' (' + common_name + ')'
					gene_list.append(common_name)
					plt.plot(time / 60., ip_monomer_counts[:, m],
							 label=name, alpha=0.5)

			plt.xlabel("Time (min)"); plt.ylabel("Free Monomer Counts")
			plt.title(f"Free Monomer Counts in Seed {seed}")
			plt.legend()

			# mRNA Counts
			plt.subplot(2, 1, 2)
			for x in dts:
				plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
			if len(cistron_ids) == 1:

				plt.plot(time / 60., ip_mRNA_counts,
						 label=cistron_ids[0])
			else:
				for r in range(len(cistron_ids)):
					common_name = get_common_names(cistron_monomer_ids[r])
					name = cistron_ids[r] + ' (' + common_name + ')'
					plt.plot(time / 60., ip_mRNA_counts[:, r],
							 label=name, alpha=0.5)
			plt.xlabel("Time (min)")
			plt.ylabel("Cistron Counts")
			plt.title(f"mRNA Counts for Proteins of Interest in Seed {seed}")
			plt.legend()

			plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
			plotOutPath = os.path.join(plotOutDir, f'00000{seed}')
			exportFigure(plt, plotOutPath, plotOutFileName + '_cohortPlot_seed_'
						 + str(seed) + '_geneIDs_' + str(gene_list), metadata)
			plt.close("all")

		def get_gene_symbols_for_monomer_ids():
			# code adapted from convert_to_flat.py
			RNAS_FILE = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli', 'flat', 'rnas.tsv')
			with io.open(RNAS_FILE, 'rb') as f:
				reader = tsv.reader(f, delimiter='\t')
				headers = next(reader)
				while headers[0].startswith('#'):
					headers = next(reader)

				gene_symbol_index = headers.index('common_name')
				protein_id_index = headers.index('monomer_ids')
				monomer_ids_to_gene_symbols = {}
				for line in reader:
					gene_symbol = line[gene_symbol_index]
					protein_id = list(line[protein_id_index][2:-2].split('", "'))[0]
					monomer_ids_to_gene_symbols[protein_id] = gene_symbol

				return monomer_ids_to_gene_symbols

		# function to obtain common names for the proteins of interest
		def get_common_names(protein_id):
			"""
			Get the common names for each protein id
			Args:
				protein_id: the id for the proteins of interest
			Returns: the common name for the protein id
			"""
			protein = protein_id[:-3]
			gene_symbol = get_gene_symbols_for_monomer_ids()[protein]
			return gene_symbol


		# make a plot for the counts of each seed:
		def plot_counts_per_protein(simOutDir, protein, cistron,):
			# get the common name for the protein
			print(protein)
			common_name = get_common_names(protein[0])

			# Plotting
			plt.figure(figsize=(8.5, 11))

			# Protein Counts
			plt.subplot(2, 1, 1)
			# loop through the seeds:
			for seed in seeds:
				# cell paths for the seed:
				cell_paths = (
					self.ap.get_cells(seed=[seed]))
				# Extract data for the seed
				time, ip_monomer_counts, ip_mRNA_counts = (
					extract_data(simOutDir, cistron, cell_paths, protein))
				# plot the data
				plt.plot(time / 60., ip_monomer_counts,
						 label=f'seed {seed}', alpha=0.5)

			plt.xlabel("Time (min)")
			plt.ylabel("Free Monomer Counts")
			plt.title(f"Free Monomer Counts for {protein[0]} ({common_name})")
			plt.legend()

			# mRNA Counts
			plt.subplot(2, 1, 2)
			# loop through the seeds:
			for seed in seeds:
				# cell paths for the seed:
				cell_paths = (
					self.ap.get_cells(seed=[seed]))
				# Extract data for the seed
				time, ip_monomer_counts, ip_mRNA_counts = (
					extract_data(simOutDir, cistron, cell_paths, protein))
				# plot the data
				plt.plot(time / 60., ip_mRNA_counts,
						 label=f'seed {seed}', alpha=0.5)

			plt.xlabel("Time (min)")
			plt.ylabel("Cistron Counts")
			plt.title(f"mRNA Counts for {cistron[0]} (common name: {common_name})")
			plt.legend()

			# export plot:
			plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
			exportFigure(plt, plotOutDir, plotOutFileName +
						 '_cohortPlot_geneID_' + common_name , metadata)
			plt.close("all")



		# todo: use unedited plotOutDir to save plots in the wildtype folder
		#  when doing the singular protein with the average counts


		# extract data paths
		cell_paths = self.ap.get_cells()
		seed_path = self.ap.get_cells(seed=[1])
		print(seed_path)
		seeds = self.ap.get_seeds()
		print(seeds)
		print(plotOutDir)

		print(cell_paths)
		#import ipdb
		#ipdb.set_trace()
		sim_dir = cell_paths[
			0]  # /Users/miagrahn/wcEcoli/out/CLClim3dNE/wildtype_000000/000002/generation_000000/000000
		print(sim_dir)
		simOutDir = os.path.join(sim_dir, 'simOut') # does not matter what this is,
		# just need to get the simOut directory to find the indexes (which should be the same for each)

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_cistron_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array

		# extract info about the protein(s) from the monomer data:
		monomer_data_idxs = []
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			monomer_data_idxs.append(monomer_idx)
		ip_monomer_data = monomer_sim_data[monomer_data_idxs]
		ip_monomer_ids = ip_monomer_data['id']
		ip_cistron_ids = ip_monomer_data['cistron_id']
		print('ip_cistron_ids:' ); print(ip_cistron_ids)

		# extract info about the protein(s) from the mRNA/cistron data:
		mRNA_data_idxs = []
		for cistron in ip_cistron_ids:
			mRNA_idx = np.where(mRNA_cistron_sim_data['id'] == cistron)
			mRNA_idx = mRNA_idx[0][0]
			mRNA_data_idxs.append(mRNA_idx)
		ip_mRNA_data = mRNA_cistron_sim_data[mRNA_data_idxs]
		# todo: check if this section is even nessessary, bc I am not sure any of it is used
		ip_gene_ids = np.array(ip_mRNA_data['gene_id'])
		cistron_ids = ip_cistron_ids
		print('cistron_ids: ');
		print(cistron_ids)


		cistron_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										   monomer_sim_data['id']))
		cistron_monomer_ids = [cistron_monomer_id_dict.get(mRNA_id)
							   for mRNA_id in cistron_ids]


		# iterate through the seeds to generate graphs for each seed and
		# consolidates average counts and all on one graph
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed])) # seed needs to be in brackets

			# Extract data for the seed
			time, ip_monomer_counts, ip_mRNA_counts = (
				extract_data(simOutDir, cistron_ids, cell_paths,
							 cistron_monomer_ids))

			# Plot the counts for the seed
			plot_counts_per_seed(seed, cell_paths, cistron_monomer_ids, time,
								 ip_monomer_counts, ip_mRNA_counts)


		# Now, plot for each individual protein across all seeds
		for protein in interest_proteins:

			# get the the cistron id for the protein
			cistron_id = \
				[key for key, v in cistron_monomer_id_dict.items() if v == protein]
			protein_id = [protein]  # fomrat the protein id for extract_data
			print("cistron_id: ")
			print(cistron_id)

			# plot the counts for the protein
			print("plotting for protein: ")
			print(protein)
			print("cistron_monomer_ids: ")
			print(cistron_monomer_ids)
			print("cistron_ids: ")
			print(cistron_ids)
			plot_counts_per_protein(simOutDir, protein_id, cistron_id)


		# todo: figure out the difference between cistron ids and cistron_monomer_ids














		# # Extract mRNA indexes for each gene of interest
		# mRNA_counts_reader = TableReader(os.path.join(simOutDir,
		# 											  'RNACounts'))
		# mRNA_idx_dict = {rna: i for i, rna in enumerate(
		# 	mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		# gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
		# 						 cistron_ids]
		#
		# # Load data
		# time = read_stacked_columns(cell_paths, 'Main', 'time', ignore_exception=True)
		# (ip_monomer_counts,) = read_stacked_bulk_molecules(
		# 	cell_paths, cistron_monomer_ids, ignore_exception=True)
		# ip_mRNA_counts = read_stacked_columns(
		# 	cell_paths, 'RNACounts', 'mRNA_cistron_counts', ignore_exception=True)[:, gene_mRNA_indexes]





		# # Plotting
		# plt.figure(figsize = (8.5, 11))
		# # Get doubling times from cells with this variant index
		# dt = read_stacked_columns(
		# 	cell_paths, 'Main', 'time',
		# 	fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
		# dts = np.zeros(len(dt))
		# for i in range(len(dt)):
		# 	if i == 0:
		# 		gt = dt[i]
		# 		dts[i] = gt
		# 	else:
		# 		gt = dt[i] + dts[i-1]
		# 		dts[i] = gt
		#
		# # Protein Counts
		# plt.subplot(2, 1, 1)
		# for x in dts:
		# 	plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		# if len(cistron_monomer_ids) == 1:
		# 	plt.plot(time / 60., ip_monomer_counts,
		# 			 label = cistron_monomer_ids[0])
		# else:
		# 	for m in range(len(cistron_monomer_ids)):
		# 		plt.plot(time / 60., ip_monomer_counts[:,m],
		# 				 label = cistron_monomer_ids[m])
		#
		# plt.xlabel("Time (min)")
		# plt.ylabel("Protein Counts")
		# plt.title(f"Protein Counts for Proteins of Interest in variant"
		# 		  f" {variant}, seed {seed}")
		# plt.legend()
		#
		# # mRNA Counts
		# plt.subplot(2, 1, 2)
		# for x in dts:
		# 	plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		# if len(cistron_ids) == 1:
		# 	plt.plot(time / 60., ip_mRNA_counts,
		# 			 label=cistron_ids[0])
		# else:
		# 	for r in range(len(cistron_ids)):
		# 		plt.plot(time / 60., ip_mRNA_counts[:,r],
		# 				 label = cistron_ids[r])
		# plt.xlabel("Time (min)")
		# plt.ylabel("Cistron Counts")
		# plt.title(f"mRNA Counts for Proteins of Interest in variant {cell_paths},"
		# 		  f" seed {seed}")
		# plt.legend()
		#
		# plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		# exportFigure(plt, plotOutDir, plotOutFileName + '_cohort_seed_testing_' + str(seed), metadata)
		# plt.close("all")

if __name__ == '__main__':
	Plot().cli()
