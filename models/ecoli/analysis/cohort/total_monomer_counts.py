"""
Plot mRNA and total protein counts for monomers across seeds and generations

NOTE: it is tough to do this plot with out splitting up the seeds because the
doubling times are different across different seeds (so they will not line up
nicely)

See free_monomer_counts.py for a similar plot that plots free monomer counts,
and see complexed_monomer_counts.py for a similar plot that plots the predicted
counts for monomers in complexes as well.
"""

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

# List the proteins you would like to have plotted here:
interest_proteins = np.array([
	#"EG10159-MONOMER[c]", #"ATP-dependent Clp protease ATP-binding subunit ClpX"
	'EG10542-MONOMER[c]', # lon
	#'PD03938[c]', # metR
	#'RPOS-MONOMER[c]', # rpoS
	#'EG11171-MONOMER[c]', # tsaD, has interesting behavior
	#"BASR-MONOMER[c]", # basR, also interesting behavior

])

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):

	def extract_data(self, simOutDir, cell_paths, monomer_ids, cistron_ids):
		"""
        Extracts the time, total monomer counts, and mRNA counts data for
        a given set of cells (typically within a seed) for the protein(s)
        of interest.

        Args:
            simOutDir: output directory for the simulation
            cell_paths: paths to the cells of interest
            monomer_ids: ids for the proteins of interest
            cistron_ids: ids for the genes of interest

        Returns: The time, total monomer counts, and mRNA counts for the
        proteins of interest.
        """
		# Extract monomer indexes for each protein of interest
		monomer_counts_reader = TableReader(os.path.join(simOutDir,
														 'MonomerCounts'))
		monomer_idx_dict = {monomer: i for i, monomer in enumerate(
			monomer_counts_reader.readAttribute('monomerIds'))}
		monomer_indexes = [monomer_idx_dict.get(monomer_id) for
						   monomer_id in monomer_ids]

		# Extract mRNA indexes for each gene of interest
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'RNACounts'))
		mRNA_idx_dict = {rna: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in cistron_ids]

		# Load the time data across the duration of all simulated cells
		time = read_stacked_columns(cell_paths, 'Main','time',
									ignore_exception=True)
		# Get the total counts for each protein
		monomer_counts = (read_stacked_columns(
			cell_paths,'MonomerCounts','monomerCounts',
			ignore_exception=True))[:, monomer_indexes]
		# Get the mRNA counts for each gene
		mRNA_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts',
			ignore_exception=True)[:, mRNA_indexes]

		return time, monomer_counts, mRNA_counts

	def extract_doubling_times(self, cell_paths):
		"""
        Extracts doubling time data for a given set of cells.
        Args:
            cell_paths: paths to cells within a seed.

        Returns: the doubling time for each generation within the seed and
        the end time for each generation (within the total time spanned by all
        cell generations).
        """
		# determine the doubling times for each cell
		dts = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()

		# determine the end time of each cell generation
		generation_ends = np.zeros(len(dts))
		for i in range(len(dts)):
			if i == 0:
				gt = dts[i]
				generation_ends[i] = gt
			else:
				gt = dts[i] + generation_ends[i - 1]
				generation_ends[i] = gt

		return generation_ends, dts

	def plot_counts_per_seed(self, seed, cell_paths, monomer_ids, cistron_ids,
							 time, monomer_counts, mRNA_counts, plotOutDir,
							 plotOutFileName, metadata):
		"""
        Plots the total monomer counts and mRNA counts for a given seed.

        Args:
            seed: the seed index being plotted
            cell_paths: paths to the cells within the seed
            monomer_ids: monomer ids for the proteins of interest
            cistron_ids: cistron ids for the proteins of interest
            time: the total time duration for the seed
            monomer_counts: total monomer counts for proteins of interest
            over the seed duration
            mRNA_counts: mRNA counts for the proteins of interest over the seed
            plotOutDir: output directory for the plot
            plotOutFileName: name of the plot file
            metadata: metadata for the simulation

        Returns: A plot export for the total monomer counts and mRNA counts
        in the particular seed.
        """
		plt.figure(figsize=(8.5, 11))

		# Get doubling times for the cells with this seed index
		dts, dt_duration = self.extract_doubling_times(cell_paths)

		# Protein Counts Plot
		plt.subplot(2, 1, 1)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)

		gene_list = []
		if len(monomer_ids) == 1:
			common_name = self.get_common_names(monomer_ids[0])
			name = monomer_ids[0] + ' (' + common_name + ')'
			gene_list.append(common_name)
			plt.plot(time / 60., monomer_counts,
					 label=name)
		else:
			for m in range(len(monomer_ids)):
				common_name = self.get_common_names(monomer_ids[m])
				name = monomer_ids[m] + ' (' + common_name + ')'
				gene_list.append(common_name)
				plt.plot(time / 60., monomer_counts[:, m],
						 label=name, alpha=0.5)

		plt.xlabel("Time (min)");
		plt.ylabel("Total Monomer Counts")
		plt.title(f"Total Monomer Counts in Seed {seed}")
		plt.legend()

		# mRNA Counts Plot
		plt.subplot(2, 1, 2)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		if len(cistron_ids) == 1:
			common_name = self.get_common_names(monomer_ids[0])
			name = cistron_ids[0] + ' (' + common_name + ')'
			plt.plot(time / 60., mRNA_counts,
					 label=name)
		else:
			for r in range(len(cistron_ids)):
				common_name = self.get_common_names(monomer_ids[r])
				name = cistron_ids[r] + ' (' + common_name + ')'
				plt.plot(time / 60., mRNA_counts[:, r],
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

	def get_gene_symbols_for_monomer_ids(self):
		"""
        Extracts the gene symbols for each monomer id in the model.
        Returns: a dictionary mapping monomer ids to gene symbols.
        Code adapted from convert_to_flat.py.
        """
		RNAS_FILE = os.path.join(
			ROOT_PATH, 'reconstruction', 'ecoli', 'flat', 'rnas.tsv')
		with (io.open(RNAS_FILE, 'rb') as f):
			reader = tsv.reader(f, delimiter='\t')
			headers = next(reader)
			while headers[0].startswith('#'):
				headers = next(reader)

			# extract relevant information
			gene_symbol_index = headers.index('common_name')
			protein_id_index = headers.index('monomer_ids')
			monomer_ids_to_gene_symbols = {}
			for line in reader:
				gene_symbol = line[gene_symbol_index]
				protein_id = list(
					line[protein_id_index][2:-2].split('", "'))[0]
				monomer_ids_to_gene_symbols[protein_id] = gene_symbol

			return monomer_ids_to_gene_symbols

	def get_common_names(self, protein_id):
		"""
        Obtains the common names for each protein of interest
        Args:
            protein_id: the name of the protein(s) of interest
        Returns: the common name for the protein(s) of interest
        """
		protein = protein_id[:-3]  # subtract the compartment
		common_name = self.get_gene_symbols_for_monomer_ids()[protein]

		return common_name

	def plot_counts_per_protein(self, simOutDir, protein, cistron, seeds,
								plotOutDir, plotOutFileName, metadata):
		"""
        Plot the total monomer counts and mRNA counts for a given protein
        across all seeds in one plot.
        Args:
            simOutDir: directory for the simulation data
            protein: protein ID for the protein of interest
            cistron: cistron ID for the protein of interest
            seeds: seed indices for the simulation
            plotOutDir: output directory for the plot
            plotOutFileName: name of the plot file
            metadata: metadata for the simulation


        Returns: A plot of the total protein counts and mRNA counts for the
        protein of interest across all seeds in the simulation.

        """
		# get the common name for the protein
		common_name = self.get_common_names(protein[0])

		# build these for the average value table:
		time_values = [];
		dts_values = [];
		dt_duration_values = []
		monomer_values = [];
		mRNA_values = []

		plt.figure(figsize=(8.5, 11))

		# Protein Counts Plot
		plt.subplot(2, 1, 1)

		# loop through all the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], only_successful=True))
			# Extract doubling time for the seed
			dts, dt_durations = self.extract_doubling_times(cell_paths)
			# Extract data for the seed
			time, monomer_counts, mRNA_counts = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))

			# save the data for the average value table:
			dts_values.append(dts)
			dt_duration_values.append(dt_durations)
			time_values.append(time)
			monomer_values.append(monomer_counts)
			mRNA_values.append(mRNA_counts)

			# plot the data
			plt.plot(time / 60., monomer_counts,
					 label=f'seed {seed}', alpha=0.5)

		# find the average quanities across all seeds
		average_generation_durations = np.mean(dt_duration_values, axis=0)

		# get the average counts for each generation across all seeds:
		average_monomer_cts = []
		average_mRNA_cts = []
		for seed in range(len(seeds)):  # range(len(seeds)) gets seed idx
			seed_dts = dts_values[seed]
			seed_time = time_values[seed]
			seed_monomer = monomer_values[seed]
			seed_mRNA = mRNA_values[seed]

			# get the average monomer and mRNA counts for each generation:
			all_gen_monomer_cts = []
			all_gen_mRNA_cts = []
			for gen in range(len(average_generation_durations)):
				if gen == 0:
					time_range = ((seed_time >= 0) &
								  (seed_time <= (seed_dts[0] * 60.)))
				else:
					time_range = ((seed_time >= ((seed_dts[gen - 1]) * 60.))
								  & (seed_time <= ((seed_dts[gen]) * 60.)))

				# get the avg. monomer and mRNA counts for each generation:
				time_idxs = np.where(time_range)
				all_gen_monomer_cts.append(
					np.mean(seed_monomer[time_idxs[0]]))
				all_gen_mRNA_cts.append(np.mean(seed_mRNA[time_idxs[0]]))

			average_monomer_cts.append(all_gen_monomer_cts)
			average_mRNA_cts.append(all_gen_mRNA_cts)

		# average the monomer and mRNA counts across all seeds:
		avg_monomer_cts = np.mean(average_monomer_cts, axis=0)
		avg_mRNA_cts = np.mean(average_mRNA_cts, axis=0)

		# plot at the middle of each average generation duration:
		avg_time_points = []
		for i in range(len(average_generation_durations)):
			if i == 0:
				avg_time_points.append(
					average_generation_durations[i] / 2.)
			else:
				summed_time = np.sum(average_generation_durations[:i])
				avg_time_points.append(
					summed_time + ((average_generation_durations[i]) / 2.))

		# plot the average counts for the protein across all seeds
		plt.scatter(avg_time_points, avg_monomer_cts,
					label='average', color='black', marker='.')

		# Plot specs for the monomer counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Total Monomer Counts")
		plt.title(f"Total Monomer Counts for {protein[0]} ({common_name})")
		plt.legend()

		# mRNA Counts Plot
		plt.subplot(2, 1, 2)
		# loop through the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], only_successful=True))
			# Extract data for the seed
			time, monomer_counts, mRNA_counts = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))
			# plot the data
			plt.plot(time / 60., mRNA_counts,
					 label=f'seed {seed}', alpha=0.5)

		# plot the average counts for the protein across all seeds
		plt.scatter(avg_time_points, avg_mRNA_cts, label='average',
					color='black', marker='.')

		# Plot specs for the mRNA counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Cistron Counts")
		plt.title(f"mRNA Counts for {cistron[0]} "
				  f"(common name: {common_name})")
		plt.legend()

		# generate a table below the plots:
		columns = ('avg. cycle duration', 'avg. monomer counts',
				   'avg. mRNA counts')
		rows = ['Generation %d' % p for p in
				(1, 2, 3, 4, 5, 6, 7, 8)] + ['avg. of Averages']

		# add the overall average to each:
		avg_generation_durations_all = np.mean(average_generation_durations)
		avg_monomer_cts_all = np.mean(avg_monomer_cts)
		avg_mRNA_cts_all = np.mean(avg_mRNA_cts)

		# add the overall average to the end of the arrays:
		avg_generation_durations = np.append(average_generation_durations,
											 avg_generation_durations_all)
		avg_monomer_cts = np.append(avg_monomer_cts, avg_monomer_cts_all)
		avg_mRNA_cts = np.append(avg_mRNA_cts, avg_mRNA_cts_all)

		# create the data for the table:
		table_data = np.array([avg_generation_durations,
							   avg_monomer_cts, avg_mRNA_cts]).T

		# create the table:
		plt.table(cellText=table_data, rowLabels=rows, colLabels=columns,
				  loc='bottom', bbox=[0.1, -0.75, 1, 0.5])

		# export plot:
		plt.subplots_adjust(hspace=0.2, top=0.95, bottom=0.25)
		exportFigure(plt, plotOutDir, plotOutFileName +
					 '_cohortPlot_geneID_' + common_name, metadata)
		plt.close("all")

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# extract paths to cells and the seed IDs
		cell_paths = self.ap.get_cells(only_successful=True)
		seeds = self.ap.get_seeds()

		# obtain the simOut directory
		sim_dir = cell_paths[0] # this can be arbitrary, just needs to exist
		simOutDir = os.path.join(sim_dir, 'simOut')

		# load data needed to determine gene/cistron ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)

		# extract info about the protein(s) from the monomer/cistron data:
		monomer_data_idxs = []
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			monomer_data_idxs.append(monomer_idx)
		# monomer_data includes helpful info: monomer_data.dtype.names
		monomer_data = monomer_sim_data[monomer_data_idxs]
		monomer_ids = monomer_data['id']
		cistron_ids = monomer_data['cistron_id']

		# create a dictionary to map cistron ids to monomer ids
		cistron_monomer_id_dict = dict(zip(monomer_sim_data['cistron_id'],
										   monomer_sim_data['id']))

		# iterate through the seeds to total monomer count graphs for each seed:
		for seed in seeds:
			# Generate all the cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], only_successful=True))
			# Extract data for the seed
			time, monomer_counts, mRNA_counts = (self.extract_data(
				simOutDir, cell_paths, monomer_ids, cistron_ids))
			# Plot the counts for the specific seed
			self.plot_counts_per_seed(
				seed, cell_paths, monomer_ids, cistron_ids, time, monomer_counts,
				mRNA_counts, plotOutDir, plotOutFileName, metadata)

		# iterate through the proteins and plot counts across all seeds:
		for protein in interest_proteins:
			# get the the cistron id for the protein
			cistron_id = [
				key for key, v in cistron_monomer_id_dict.items() if v == protein]
			protein_id = [protein]  # format the protein id for extract_data()

			# plot the counts for the protein
			print("plotting for protein: "); print(protein)
			print("cistron_id: "); print(cistron_id[0])
			self.plot_counts_per_protein(
				simOutDir, protein_id, cistron_id, seeds,
				plotOutDir, plotOutFileName, metadata)


if __name__ == '__main__':
	Plot().cli()
