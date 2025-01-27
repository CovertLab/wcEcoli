"""
Plot mRNA and complexed monomer counts across seeds and generations. The
predicted # of monomers in complexes at a given point in time is calculated by
subtracting the free monomer counts from the total monomer counts.

NOTE: it is tough to do this plot with out splitting up the seeds because the
doubling times are different across different seeds (so they will not line up
nicely).

See free_monomer_counts.py for a similar plot that plots free monomer counts,
and see total_monomer_counts.py for a similar plot that plots the total
counts for monomers.
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

# Replace with the proteins you would like to visualize here:
interest_proteins = np.array([
	#'PD03938[c]', # metR
	#'RPOS-MONOMER[c]', # rpoS
	#"BASR-MONOMER[c]", # basR
	#"EG11171-MONOMER[c]", #tsaD
	#"EG11734-MONOMER[c]", # phoH
	#"EG10871-MONOMER[c]", #rplJ
	#"EG11534-MONOMER[c]", # ibpA
	"G6463-MONOMER[c]", # ClpS
])

# Specifiy generations to be skipped if desired:
SKIP_GENERATIONS = 2

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def extract_data(self, simOutDir, cell_paths, monomer_ids, cistron_ids):
		"""
		Extracts the time, complexed monomer counts, and mRNA counts data for
		a given set of cells (typically within a seed) for the protein(s)
		of interest.

		Args:
			simOutDir: output directory for the simulation
			cell_paths: paths to the cells of interest
			monomer_ids: ids for the proteins of interest
			cistron_ids: ids for the genes of interest

		Returns: The time, compelexed monomer counts, and mRNA counts for the
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
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}
		mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
						cistron_ids]

		# Load the time data
		time = read_stacked_columns(cell_paths, 'Main',
									'time', ignore_exception=True)
		# Get the total counts for each protein
		total_monomer_counts = (
								   read_stacked_columns(cell_paths,
														'MonomerCounts',
														'monomerCounts',
														ignore_exception=True)
							   )[:, monomer_indexes]
		# Get the free monomer counts for each protein
		(free_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, monomer_ids, ignore_exception=True)
		# Get the mRNA counts for each gene/protein
		mRNA_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_cistron_counts',
			ignore_exception=True)[:, mRNA_indexes]

		# reshape free_monomer_counts if there is only one protein:
		if len(monomer_ids) == 1:
			free_monomer_counts = free_monomer_counts.reshape(-1, 1)

		# calculate the predicted counts for the proteins in complexes
		complexed_monomer_counts = (
				total_monomer_counts - free_monomer_counts)

		return (time, total_monomer_counts, free_monomer_counts,
				complexed_monomer_counts, mRNA_counts)

	def extract_doubling_times(self, seed):
		"""
		Extracts doubling time data for a given set of cells.
		Args:
			seed: seed of interest.

		Returns: the doubling time for each generation within the seed and
		the end time for each generation.

		"""
		# redefine cell_paths to be the first seed
		cell_paths = self.ap.get_cells(seed=[seed])

		# Get doubling times for the cells with this seed index
		dt = read_stacked_columns(
			cell_paths, 'Main', 'time',
			fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()

		# Determine the end time for each generation
		generation_ends = np.zeros(len(dt))
		for i in range(len(dt)):
			if i == 0:
				gt = dt[i]
				generation_ends[i] = gt
			else:
				gt = dt[i] + generation_ends[i - 1]
				generation_ends[i] = gt

		dt = dt[SKIP_GENERATIONS:]
		generation_ends = generation_ends[SKIP_GENERATIONS:]

		return generation_ends, dt

	def plot_counts_per_seed(self, seed, cell_paths, monomer_ids, cistron_ids,
							 time, complexed_monomer_counts, mRNA_counts,
							 plotOutDir, plotOutFileName, metadata):
		"""
		Plots the complexed monomer counts and mRNA counts for a given seed.

		Args:
			seed: the seed index being plotted
			cell_paths: paths to the cells within the seed
			monomer_ids: monomer ids for the proteins of interest
			cistron_ids: cistron ids for the proteins of interest
			time: the total time duration for the seed
			complexed_monomer_counts: complexed monomer counts for proteins
			of interest over the seed duration
			mRNA_counts: mRNA counts for the proteins of interest over the
			seed duration
			plotOutDir: directory for the plot output
			plotOutFileName: name for the plot output
			metadata: metadata for the simulation

		Returns: A plot export for the complexed monomer counts and mRNA
		counts in the particular seed.

		"""
		plt.figure(figsize=(8.5, 11))

		# Get doubling times for the cells with this seed index
		dts, dt_duration = self.extract_doubling_times(seed)

		# Protein Counts Plot
		plt.subplot(2, 1, 1)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)

		gene_list = []
		if len(monomer_ids) == 1:
			common_name = self.get_common_name(monomer_ids[0])
			name = monomer_ids[0] + ' (' + common_name + ')'
			gene_list.append(common_name)
			plt.plot(time / 60., complexed_monomer_counts,
					 label=name)
		else:
			for m in range(len(monomer_ids)):
				common_name = self.get_common_name(monomer_ids[m])
				name = monomer_ids[m] + ' (' + common_name + ')'
				gene_list.append(common_name)
				plt.plot(time / 60., complexed_monomer_counts[:, m],
						 label=name, alpha=0.5)

		plt.xlabel("Time (min)");
		plt.ylabel("Complexed Monomer Counts")
		plt.title(f"Predicted Counts of Monomers within"
				  f" Complexes in Seed {seed}, starting with generation"
				  f" {SKIP_GENERATIONS}")
		plt.legend()

		# mRNA Counts Plot
		plt.subplot(2, 1, 2)
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=2)
		if len(cistron_ids) == 1:
			common_name = self.get_common_name(monomer_ids[0])
			name = cistron_ids[0] + ' (' + common_name + ')'
			plt.plot(time / 60., mRNA_counts,
					 label=name)
		else:
			for r in range(len(cistron_ids)):
				common_name = self.get_common_name(monomer_ids[r])
				name = cistron_ids[r] + ' (' + common_name + ')'
				plt.plot(time / 60., mRNA_counts[:, r],
						 label=name, alpha=0.5)
		plt.xlabel("Time (min)")
		plt.ylabel("Cistron Counts")
		plt.title(f"mRNA Counts for Proteins of Interest in Seed {seed}, "
				  f"starting with generation {SKIP_GENERATIONS}")
		plt.legend()

		plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
		plotOutPath = os.path.join(plotOutDir, f'00000{seed}')
		exportFigure(plt, plotOutPath, plotOutFileName + f'_cohortPlot_'
					f'startGen_{SKIP_GENERATIONS}_seed_'
					 + str(seed) + '_geneIDs_' + str(gene_list), metadata)
		plt.close("all")

	def plot_all_count_types_per_seed(self, seed, cell_paths, monomer_ids, time,
									  total_monomer_counts, free_monomer_counts,
									  complexed_monomer_counts, plotOutDir,
									  plotOutFileName, metadata):
		"""
		Plots the total, free, and complexed monomer counts for a given seed.
		Args:
			seed: current seed to be plotted within the cohort
			cell_paths: paths to cells within the selected seed
			monomer_ids: monomer IDs for the proteins of interest
			time: time steps spanning the cells within the seed
			total_monomer_counts: total monomer counts at each timestep
			free_monomer_counts: free monomer counts at each timestep
			complexed_monomer_counts: predicted counts of monomers within
			complexes at each timestep (calculated by subtracting the free
			monomer counts from the total).
			plotOutDir: directory for the plot output
			plotOutFileName: name for the plot output
			metadata: metadata for the simulation

		Returns: A plot export for the total, free, and predicted
		complexed monomer counts for each protein in the seed

		"""
		plt.figure(figsize=(8.5, 14))

		# Get doubling times for the cells with this seed index
		dts, dt_duration = self.extract_doubling_times(seed)

		# Total Monomer Counts Plot
		plt.subplot(3, 1, 1)

		# plot the doubling times
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=1,
						alpha=0.5)

		# plot the data
		gene_list = []
		if len(monomer_ids) == 1:
			common_name = self.get_common_name(monomer_ids[0])
			name = monomer_ids[0] + ' (' + common_name + ')'
			gene_list.append(common_name)
			plt.plot(time / 60., total_monomer_counts,
					 label=name)
		else:
			for m in range(len(monomer_ids)):
				common_name = self.get_common_name(monomer_ids[m])
				name = monomer_ids[m] + ' (' + common_name + ')'
				gene_list.append(common_name)
				plt.plot(time / 60., total_monomer_counts[:, m],
						 label=name, alpha=0.9)

		plt.xlabel("Time (min)");
		plt.ylabel("Total Monomer Counts")
		plt.title(f"Total Monomer Counts in Seed {seed}, starting with "
				  f"generation {SKIP_GENERATIONS}")
		plt.legend()

		# Free Monomer Counts Plot
		plt.subplot(3, 1, 2)

		# plot the doubling times
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=1,
						alpha=0.5)

		# plot the data
		if len(monomer_ids) == 1:
			common_name = self.get_common_name(monomer_ids[0])
			name = monomer_ids[0] + ' (' + common_name + ')'
			plt.plot(time / 60., free_monomer_counts,
					 label=name, alpha=0.5)
		else:
			for m in range(len(monomer_ids)):
				common_name = self.get_common_name(monomer_ids[m])
				name = monomer_ids[m] + ' (' + common_name + ')'
				plt.plot(time / 60., free_monomer_counts[:, m],
						 label=name, alpha=0.5)

		plt.xlabel("Time (min)")
		plt.ylabel("Free Monomer Counts")
		plt.title(f"Free Monomer Counts in Seed {seed}, starting with generation"
				  f" {SKIP_GENERATIONS}")
		plt.legend()

		# Complexed Monomer Counts Plot
		plt.subplot(3, 1, 3)

		# plot the doubling times
		for x in dts:
			plt.axvline(x=x, color='#bcbd22', linestyle='--', linewidth=1,
						alpha=0.5)

		# plot the data
		if len(monomer_ids) == 1:
			common_name = self.get_common_name(monomer_ids[0])
			name = monomer_ids[0] + ' (' + common_name + ')'
			plt.plot(time / 60., complexed_monomer_counts,
					 label=name)
		else:
			for m in range(len(monomer_ids)):
				common_name = self.get_common_name(monomer_ids[m])
				name = monomer_ids[m] + ' (' + common_name + ')'
				plt.plot(time / 60., complexed_monomer_counts[:, m],
						 label=name, alpha=0.9)

		plt.xlabel("Time (min)")
		plt.ylabel("Complexed Monomer Counts")
		plt.title(f"Predicted Counts for Monomers within Complexes"
				  f" in Seed {seed}, starting with generation {SKIP_GENERATIONS}")
		plt.legend()

		plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
		plotOutPath = os.path.join(plotOutDir, f'00000{seed}')
		exportFigure(plt, plotOutPath, plotOutFileName +
					 f'_cohortPlot_compareCounts_startGen_{SKIP_GENERATIONS}'
					 f'_seed_'+str(seed)+'_geneIDs_'+str(gene_list), metadata)
		plt.close("all")

	# function to match gene symbols to monomer ids
	def get_gene_symbols_for_monomer_ids(self):
		"""
		Extracts the gene symbols for each monomer id in the model.

		Returns: a dictionary mapping monomer ids to gene symbols.
		Code adapted from convert_to_flat.py.
		"""
		RNAS_FILE = os.path.join(ROOT_PATH, 'reconstruction', 'ecoli',
								 'flat', 'rnas.tsv')
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

	def get_common_name(self, protein_id):
		"""
        Obtains the common names for each protein of interest
        Args:
            protein_id: the name of the protein(s) of interest

        Returns: the common name for the protein(s) of interest
        """
		protein = protein_id[:-3]
		gene_symbol = self.get_gene_symbols_for_monomer_ids()[protein]
		return gene_symbol

	def plot_counts_per_protein(self, simOutDir, protein, cistron, seeds,
								generations, plotOutDir, plotOutFileName,
								metadata):
		"""
        Plot the complexed monomer counts and mRNA counts for a given protein
        across all seeds in one plot. Also includes a rough average of the
        counts per generation on the plot and a table with these values
        averaged across all seeds for each generation.
        Args:
            simOutDir: directory for the simulation data
            protein: protein ID for the protein of interest
            cistron: cistron ID for the protein of interest
            seeds: seeds in the simulation
            generations: generations in the simulation to plot
			plotOutDir: directory for the plot output
			plotOutFileName: name for the plot output
			metadata: metadata for the simulation

        Returns: A plot of the complexed protein counts and mRNA counts for the
        protein of interest across all seeds in the simulation.
        """
		# get the common name for the protein
		common_name = self.get_common_name(protein[0])

		# build these for the average value table:
		time_values = []
		dts_values = []
		dt_duration_values = []
		monomer_values = []
		mRNA_values = []

		plt.figure(figsize=(8.5, 11))

		# Protein Counts Plot
		plt.subplot(2, 1, 1)

		# loop through the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], generation=generations,
								  only_successful=True))
			# Extract doubling time for the seed
			dts, dt_durations = self.extract_doubling_times(seed)
			# Extract data for the seed
			(time, total_monomer_counts, free_monomer_counts,
			 complexed_monomer_counts, mRNA_counts) = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))

			# save the data for the average value table:
			dts_values.append(dts)
			dt_duration_values.append(dt_durations)
			time_values.append(time)
			monomer_values.append(complexed_monomer_counts)
			mRNA_values.append(mRNA_counts)

			# plot the data
			plt.plot(time / 60., complexed_monomer_counts,
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

		# Plot specs for the monomer counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Complexed Monomer Counts")
		plt.title(f"Predicted Counts for {protein[0]} ({common_name}) "
				  f"Monomers within Complexes")
		plt.legend()

		# mRNA Counts Plot
		plt.subplot(2, 1, 2)
		# loop through the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], generation=generations, only_successful=True))
			# Extract data for the seed
			(time, total_monomer_counts, free_monomer_counts,
			 complexed_monomer_counts, mRNA_counts) = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))
			# plot the data
			plt.plot(time / 60., mRNA_counts,
					 label=f'seed {seed}', alpha=0.5)

		# Plot specs for the mRNA counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Cistron Counts")
		plt.title(f"mRNA Counts for {cistron[0]} "
				  f"(common name: {common_name})")
		plt.legend()

		# generate a table below the plots:
		columns = ('avg. cycle duration', 'avg. monomer counts',
				   'avg. mRNA counts')
		rows = (['Generation %d' % p for p in
				tuple(range(SKIP_GENERATIONS, (generations[-1] + 1 )))] +
				['Overall Average'])

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
					 f'_cohortPlot_startGen_{SKIP_GENERATIONS}'
					 f'_geneID_' + common_name, metadata)
		plt.close("all")

	def plot_all_count_types_per_protein(self, simOutDir, protein, cistron,
										 seeds, generations, plotOutDir,
										 plotOutFileName, metadata):
		"""
		Plots the total, free, and complexed monomer counts for a given
		protein. Allows one to compare the make up of the total monomer
		counts between free and complexed forms.
		Args:
			simOutDir: directory for the simulation data
			protein: protein of interest
			cistron: corresponding cistron for the protein of interest
			seeds: seeds in the simulation
			generations: generations to plot per seed
			plotOutDir: directory for the plot output
			plotOutFileName: name for the plot output
			metadata: metadata for the simulation

		Returns: A plot with the total, free, and predicted complexed
		monomer counts for a protein for all seeds in one graph.

		"""
		# get the common name for the protein
		common_name = self.get_common_name(protein[0])

		# Plotting
		plt.figure(figsize=(8.5, 14))
		colors = ['darkorange', 'lightseagreen', 'blueviolet',
				  'olive', 'hotpink', ]

		# Total Counts Plot
		plt.subplot(3, 1, 1)

		# loop through the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], generation=generations,
								  only_successful=True))
			# Extract data for the seed
			(time, total_monomer_counts, free_monomer_counts,
			 complexed_monomer_counts, mRNA_counts) = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))
			c = colors[seed]

			# plot the data
			plt.plot(time / 60., total_monomer_counts,
					 label=f'seed {seed}', alpha=0.9,
					 linestyle=':', color=c)

		# Plot specs for the total monomer counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Total Monomer Counts")
		plt.title(f"Total Monomer Counts for {protein[0]} ({common_name})")
		plt.legend()

		# Free Monomer Counts Plot
		plt.subplot(3, 1, 2)

		# loop through the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], generation=generations,
								  only_successful=True))
			# Extract data for the seed
			(time, total_monomer_counts, free_monomer_counts,
			 complexed_monomer_counts, mRNA_counts) = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))
			c = colors[seed]

			# plot the data
			plt.plot(time / 60., free_monomer_counts,
					 label=f'seed {seed}', alpha=0.8,
					 linestyle='--', color=c)

		# Plot specs for the free monomer counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Free Monomer Counts")
		plt.title(f"Free Monomer Counts for {protein[0]} ({common_name})")
		plt.legend()

		# Complexed Monomer Counts Plot
		plt.subplot(3, 1, 3)

		# loop through the seeds:
		for seed in seeds:
			# cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], generation=generations,
								  only_successful=True))
			# Extract data for the seed
			(time, total_monomer_counts, free_monomer_counts,
			 complexed_monomer_counts, mRNA_counts) = (
				self.extract_data(simOutDir, cell_paths, protein, cistron))
			c = colors[seed]

			# plot the data
			plt.plot(time / 60., total_monomer_counts,
					 alpha=1, linestyle=':', color=c, linewidth=0.5)
			plt.plot(time / 60., complexed_monomer_counts,
					 label=f'seed {seed}', alpha=0.9, color=c,
					 linewidth=1)
			plt.plot(time / 60., free_monomer_counts,
					 alpha=0.5, linestyle='--', color=c)

		# Plot specs for the complexed monomer counts graph
		plt.xlabel("Time (min)");
		plt.ylabel("Complexed Monomer Counts")
		plt.title(f"Predicted Counts for {protein[0]} ({common_name}) "
				  f"Monomers within Complexes")
		plt.legend()

		# export plot:
		plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
		exportFigure(plt, plotOutDir, plotOutFileName +
					 f'_cohortPlot_compareCounts_startGen_{SKIP_GENERATIONS}'
					 f'_geneID_' + common_name,
					 metadata)
		plt.close("all")

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):

		# extract data paths
		seeds = self.ap.get_seeds()
		generations = self.ap.get_generations()[SKIP_GENERATIONS:]
		cell_paths = self.ap.get_cells(generation=generations,
									   only_successful=True)

		# get the simOut directory
		sim_dir = cell_paths[0] # this can be arbitrary, just needs to exist
		simOutDir = os.path.join(sim_dir, 'simOut')

		# load data needed to determine gene/cistron ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		monomer_sim_data = (
			sim_data.process.translation.monomer_data.struct_array)

		# extract info about the protein(s) from the monomer data:
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

		# iterate through the seeds to generate graphs for each seed
		for seed in seeds:
			# Generate all the cell paths for the seed:
			cell_paths = (
				self.ap.get_cells(seed=[seed], generation=generations,
								  only_successful=True))
			# Extract data for the seed
			(time, total_monomer_counts, free_monomer_counts,
			 complexed_monomer_counts, mRNA_counts) = (
				self.extract_data(simOutDir, cell_paths, monomer_ids,
							 cistron_ids))
			# Plot the counts for the seed
			self.plot_counts_per_seed(seed, cell_paths, monomer_ids, cistron_ids,
								 time, complexed_monomer_counts , mRNA_counts,
									  plotOutDir, plotOutFileName, metadata)

			# Plot comparisons each type of count per seed
			self.plot_all_count_types_per_seed(seed, cell_paths, monomer_ids,
												time, total_monomer_counts,
												free_monomer_counts,
												complexed_monomer_counts,
											   plotOutDir, plotOutFileName,
											   metadata)

		# Now, plot for each individual protein across all seeds
		for protein in interest_proteins:
			# get the the cistron id for the protein
			cistron_id = \
				[key for key,
				v in cistron_monomer_id_dict.items() if v == protein]
			protein_id = [protein]  # fomrat the protein id for extract_data

			# plot the counts for the protein
			print("plotting for protein: "); print(protein)
			print("cistron_id: "); print(cistron_id[0])
			self.plot_counts_per_protein(simOutDir, protein_id, cistron_id, seeds,
										 generations,
										 plotOutDir, plotOutFileName, metadata)
			self.plot_all_count_types_per_protein(simOutDir, protein_id,
												  cistron_id, seeds,
												  generations, plotOutDir,
												  plotOutFileName, metadata)

if __name__ == '__main__':
	Plot().cli()
