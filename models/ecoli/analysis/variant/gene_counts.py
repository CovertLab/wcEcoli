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
	# 'G7634-MONOMER[i]',
	#'EG11854-MONOMER[c]',
	#'G6606-MONOMER[c]',
	#'MONOMER0-2678[c]',
	#'EG10037-MONOMER[c]',
	'NG-GFP-MONOMER[c]',
	#'TRYPSYN-APROTEIN[c]',
	#"ANTHRANSYNCOMPI-MONOMER[c]"
	#'PD00519[c]',
])

show_seed_plot = 1
show_gen_plot = 1

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def seed_plot(self, control_var, experimental_vars):
		# retrieve the seed data
		seeds = self.ap.get_seeds()
		# Plotting
		plt.figure(figsize=(8.5, 11))
		colors = ["turquoise", "yellowgreen", "mediumpurple", "deeppink",
				  "deepskyblue","lightcoral","gold", "darkorange", "darkred",
				  "darkgreen", "darkblue", "darkviolet", "darkturquoise",]
		# Protein Counts
		plt.subplot(2, 1, 1)
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
				plt.plot(ctime / 60., cip_monomer_counts,
						 linestyle=LS[seed],
						 label=name, color="#FF796C", linewidth=.5)
				for variant in experimental_vars:
					seed_dir = self.ap.get_cells(variant=[variant], seed=[seed],
												 only_successful=True)
					# Load data for the seed
					time = read_stacked_columns(seed_dir, 'Main', 'time',
												ignore_exception=True)
					(ip_monomer_counts,) = read_stacked_bulk_molecules(
						seed_dir, self.cistron_monomer_ids, ignore_exception=True)
					name = 'seed ' + str(seed) + ' exp  var ' + str(variant)
					plt.plot(time / 60., ip_monomer_counts, label=name,
							 color=colors[variant-1], linestyle=LS[seed], linewidth=.5)
				# plot specs
				plt.title(f"Protein Counts for {self.cistron_monomer_ids[0]} in "
						  f"\nthe control variant"
						  f" and {len(experimental_vars)} experimental variant(s)")

			else:
				LS = ['-', ':', '-.', '--']
				for m in range(len(self.cistron_monomer_ids)):
					ls = LS[m]
					name = (self.cistron_monomer_ids[m] + ' seed ' +
							str(seed) + ' cntrl var ')
					plt.plot(ctime / 60., cip_monomer_counts[:, m],
							 linestyle=ls, label=name, color="#FF796C",
							 linewidth=.5)
					for variant in experimental_vars:
						seed_dir = self.ap.get_cells(variant=[variant], seed=[seed],
													 only_successful=True)
						time = read_stacked_columns(seed_dir, 'Main', 'time',
													ignore_exception=True)
						(ip_monomer_counts,) = read_stacked_bulk_molecules(
							seed_dir, self.cistron_monomer_ids, ignore_exception=True)
						name = (self.cistron_monomer_ids[m] + ' seed ' +
								str(seed) + ' exp  var ' + str(variant))
						plt.plot(time / 60., ip_monomer_counts[:, m],
								 label=name, linestyle=ls, color=colors[variant-1],
								 linewidth=.5)
				# plot specs
				plt.title(f"Protein Counts for {self.cistron_monomer_ids} in the"
						  f"\n control variant"
						  f" and {len(experimental_vars)} experimental variant(s)")

		# finish protein count plot
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.legend()

		# mRNA Counts
		plt.subplot(2, 1, 2)
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
				plt.plot(ctime / 60., cip_mRNA_counts, linestyle=LS[seed],
						 label=name, color="#FF796C", linewidth=.5)
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
					plt.plot(time / 60., ip_mRNA_counts, label=name,
							 color=colors[variant-1], linestyle=LS[seed], linewidth=.5)
				# plot specs
				plt.title(f"mRNA Counts for {self.cistron_ids[0]} in the \n "
						  f"control variant"
						  f" and {len(experimental_vars)} experimental variant(s)")
			else:
				LS = ['-', ':', '-.', '--']
				for m in range(len(self.cistron_ids)):
					ls = LS[m]
					name = (self.cistron_ids[m] + ' seed ' + str(seed) +
							' cntrl var ')
					plt.plot(ctime / 60., cip_mRNA_counts[:, m],
							 linestyle=ls, label=name, color="#FF796C",
							 linewidth=.4)
					for variant in experimental_vars:
						seed_dir = self.ap.get_cells(variant=[variant],
													 seed=[seed],
													 only_successful=True)
						time = read_stacked_columns(seed_dir, 'Main',
													'time',
													ignore_exception=True)
						ip_mRNA_counts = read_stacked_columns(
							seed_dir, 'RNACounts', 'mRNA_cistron_counts',
							ignore_exception=True)[:, self.new_gene_mRNA_indexes]
						name = (self.cistron_ids[m] + ' seed ' + str(seed) +
								' exp  var ' + str(variant))
						plt.plot(time / 60., ip_mRNA_counts[:, m],
								 label=name, linestyle=ls, color=colors[variant-1],
								 linewidth=.4)
				# plot specs
				plt.title(f"mRNA Counts for {self.cistron_ids} in the\n"
						  f" control variant"
						  f" and {len(experimental_vars)} experimental variant(s)")

		plt.xlabel("Time (min)")
		plt.ylabel("Cistron Counts")
		plt.legend()
		plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)

	def generate_gen_data(self, experimental_var):
		# get the seeds and generations for the experimental variant:
		generations = self.ap.n_generation # this is the max number of gens
		seeds = self.ap.get_seeds(variant=[experimental_var])
		gen_data = np.zeros((generations, 4), dtype=object)

		for gen in range(0, generations, 1):
			time_data = np.zeros((len(seeds), 1), dtype=object) - 1.
			nanarray = []
			for i in range(len(interest_proteins)):
				nanarray.append(np.nan)
			monomer_data = np.zeros((len(seeds), 1), dtype=object)
			mRNA_data = np.zeros((len(seeds), 1), dtype=object)
			for s in seeds:
				seed_dir = self.ap.get_cells(variant=[experimental_var],
											 seed=[s], generation=[gen],
											 only_successful=True)
				if seed_dir == []:
					time_data[s] = -1
					monomer_data[s] = nanarray
					mRNA_data[s] = nanarray
				else:
					# Load data for the seed
					time = read_stacked_columns(seed_dir, 'Main',
												'time',
												ignore_exception=True)
					(ip_monomer_counts,) = read_stacked_bulk_molecules(
						seed_dir, self.cistron_monomer_ids, ignore_exception=True)
					ip_mRNA_counts = read_stacked_columns(
						seed_dir, 'RNACounts', 'mRNA_cistron_counts',
						ignore_exception=True)[:, self.new_gene_mRNA_indexes]
					time_data[s] = (time[-1] - time[0]) /2. + time[0]
					monomer_data[s] = [np.mean(ip_monomer_counts, axis=0)]
					mRNA_data[s] = [np.mean(ip_mRNA_counts, axis=0)]
			time = time_data[time_data >= 0] # get rid of the -1 values
			monomer_d = []
			mRNA_d = []
			for m in range(len(interest_proteins)):
				m_mono = []
				m_mRNA = []
				for s in range(len(seeds)):
					if monomer_data[s][m] >= 0:
						m_mono.append(monomer_data[s][m])
					if mRNA_data[s][m] >= 0:
						m_mRNA.append(mRNA_data[s][m])
				m_mono = np.array(m_mono); m_mRNA = np.array(m_mRNA)
				m_mono = np.mean(m_mono, axis=0); m_mRNA = np.mean(m_mRNA,
																   axis=0)
				monomer_d.append(m_mono); mRNA_d.append(m_mRNA)
			time = np.mean(time, axis=0)
			time = np.array(time)
			num_gens = len(time_data) # number of generations averaged over
			gen_data[gen] = (time, monomer_d, mRNA_d, num_gens)

		return gen_data

	def generation_plot(self, control_gen_data, experimental_gen_data,
						experimental_vars):
		# Plotting
		plt.figure(figsize=(8.5, 11))
		# Protein Counts
		plt.subplot(2, 1, 1)
		if len(self.cistron_monomer_ids) == 1:
			time = control_gen_data[:, 0]
			ip_monomer_counts = control_gen_data[:, 1]
			monomer_counts = np.zeros((len(ip_monomer_counts), 1))
			for gen in range(len(ip_monomer_counts)):
				monomer_counts[gen] = ip_monomer_counts[gen][0]
			plt.scatter(time / 60., monomer_counts, label='cntrl var')
			plt.plot(time / 60., monomer_counts, linestyle=':')
			for var in range(len(experimental_vars)):
				time = experimental_gen_data[var][:, 0]
				ip_monomer_counts = experimental_gen_data[var][:, 1]
				monomer_counts = np.zeros((len(ip_monomer_counts), 1))
				for gen in range(len(ip_monomer_counts)):
					monomer_counts[gen] = ip_monomer_counts[gen][0]
				plt.scatter(time / 60., monomer_counts, label='exp var ' +
																 str(var))
				plt.plot(time / 60., monomer_counts, linestyle='-')
			plt.title(f"Protein Counts for {self.cistron_monomer_ids[0]} in the"
					  f"\n control variant"
					  f" and {len(experimental_vars)} experimental variant(s)"
					  f"\n averaged over all seeds at each generation")
		else:
			for m in range(len(self.cistron_monomer_ids)):
				ip_monomer_counts = control_gen_data[:, 1]
				monomer_counts = np.zeros((len(ip_monomer_counts), 1))
				for gen in range(len(ip_monomer_counts)):
					monomer_counts[gen] = ip_monomer_counts[gen][m]
				time = control_gen_data[:, 0]
				name = self.cistron_monomer_ids[m] + ' cntrl var '
				plt.scatter(time / 60., monomer_counts, label=name)
				plt.plot(time / 60., monomer_counts, linestyle=':')
				for var in range(len(experimental_vars)):
					time = experimental_gen_data[var][:, 0]
					ip_monomer_counts = experimental_gen_data[var][:, 1]
					monomer_counts = np.zeros((len(ip_monomer_counts), 1))
					for gen in range(len(ip_monomer_counts)):
						monomer_counts[gen] = ip_monomer_counts[gen][m]
					name = self.cistron_monomer_ids[m] + ' exp var ' + str(var)
					plt.scatter(time / 60., monomer_counts, label=name)
					plt.plot(time / 60., monomer_counts, linestyle='-')
			plt.title(f"Protein Counts for {self.cistron_monomer_ids} in the"
					  f"\n control variant and {len(experimental_vars)} "
					  f"experimental variant(s)"
					  f"\n averaged over all seeds at each generation")

		plt.xlabel("Time (min)")
		plt.ylabel("Protein Counts")
		plt.legend()

		# mRNA Counts
		plt.subplot(2, 1, 2)
		if len(self.cistron_ids) == 1:
			time = control_gen_data[:, 0]
			ip_mRNA_counts = control_gen_data[:, 2]
			mRNA_counts = np.zeros((len(ip_mRNA_counts), 1))
			for gen in range(len(ip_mRNA_counts)):
				mRNA_counts[gen] = ip_mRNA_counts[gen][0]
			plt.scatter(time / 60., mRNA_counts, label='cntrl var')
			plt.plot(time / 60., mRNA_counts, linestyle=':')
			for var in range(len(experimental_vars)):
				time = experimental_gen_data[var][:, 0]
				ip_mRNA_counts = experimental_gen_data[var][:, 2]
				mRNA_counts = np.zeros((len(ip_mRNA_counts), 1))
				for gen in range(len(ip_mRNA_counts)):
					mRNA_counts[gen] = ip_mRNA_counts[gen][0]
				plt.scatter(time / 60., mRNA_counts, label='exp var ' +
															  str(var))
				plt.plot(time / 60., mRNA_counts, linestyle='-')
			plt.title(f"mRNA Counts for {self.cistron_ids[0]} in the control"
					  f" variant"
					  f" and {len(experimental_vars)} experimental variant(s)"
					  f"\n averaged over all seeds at each generation")
		else:
			for m in range(len(self.cistron_ids)):
				time = control_gen_data[:, 0]
				ip_mRNA_counts = control_gen_data[:, 2]
				mRNA_counts = np.zeros((len(ip_mRNA_counts), 1))
				for gen in range(len(ip_mRNA_counts)):
					mRNA_counts[gen] = ip_mRNA_counts[gen][m]
				name = self.cistron_ids[m] + ' cntrl var '
				plt.scatter(time / 60., mRNA_counts, label=name)
				plt.plot(time / 60., mRNA_counts, linestyle=':')
				for var in range(len(experimental_vars)):
					time = experimental_gen_data[var][:, 0]
					ip_mRNA_counts = experimental_gen_data[var][:, 2]
					mRNA_counts = np.zeros((len(ip_mRNA_counts), 1))
					for gen in range(len(ip_mRNA_counts)):
						mRNA_counts[gen] = ip_mRNA_counts[gen][m]
					name = self.cistron_ids[m] + ' exp var ' + str(var)
					plt.scatter(time / 60., mRNA_counts, label=name)
					plt.plot(time / 60., mRNA_counts, linestyle='-')
			plt.title(f"mRNA Counts for {self.cistron_ids} in the \n control"
					  f" variant"
					  f" and {len(experimental_vars)} experimental variant(s)"
					  f"\n averaged over all seeds at each generation")
		plt.xlabel("Time (min)")
		plt.ylabel("Cistron Counts")
		plt.legend()
		plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)

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
		monomer_data_idxs = []
		for protein in interest_proteins:
			monomer_idx = np.where(monomer_sim_data['id'] == protein)
			monomer_idx = monomer_idx[0][0]
			monomer_data_idxs.append(monomer_idx)
		ip_monomer_data = monomer_sim_data[monomer_data_idxs]
		ip_cistron_ids = ip_monomer_data['cistron_id']

		# extract info about the protein(s) from the mRNA/cistron data:
		mRNA_data_idxs = []
		for cistron in ip_cistron_ids:
			mRNA_idx = np.where(mRNA_cistron_sim_data['id'] == cistron)
			mRNA_idx = mRNA_idx[0][0]
			mRNA_data_idxs.append(mRNA_idx)
		self.cistron_ids = ip_cistron_ids
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

		# specifiy the control and variants:
		control_var = all_variants[0]
		experimental_vars = all_variants[1:]
		#experimental_vars = [16, 17, 18, 19, 20]

		# plot for seeds:
		if show_seed_plot == 1:
			self.seed_plot(control_var, experimental_vars)
			plt.subplots_adjust(hspace=0.5, top=0.95, bottom=0.05)
			exportFigure(plt, plotOutDir, plotOutFileName +
						 str(interest_proteins)
					 	+'_all_variants', metadata)

		# plot for generations:
		if show_gen_plot == 1:
			# generate the control variant's data:
			control_gen_data = self.generate_gen_data(control_var)
			# generate the experimental variant's data:
			exp_gen_data = []
			for variant_2 in experimental_vars:
				self.variant_pair = [control_var, variant_2]
				experimental_var = self.variant_pair[1]
				generation_data = self.generate_gen_data(experimental_var)
				exp_gen_data.append(generation_data)
			# make the plot:
			self.generation_plot(control_gen_data, exp_gen_data, experimental_vars)
			exportFigure(plt, plotOutDir, plotOutFileName + str(interest_proteins)+
					 	'_all_variants_over_all_generations', metadata)

		plt.close("all")

if __name__ == '__main__':
	Plot().cli()
