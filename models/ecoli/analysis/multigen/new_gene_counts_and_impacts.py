"""
Plot mRNA and protein counts for new genes across multiple generations, as well
as plots to analyze the impact of new gene expression, including growth rate,
RNAP and ribosome counts, and ppGpp concentration.
"""
# TODO: update file header comment

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile,
				validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()
		sim_dir = cell_paths[0]
		simOutDir = os.path.join(sim_dir, 'simOut')

		# Determine new gene ids
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		mRNA_sim_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_sim_data = sim_data.process.translation.monomer_data.struct_array
		new_gene_mRNA_ids = mRNA_sim_data[mRNA_sim_data['is_new_gene']]['id'].tolist()
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
		assert len(new_gene_monomer_ids) == len(new_gene_mRNA_ids), \
			'number of new gene monomers and mRNAs should be equal'

		# TODO: make a option for this plot that you could run without new genes

		# Extract mRNA indexes for each new gene
		mRNA_counts_reader = TableReader(os.path.join(simOutDir,
													  'RNACounts'))
		mRNA_idx_dict = {rna[:-3]: i for i, rna in enumerate(
			mRNA_counts_reader.readAttribute('mRNA_ids'))}
		new_gene_mRNA_indexes = [mRNA_idx_dict.get(mRNA_id) for mRNA_id in
								 new_gene_mRNA_ids]

		# Extract protein indexes for each new gene
		monomer_counts_reader = TableReader(
			os.path.join(simOutDir, "MonomerCounts"))
		monomer_idx_dict = {monomer: i for i, monomer in
							enumerate(monomer_counts_reader.readAttribute(
								'monomerIds'))}
		new_gene_monomer_indexes = [monomer_idx_dict.get(monomer_id) for
									monomer_id in new_gene_monomer_ids]

		# Extract RNA indexes for each new gene
		rnap_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))
		RNA_idx_dict = {
			rna[:-3]: i for i, rna in
			enumerate(rnap_reader.readAttribute('rnaIds'))}
		new_gene_RNA_indexes = [
			RNA_idx_dict.get(mRNA_id) for mRNA_id in new_gene_mRNA_ids]

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', ignore_exception=True)
		time_no_first = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True, ignore_exception=True)
		(new_gene_monomer_counts,) = read_stacked_bulk_molecules(
			cell_paths, new_gene_monomer_ids, ignore_exception=True)
		all_mRNA_stacked_counts = read_stacked_columns(
			cell_paths, 'RNACounts', 'mRNA_counts', ignore_exception=True)
		new_gene_mRNA_counts = all_mRNA_stacked_counts[:,new_gene_mRNA_indexes]
		new_gene_copy_numbers = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'gene_copy_number',
			ignore_exception=True)[:,new_gene_RNA_indexes]

		plot_suffixes = ["", "_standard_axes_y", "_standard_axes_both"]
		standard_xlim = (0,2000)
		total_plots = 9 # TODO Modularize and get rid of this magic number

		for i in range(len(plot_suffixes)):

			plot_suffix = plot_suffixes[i]

			# Plotting
			plt.figure(figsize = (8.5, 22))
			plot_num = 1

			# Growth Rate
			ax1 = plt.subplot(total_plots, 1, plot_num)
			growth_rate = np.ravel(read_stacked_columns(
				cell_paths, "Mass", "instantaneous_growth_rate",
				ignore_exception=True))
			moving_window = min(300, len(growth_rate))
			convolution_array = (np.ones(moving_window) / moving_window)
			growth_rate_convolved = np.convolve(
				convolution_array, growth_rate, mode='same')
			plt.plot(time.flatten() / 60., growth_rate_convolved)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(-.0001,.0005)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Growth Rate")
			plt.title("Growth Rate")
			plot_num += 1

			# Mass
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			mass = read_stacked_columns(
				cell_paths, "Mass", "cellMass", ignore_exception=True)
			plt.plot(time / 60., mass)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim(0,4000)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Mass (fg)")
			plt.title("Cell Mass")
			plot_num += 1

			# Gene Copy Number
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if len(new_gene_mRNA_ids) == 1:
				plt.plot(time / 60., new_gene_copy_numbers,
						 label=new_gene_mRNA_ids[0])
			else:
				for r in range(len(new_gene_mRNA_ids)):
					plt.plot(time / 60., new_gene_copy_numbers[:,r],
							 label = new_gene_mRNA_ids[r])
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-1,10))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Gene Copy Number")
			plt.title("New Gene Copy Number")
			plt.legend()
			plot_num += 1

			# mRNA Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if plot_suffix == "":
				if len(new_gene_mRNA_ids) == 1:
					plt.plot(time / 60., new_gene_mRNA_counts,
							 label=new_gene_mRNA_ids[0])
				else:
					for r in range(len(new_gene_mRNA_ids)):
						plt.plot(time / 60., new_gene_mRNA_counts[:,r],
								 label = new_gene_mRNA_ids[r])
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# plot on log scale instead
				if len(new_gene_mRNA_ids) == 1:
					plt.plot(time / 60., np.log10(new_gene_mRNA_counts + 1),
							 label=new_gene_mRNA_ids[0])
				else:
					for r in range(len(new_gene_mRNA_ids)):
						plt.plot(time / 60., np.log10(new_gene_mRNA_counts[:,r] + 1),
								 label = new_gene_mRNA_ids[r])
				plt.ylim((-1,4.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA Counts")
			plt.title("New Gene mRNA Counts")
			plt.legend()
			plot_num += 1

			# Protein Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			if plot_suffix == "":
				if len(new_gene_monomer_ids) == 1:
					plt.plot(time / 60., new_gene_monomer_counts,
							 label=new_gene_monomer_ids[0])
				else:
					for m in range(len(new_gene_monomer_ids)):
						plt.plot(time / 60., new_gene_monomer_counts[:,m],
								 label = new_gene_monomer_ids[m])
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				# plot on log scale instead
				if len(new_gene_monomer_ids) == 1:
					plt.plot(time / 60., np.log10(new_gene_monomer_counts + 1),
							 label=new_gene_monomer_ids[0])
				else:
					for m in range(len(new_gene_monomer_ids)):
						plt.plot(time / 60., np.log10(new_gene_monomer_counts[:,m] + 1),
								 label = new_gene_monomer_ids[m])
				plt.ylim((-1,7.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Protein Counts")
			plt.title("New Gene Protein Counts")
			plt.legend()
			plot_num += 1

			# ppGpp
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			ppGpp_concentration = read_stacked_columns(
				cell_paths, "GrowthLimits", "ppgpp_conc", remove_first=True,
				ignore_exception=True)
			plt.plot(time_no_first / 60., ppGpp_concentration)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,150))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("ppGpp Concentration (uM)")
			plt.title("ppGpp Concentration")
			plot_num += 1

			# Total RNAP Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			rnap_id = [sim_data.molecule_ids.full_RNAP]
			(inactive_rnap_counts,) = read_stacked_bulk_molecules(
				cell_paths, (rnap_id,), ignore_exception=True)
			# Active
			uniqueMoleculeCounts = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			active_rnap_index = uniqueMoleculeCounts.readAttribute(
				"uniqueMoleculeIds").index('active_RNAP')
			active_rnap_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts',
				ignore_exception=True)[:, active_rnap_index]
			# Total
			total_rnap_counts = inactive_rnap_counts + active_rnap_counts
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,10000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.plot(time / 60., total_rnap_counts)
			plt.xlabel("Time (min)")
			plt.ylabel("Total RNAP Counts")
			plt.title("Total RNAP Counts")
			plot_num += 1

			# Total Ribosome Counts
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			# Inactive
			complex_id_30s = [sim_data.molecule_ids.s30_full_complex]
			complex_id_50s = [sim_data.molecule_ids.s50_full_complex]
			(complex_counts_30s, complex_counts_50s) = read_stacked_bulk_molecules(
				cell_paths, (complex_id_30s, complex_id_50s), ignore_exception=True)
			inactive_ribosome_counts = np.minimum(
				complex_counts_30s, complex_counts_50s)
			# Active
			unique_molecule_counts_table = TableReader(
				os.path.join(simOutDir, "UniqueMoleculeCounts"))
			ribosome_index = unique_molecule_counts_table.readAttribute(
				"uniqueMoleculeIds").index('active_ribosome')
			active_ribosome_counts = read_stacked_columns(
				cell_paths, 'UniqueMoleculeCounts',
				'uniqueMoleculeCounts', ignore_exception=True)[:, ribosome_index]
			# Total
			total_ribosome_counts = active_ribosome_counts + inactive_ribosome_counts
			plt.plot(time / 60., total_ribosome_counts)
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((0,40000))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("Total Ribosome Counts")
			plt.title("Total Ribosome Counts")
			plot_num += 1

			# TODO Yield

			# New Gene RNA Synthesis Probability
			plt.subplot(total_plots, 1, plot_num, sharex=ax1)
			new_gene_target_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'target_rna_synth_prob',
				ignore_exception=True)[:, new_gene_RNA_indexes]
			new_gene_actual_rna_synth_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
				ignore_exception=True)[:, new_gene_RNA_indexes]

			if len(new_gene_mRNA_ids) == 1:
				plt.plot(time / 60., new_gene_target_rna_synth_prob,
						 label=new_gene_mRNA_ids[0] + ": Target")
				plt.plot(time / 60., new_gene_actual_rna_synth_prob,
						 label=new_gene_mRNA_ids[0] + ": Actual")
			else:
				for r in range(len(new_gene_mRNA_ids)):
					plt.plot(time / 60., new_gene_target_rna_synth_prob[:,r],
							 label = new_gene_mRNA_ids[r] + ": Target")
					plt.plot(time / 60., new_gene_actual_rna_synth_prob[:, r],
							 label=new_gene_mRNA_ids[r] + ": Actual")
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_y":
				plt.ylim((-0.1,0.5))
			if plot_suffix == "_standard_axes_both" or plot_suffix == "_standard_axes_x":
				plt.xlim(standard_xlim)
			plt.xlabel("Time (min)")
			plt.ylabel("RNA Synthesis Probability")
			plt.title("New Gene RNA Synthesis Probability")
			plt.legend()
			plot_num += 1

			# TODO Target Protein Init Prob

			# TODO Actual Protein Init Prob

			plt.subplots_adjust(hspace = 0.7, top = 0.95, bottom = 0.05)
			exportFigure(plt, plotOutDir, plotOutFileName + plot_suffix, metadata)
			plt.close("all")

if __name__ == '__main__':
	Plot().cli()
