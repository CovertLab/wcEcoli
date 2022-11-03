"""
Plot expression and related properties for rRNA operons
"""

from __future__ import absolute_import, division, print_function

import os
import numpy as np

from matplotlib import pyplot as plt
import pickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot

MA_WINDOWS = [200, 400] # moving average window to plot timeseries, in timesteps
COLORS = ['r', 'b']
LABELPAD = 10
FONTSIZE = 15
LEGEND_FONT = 10

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	# TODO: add listener that tracks partial ribosome formation, in order to get
	#  how many RNAPs are on rRNAs at this time point
	def make_plot(self, axs, time, var, titles, ma, ylim=None):
		if ma:
			var = self.moving_average(var)
			time = self.moving_average(time)
			for (i, mv_avg) in enumerate(var):
				if len(np.shape(mv_avg)) > 1:
					for (j, subvar) in enumerate(mv_avg):
						axs[j].plot(time[i], subvar, c=COLORS[i], label="Mov avg window %f timesteps" % MA_WINDOWS[i])
						if i==0 and titles:
							axs[j].set_title(titles[j], pad=LABELPAD, fontsize=FONTSIZE)
						if ylim:
							axs[j].set_ylim(ylim)
						if i == len(MA_WINDOWS) - 1:
							axs[j].legend(fontsize=LEGEND_FONT)
				else:
					axs.plot(time[i], mv_avg, c=COLORS[i], label="Mov avg window %f timesteps" % MA_WINDOWS[i])
					if titles:
						axs.set_title(titles, pad=LABELPAD, fontsize=FONTSIZE)
					if ylim:
						axs.set_ylim(ylim)
					if i == len(MA_WINDOWS) - 1:
						axs.legend(fontsize=LEGEND_FONT)
		else:
			if len(np.shape(var)) > 1:
				var = var.T
				for (j, subvar) in enumerate(var):
					axs[j].plot(time, subvar, c=COLORS[0])
					if titles:
						axs[j].set_title(titles[j], pad=LABELPAD, fontsize=FONTSIZE)
					if ylim:
						axs.set_ylim(ylim)
			else:
				axs.plot(time, var, c=COLORS[0])
				if titles:
					axs.set_title(titles, pad=LABELPAD, fontsize=FONTSIZE)
				if ylim:
					axs.set_ylim(ylim)

	def moving_average(self, timeseries):
		# Takes the moving average to show the timescale of dynamics
		convs = [np.ones(window) / window for window in MA_WINDOWS]
		if len(np.shape(timeseries)) > 1:
			timeseries = timeseries.T
			moving_averages = []
			for conv in convs:
				mov_avg = [np.convolve(series, conv, mode='valid')
					for series in timeseries]
				moving_averages.append(np.array(mov_avg))
			return moving_averages
		else:
			return [np.convolve(timeseries, conv, mode='valid') for conv in convs]

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = pickle.load(open(simDataFile, "rb"))
		transcription = sim_data.process.transcription
		rRNA_ids = transcription.rna_data['id'][transcription.rna_data['is_rRNA']]
		num_rRNAs = len(rRNA_ids)

		# Time
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initial_time = main_reader.readAttribute("initialTime")
		time = (main_reader.readColumn("time") - initial_time)/60

		# Transcription initiation events
		rnap_data = TableReader(os.path.join(simOutDir, "RnapData"))
		rna_init_event = rnap_data.readColumn("rnaInitEvent")
		rna_ids = np.array(rnap_data.readAttribute("rnaIds"))
		rRNA_idxs = [np.where(rna_ids == rRNA_id)[0][0] for rRNA_id in rRNA_ids]
		rRNA_init_events = rna_init_event[:, rRNA_idxs]

		ribosome_data = TableReader(os.path.join(simOutDir, "RibosomeData"))
		total_rRNA_initiated = ribosome_data.readColumn("total_rRNA_initiated", squeeze=True)
		rRNA5S_initiated = ribosome_data.readColumn("rRNA5S_initiated", squeeze=True)
		rRNA16S_initiated = ribosome_data.readColumn("rRNA16S_initiated", squeeze=True)
		rRNA23S_initiated = ribosome_data.readColumn("rRNA23S_initiated", squeeze=True)

		# For each TU, get a cistron from it and calculate its gene dosage
		rna_synth_prob = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		gene_ids = np.array(rna_synth_prob.readAttribute("gene_ids"))
		gene_copy_number = rna_synth_prob.readColumn("gene_copy_number")
		rRNA_cistron_idxs = np.array([transcription.rna_id_to_cistron_indexes(rRNA_id)[0]
			for rRNA_id in rRNA_ids])
		rRNA_cistron_ids = transcription.cistron_data['gene_id'][rRNA_cistron_idxs]
		rRNA_gene_idxs = []
		for cistron in rRNA_cistron_ids:
			rRNA_gene_idxs.extend(np.where(gene_ids == cistron)[0])
		rRNA_gene_dosage = gene_copy_number[:, rRNA_gene_idxs]
		total_rRNA_gene_dosage = [np.sum(rRNA, axis=0) for rRNA in rRNA_gene_dosage]

		# Initiations per gene dosage
		rRNA_initiation_per_gene_dosage = rRNA_init_events / rRNA_gene_dosage
		total_rRNA_initiation_per_gene_dosage = total_rRNA_initiated / total_rRNA_gene_dosage

		# Counts of various ribosome-related molecules
		ribosome_30S_id = [sim_data.molecule_ids.s30_full_complex]
		ribosome_50S_id = [sim_data.molecule_ids.s50_full_complex]
		rRNA_5S_ids = sim_data.molecule_groups.s50_5s_rRNA
		rRNA_16S_ids = sim_data.molecule_groups.s30_16s_rRNA
		rRNA_23S_ids = sim_data.molecule_groups.s50_23s_rRNA

		(ribosome_30S_counts, ribosome_50S_counts, rRNA_5S_counts,
		rRNA_16S_counts, rRNA_23S_counts) = read_bulk_molecule_counts(
			simOutDir, (ribosome_30S_id, ribosome_50S_id, rRNA_5S_ids,
			rRNA_16S_ids, rRNA_23S_ids))
		rRNA_5S_counts = np.sum(rRNA_5S_counts, axis=1)
		rRNA_16S_counts = np.sum(rRNA_16S_counts, axis=1)
		rRNA_23S_counts = np.sum(rRNA_23S_counts, axis=1)

		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index('active_ribosome')
		active_ribosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
		total_30S_counts = ribosome_30S_counts + active_ribosome
		total_50S_counts = ribosome_50S_counts + active_ribosome

		## Make figure
		# Ncols is number of rRNA TUs, 1 for total, and 3 for 5S/16S/23S
		ncols = num_rRNAs + 1 + 3
		nrows=4

		fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 5 * nrows))
		# Plot operon-level rRNA initiation events
		axs[0, 0].set_ylabel("Initiations per second", labelpad=LABELPAD)
		self.make_plot(axs=[axs[0, i] for i in range(num_rRNAs)], time=time, var=rRNA_init_events,
					   titles=[rRNA for rRNA in rRNA_ids], ma=True, ylim=[0, 2])
		self.make_plot(axs=axs[0, num_rRNAs], time=time, var=total_rRNA_initiated,
					   titles="Total rRNAs", ma=True, ylim=[4, 8])
		for i in range(ncols):
			axs[0, i].set_xlabel("Time (min)", labelpad=LABELPAD, fontsize=FONTSIZE)

		# Plot gene dosage of each operon
		axs[1, 0].set_ylabel("Gene dosage", labelpad=LABELPAD, fontsize=FONTSIZE)
		self.make_plot(axs=[axs[1, i] for i in range(num_rRNAs)], time=time, var=rRNA_gene_dosage,
					   titles=None, ma=False)
		self.make_plot(axs=axs[1, num_rRNAs], time=time, var=total_rRNA_gene_dosage,
					   titles=None, ma=False)

		# Plot the rRNA initiation events per gene dosage
		axs[2, 0].set_ylabel("Initiations per second per gene dosage", labelpad=LABELPAD, fontsize=FONTSIZE)
		self.make_plot(axs=[axs[2,i] for i in range(num_rRNAs)], time=time, var=rRNA_initiation_per_gene_dosage,
					   titles=None, ma=True, ylim=[0, 1])
		self.make_plot(axs=axs[2, num_rRNAs], time=time, var=total_rRNA_initiation_per_gene_dosage,
					   titles=None, ma=True, ylim=[0, 1])

		# Plot effective initiation events of 5S, 16S, and 23S rRNA
		self.make_plot(axs=axs[0, num_rRNAs+1], time=time, var=rRNA5S_initiated,
					   titles="5S rRNA inits", ma=True, ylim=[4, 8])
		self.make_plot(axs=axs[0, num_rRNAs+2], time=time, var=rRNA16S_initiated,
					   titles="16S rRNA inits", ma=True, ylim=[4, 8])
		self.make_plot(axs=axs[0, num_rRNAs+3], time=time, var=rRNA23S_initiated,
					   titles="23S rRNA inits", ma=True, ylim=[4, 8])

		# Plot the amount of (active or inactive) ribosomal subunits
		axs[3, 0].set_ylabel("Counts", labelpad=LABELPAD, fontsize=FONTSIZE)
		self.make_plot(axs=axs[3, 0], time=time, var=total_30S_counts,
					   titles="Total 30S", ma=True)
		self.make_plot(axs=axs[3, 1], time=time, var=total_50S_counts,
					   titles="Total 50S", ma=True)

		# Plot the amount of uncomplexed rRNAs
		self.make_plot(axs=axs[3, 2], time=time, var=rRNA_5S_counts,
					   titles="Free 5S rRNA", ma=True)
		self.make_plot(axs=axs[3, 3], time=time, var=rRNA_16S_counts,
					   titles="Free 16S rRNA", ma=True)
		self.make_plot(axs=axs[3, 4], time=time, var=rRNA_23S_counts,
					   titles="Free 23S rRNA", ma=True)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
