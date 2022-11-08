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
	def make_plot(self, axs, time, var, titles, ma=None, ylim=None, c=COLORS):
		if ma:
			var = self.moving_average(var, ma)
			time = self.moving_average(time, ma)
			for (i, mv_avg) in enumerate(var):
				if len(np.shape(mv_avg)) > 1:
					for (j, subvar) in enumerate(mv_avg):
						axs[j].plot(time[i], subvar, c=c[i], label="Mov avg window %fs" % int(ma[i] * self.avg_timestep))
						if i==0 and titles:
							axs[j].set_title(titles[j], pad=LABELPAD, fontsize=FONTSIZE)
						if ylim:
							axs[j].set_ylim(ylim)
						if i == len(ma) - 1:
							axs[j].legend(fontsize=LEGEND_FONT)
				else:
					axs.plot(time[i], mv_avg, c=c[i], label="Mov avg window %fs" % int(ma[i] * self.avg_timestep))
					if titles:
						axs.set_title(titles, pad=LABELPAD, fontsize=FONTSIZE)
					if ylim:
						axs.set_ylim(ylim)
					if i == len(ma) - 1:
						axs.legend(fontsize=LEGEND_FONT)
		else:
			if len(np.shape(var)) > 1:
				var = var.T
				for (j, subvar) in enumerate(var):
					axs[j].plot(time, subvar, c=COLORS[0])
					if titles:
						axs[j].set_title(titles[j], pad=LABELPAD, fontsize=FONTSIZE)
					if ylim:
						axs[j].set_ylim(ylim)
			else:
				axs.plot(time, var, c=COLORS[0])
				if titles:
					axs.set_title(titles, pad=LABELPAD, fontsize=FONTSIZE)
				if ylim:
					axs.set_ylim(ylim)

	def moving_average(self, timeseries, ma):
		# Takes the moving average to show the timescale of dynamics
		convs = [np.ones(window) / window for window in ma]
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
		timestep = main_reader.readColumn("timeStepSec")
		self.avg_timestep = np.mean(timestep)

		# Transcriptional probabilities
		rna_synth_prob = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		actual_synth_prob = rna_synth_prob.readColumn("actual_rna_synth_prob")
		target_synth_prob = rna_synth_prob.readColumn("target_rna_synth_prob")
		is_overcrowded = rna_synth_prob.readColumn("tu_is_overcrowded")
		rna_ids = np.array(rna_synth_prob.readAttribute("rnaIds"))
		rRNA_idxs = [np.where(rna_ids == rRNA_id)[0][0] for rRNA_id in rRNA_ids]
		rRNA_actual_synth_prob = actual_synth_prob[:, rRNA_idxs]
		rRNA_target_synth_prob = target_synth_prob[:, rRNA_idxs]
		rRNA_is_overcrowded_fraction = np.mean(is_overcrowded[:, rRNA_idxs], axis=0)
		total_rRNA_is_overcrowded_fraction = np.mean([is_overcrowded[i, rRNA_idxs].any()
														 for i in range(np.shape(rRNA_is_overcrowded_fraction)[0])])
		total_rRNA_actual_synth_prob = np.sum(rRNA_actual_synth_prob, axis=1)
		total_rRNA_target_synth_prob = np.sum(rRNA_target_synth_prob, axis=1)

		# Partial rRNAs, or RNAPs on rrn genes
		rna_counts = TableReader(os.path.join(simOutDir, "mRnaCounts"))
		partial_rRNA_counts = rna_counts.readColumn("partial_rRNA_counts")
		rna_ids = np.array(rna_counts.readAttribute("rRNA_ids"))
		rRNA_idxs = [np.where(rna_ids == rRNA_id)[0][0] for rRNA_id in rRNA_ids]
		RNAP_on_rRNA_counts = partial_rRNA_counts[:, rRNA_idxs]
		total_RNAP_on_rRNA_counts = np.sum(RNAP_on_rRNA_counts, axis=1)

		# Transcription initiation events
		rnap_data = TableReader(os.path.join(simOutDir, "RnapData"))
		rna_init_event = rnap_data.readColumn("rnaInitEvent")
		rna_ids = np.array(rnap_data.readAttribute("rnaIds"))
		rRNA_idxs = [np.where(rna_ids == rRNA_id)[0][0] for rRNA_id in rRNA_ids]
		rRNA_init_events = rna_init_event[:, rRNA_idxs]
		rRNA_init_per_sec = np.array([rna_init / timestep for rna_init in rRNA_init_events.T]).T

		ribosome_data = TableReader(os.path.join(simOutDir, "RibosomeData"))
		total_rRNA_initiated = ribosome_data.readColumn("total_rRNA_initiated", squeeze=True)
		total_rRNA_init_per_sec = total_rRNA_initiated / timestep
		rRNA5S_initiated = ribosome_data.readColumn("rRNA5S_initiated", squeeze=True)
		rRNA5S_init_per_sec = rRNA5S_initiated / timestep
		rRNA16S_initiated = ribosome_data.readColumn("rRNA16S_initiated", squeeze=True)
		rRNA16S_init_per_sec = rRNA16S_initiated / timestep
		rRNA23S_initiated = ribosome_data.readColumn("rRNA23S_initiated", squeeze=True)
		rRNA23S_init_per_sec = rRNA23S_initiated / timestep

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
		rRNA_initiation_per_gene_dosage = rRNA_init_per_sec / rRNA_gene_dosage
		total_rRNA_initiation_per_gene_dosage = total_rRNA_init_per_sec / total_rRNA_gene_dosage

		# Counts of various ribosome-related molecules and RNAPs
		ribosome_30S_id = [sim_data.molecule_ids.s30_full_complex]
		ribosome_50S_id = [sim_data.molecule_ids.s50_full_complex]
		rRNA_5S_ids = sim_data.molecule_groups.s50_5s_rRNA
		rRNA_16S_ids = sim_data.molecule_groups.s30_16s_rRNA
		rRNA_23S_ids = sim_data.molecule_groups.s50_23s_rRNA
		inactive_rnap_id = [sim_data.molecule_ids.full_RNAP]

		(ribosome_30S_counts, ribosome_50S_counts, rRNA_5S_counts,
		rRNA_16S_counts, rRNA_23S_counts, inactive_RNAP_counts) = read_bulk_molecule_counts(
			simOutDir, (ribosome_30S_id, ribosome_50S_id, rRNA_5S_ids,
			rRNA_16S_ids, rRNA_23S_ids, inactive_rnap_id))
		rRNA_5S_counts = np.sum(rRNA_5S_counts, axis=1)
		rRNA_16S_counts = np.sum(rRNA_16S_counts, axis=1)
		rRNA_23S_counts = np.sum(rRNA_23S_counts, axis=1)

		unique_molecules = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		unique_molecule_ids = unique_molecules.readAttribute("uniqueMoleculeIds")
		ribosome_idx = unique_molecule_ids.index('active_ribosome')
		rnap_idx = unique_molecule_ids.index('active_RNAP')
		unique_molecule_counts = unique_molecules.readColumn("uniqueMoleculeCounts")
		active_ribosome = unique_molecule_counts[:, ribosome_idx]
		active_RNAP_counts = unique_molecule_counts[:, rnap_idx]
		total_RNAP_counts = active_RNAP_counts + inactive_RNAP_counts
		total_30S_counts = ribosome_30S_counts + active_ribosome
		total_50S_counts = ribosome_50S_counts + active_ribosome

		RNAP_active_fraction_on_rRNA = np.array([RNAP_counts / active_RNAP_counts for RNAP_counts in RNAP_on_rRNA_counts.T]).T
		RNAP_total_fraction_on_rRNA = np.array([RNAP_counts / total_RNAP_counts for RNAP_counts in RNAP_on_rRNA_counts.T]).T
		total_RNAP_active_fraction_on_rRNA = total_RNAP_on_rRNA_counts / active_RNAP_counts
		total_RNAP_total_fraction_on_rRNA = total_RNAP_on_rRNA_counts / total_RNAP_counts

		# Get rRNA common names for operons
		rRNA_operon_to_name = {"TU0-1181[c]": "rrnA", "TU0-1182[c]": "rrnB",
			"TU0-1183[c]": "rrnC", "TU0-1186[c]": "rrnE", "TU0-1187[c]": "rrnG",
			"TU0-1189[c]": "rrnH", "TU0-1191[c]": "rrnD"
			}

		## Make figure
		# Ncols is number of rRNA TUs, 1 for total, and 3 for 5S/16S/23S
		ncols = num_rRNAs + 1 + 3
		nrows=7

		fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8 * ncols, 8 * nrows))
		# Plot operon-level rRNA initiation events
		axs[0, 0].set_ylabel("Initiations per second", labelpad=LABELPAD)
		rRNA_names = [rRNA_operon_to_name[rRNA] if rRNA in rRNA_operon_to_name
			else rRNA for rRNA in rRNA_ids]
		self.make_plot(axs=[axs[0, i] for i in range(num_rRNAs)], time=time, var=rRNA_init_per_sec,
					   titles=rRNA_names, ma=MA_WINDOWS, ylim=[0, 4])
		self.make_plot(axs=axs[0, num_rRNAs], time=time, var=total_rRNA_init_per_sec,
					   titles="Total rRNAs", ma=MA_WINDOWS, ylim=[0, 15])
		for i in range(ncols):
			axs[0, i].set_xlabel("Time (min)", labelpad=LABELPAD, fontsize=FONTSIZE)

		# Plot gene dosage of each operon
		axs[1, 0].set_ylabel("Gene dosage", labelpad=LABELPAD, fontsize=FONTSIZE)
		self.make_plot(axs=[axs[1, i] for i in range(num_rRNAs)], time=time, var=rRNA_gene_dosage,
					   titles=None, ylim=[0, 5])
		self.make_plot(axs=axs[1, num_rRNAs], time=time, var=total_rRNA_gene_dosage,
					   titles=None, ylim=[0, 5])

		# Plot the rRNA initiation events per gene dosage
		axs[2, 0].set_ylabel("Initiations per second per gene dosage", labelpad=LABELPAD, fontsize=FONTSIZE)
		self.make_plot(axs=[axs[2,i] for i in range(num_rRNAs)], time=time, var=rRNA_initiation_per_gene_dosage,
					   titles=None, ma=MA_WINDOWS, ylim=[0, 1])
		self.make_plot(axs=axs[2, num_rRNAs], time=time, var=total_rRNA_initiation_per_gene_dosage,
					   titles=None, ma=MA_WINDOWS, ylim=[0, 1])

		# Plot effective initiation events of 5S, 16S, and 23S rRNA
		self.make_plot(axs=axs[0, num_rRNAs+1], time=time, var=rRNA5S_init_per_sec,
					   titles="5S rRNA inits", ma=MA_WINDOWS, ylim=[0, 15])
		self.make_plot(axs=axs[0, num_rRNAs+2], time=time, var=rRNA16S_init_per_sec,
					   titles="16S rRNA inits", ma=MA_WINDOWS, ylim=[0, 15])
		self.make_plot(axs=axs[0, num_rRNAs+3], time=time, var=rRNA23S_init_per_sec,
					   titles="23S rRNA inits", ma=MA_WINDOWS, ylim=[0, 15])

		# Plot the target and actual synthesis probabilities of rRNA operons
		axs[3, 0].set_ylabel("Target and actual rRNA synthesis probabilities")
		self.make_plot(axs=[axs[3, i] for i in range(num_rRNAs)], time=time, var=rRNA_actual_synth_prob,
					   titles=None, ma=[200], c=['r'])
		self.make_plot(axs=[axs[3, i] for i in range(num_rRNAs)], time=time, var=rRNA_target_synth_prob,
					   titles=None, ma=[200], c=['g'],  ylim=[0, 0.04])
		for i in range(num_rRNAs):
			axs[3, i].text(np.mean(time), axs[3, i].get_ylim()[1],
						"Overcrowded time fraction: %f"%rRNA_is_overcrowded_fraction[i])
		self.make_plot(axs=axs[3, num_rRNAs], time=time, var=total_rRNA_actual_synth_prob,
					   titles=None, ma=[200], c=['r'])
		self.make_plot(axs=axs[3, num_rRNAs], time=time, var=total_rRNA_target_synth_prob,
					   titles=None, ma=[200], c=['g'], ylim=(0, None))
		axs[3, num_rRNAs].text(np.mean(time), axs[3, num_rRNAs].get_ylim()[1],
							"Overcrowded time fraction: %f"%total_rRNA_is_overcrowded_fraction)

		# Plot the number of RNAPs on rrn genes
		axs[4, 0].set_ylabel("Counts of RNAPs on rrn genes")
		self.make_plot(axs=[axs[4, i] for i in range(num_rRNAs)], time=time, var=RNAP_on_rRNA_counts,
					   titles=None, ma=MA_WINDOWS, c=['r', 'b'], ylim=[0, 200])
		self.make_plot(axs=axs[4, num_rRNAs], time=time, var=total_RNAP_on_rRNA_counts,
					   titles=None, ma=MA_WINDOWS, c=['r', 'b'], ylim=[0, 1000])

		# Plot the fraction of active or total RNAPs on rrn genes
		axs[5, 0].set_ylabel("Fraction of total/active RNAPs on rrn genes")
		self.make_plot(axs=[axs[5, i] for i in range(num_rRNAs)], time=time,
					   var=RNAP_total_fraction_on_rRNA,
					   titles=None, ma=MA_WINDOWS, c=['g', 'k'], ylim=[0, 0.2])
		self.make_plot(axs=[axs[5, i] for i in range(num_rRNAs)], time=time,
					   var=RNAP_active_fraction_on_rRNA,
					   titles=None, ma=MA_WINDOWS, c=['r','b'])
		self.make_plot(axs=axs[5, num_rRNAs], time=time, var=total_RNAP_active_fraction_on_rRNA,
					   titles=None, ma=MA_WINDOWS, c=['m', 'y'],  ylim=[0, None])
		self.make_plot(axs=axs[5, num_rRNAs], time=time, var=total_RNAP_total_fraction_on_rRNA,
					   titles=None, ma=MA_WINDOWS, c=['g', 'k'])

		# Plot the amount of (active or inactive) ribosomal subunits
		axs[6, 0].set_ylabel("Counts", labelpad=LABELPAD, fontsize=FONTSIZE)
		self.make_plot(axs=axs[6, 0], time=time, var=total_30S_counts,
					   titles="Total 30S", ma=MA_WINDOWS,  ylim=(0, None))
		self.make_plot(axs=axs[6, 1], time=time, var=total_50S_counts,
					   titles="Total 50S", ma=MA_WINDOWS,  ylim=(0, None))

		# Plot the amount of uncomplexed rRNAs
		self.make_plot(axs=axs[6, 2], time=time, var=rRNA_5S_counts,
					   titles="Free 5S rRNA", ma=MA_WINDOWS,  ylim=(0, None))
		self.make_plot(axs=axs[6, 3], time=time, var=rRNA_16S_counts,
					   titles="Free 16S rRNA", ma=MA_WINDOWS,  ylim=(0, None))
		self.make_plot(axs=axs[6, 4], time=time, var=rRNA_23S_counts,
					   titles="Free 23S rRNA", ma=MA_WINDOWS,  ylim=(0, None))

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
