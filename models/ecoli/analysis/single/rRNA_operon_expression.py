"""
Plot various expression-related properties for rRNA operons
"""

import os
import numpy as np

from matplotlib import pyplot as plt
import pickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot

LABELPAD = 10
FONTSIZE = 15
COLORS = ['C0', 'C1', 'C2', 'royalblue']

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = self.read_pickle_file(simDataFile)
		transcription = sim_data.process.transcription
		is_rRNA_cistron = transcription.cistron_data['is_rRNA']
		cistron_rRNA_ids = transcription.cistron_data['id'][is_rRNA_cistron]
		gene_rRNA_ids = transcription.cistron_data['gene_id'][is_rRNA_cistron]

		molecule_groups = sim_data.molecule_groups
		cistrons_of_rRNA_TUs = [
			getattr(molecule_groups, operon) for operon in
			molecule_groups.rrn_operons]
		rRNA_operon_names = molecule_groups.rrn_operons
		rRNA_cistron_types = ["5S", "23S", "16S"]

		rRNA_TU_to_cistron_mapping_matrix = np.array([
			np.isin(cistron_rRNA_ids, cistrons)
			for cistrons in cistrons_of_rRNA_TUs])
		rRNA_cistron_to_type_mapping_matrix = np.array([
			transcription.cistron_data['is_5S_rRNA'][is_rRNA_cistron],
			transcription.cistron_data['is_16S_rRNA'][is_rRNA_cistron],
			transcription.cistron_data['is_23S_rRNA'][is_rRNA_cistron]])

		# Time
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initial_time = main_reader.readAttribute("initialTime")
		time = (main_reader.readColumn("time") - initial_time)/60
		timestep = main_reader.readColumn("timeStepSec")

		# Transcriptional probabilities
		rna_synth_prob = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		actual_cistron_synth_prob = rna_synth_prob.readColumn("actual_rna_synth_prob_per_cistron")
		target_cistron_synth_prob = rna_synth_prob.readColumn("target_rna_synth_prob_per_cistron")
		cistron_id_to_index = {cistron_id: i for (i, cistron_id) in enumerate(
			rna_synth_prob.readAttribute("cistron_ids"))}
		rRNA_idxs = np.array([cistron_id_to_index[rRNA_id] for rRNA_id in cistron_rRNA_ids])
		rRNA_actual_synth_prob = actual_cistron_synth_prob[:, rRNA_idxs]
		rRNA_target_synth_prob = target_cistron_synth_prob[:, rRNA_idxs]
		grouped_rRNA_actual_synth_prob = np.array(
			rRNA_cistron_to_type_mapping_matrix @ rRNA_actual_synth_prob.T)
		grouped_rRNA_target_synth_prob = np.array(
			rRNA_cistron_to_type_mapping_matrix @ rRNA_target_synth_prob.T)

		# Actual transcription initiation events
		rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))
		actual_rna_init_per_cistron = rnap_data_reader.readColumn("rna_init_event_per_cistron")
		rRNA_actual_rna_init = actual_rna_init_per_cistron[:, rRNA_idxs]
		rRNA_actual_rna_init_rate = np.array(
			[rna_init / timestep for rna_init in rRNA_actual_rna_init.T]).T
		grouped_rRNA_actual_rna_init_rate = np.array(
			rRNA_cistron_to_type_mapping_matrix @ rRNA_actual_rna_init_rate.T)

		# Expected transcription initiation events
		expected_rna_init_per_cistron = rna_synth_prob.readColumn("expected_rna_init_per_cistron")
		rRNA_expected_rna_init = expected_rna_init_per_cistron[:, rRNA_idxs]
		rRNA_expected_rna_init_rate = np.array([rna_init / timestep for rna_init in rRNA_expected_rna_init.T]).T
		grouped_rRNA_expected_rna_init_rate = np.array(
			rRNA_cistron_to_type_mapping_matrix @ rRNA_expected_rna_init_rate.T)

		# Gene dosage
		gene_copy_number = rna_synth_prob.readColumn("gene_copy_number")
		gene_id_to_index = {gene_id: i for (i, gene_id) in enumerate(
			rna_synth_prob.readAttribute("gene_ids"))}
		rRNA_gene_idxs = np.array([gene_id_to_index[gene_rRNA_id] for gene_rRNA_id in gene_rRNA_ids])
		rRNA_gene_dosages = gene_copy_number[:, rRNA_gene_idxs]
		grouped_rRNA_gene_dosages = np.array(rRNA_cistron_to_type_mapping_matrix @ rRNA_gene_dosages.T)

		# Initiations per gene dosage
		rRNA_initiation_per_gene_dosage = rRNA_expected_rna_init_rate / rRNA_gene_dosages
		grouped_rRNA_initiation_per_gene_dosage = np.array(
			rRNA_cistron_to_type_mapping_matrix @ rRNA_initiation_per_gene_dosage.T)

		# Partial rRNAs, or RNAPs on rrn genes
		rna_counts = TableReader(os.path.join(simOutDir, "RNACounts"))
		partial_rRNA_cistron_counts = rna_counts.readColumn("partial_rRNA_cistron_counts")
		rna_id_to_index = {rna_id: i for (i, rna_id) in enumerate(
			rna_counts.readAttribute("rRNA_cistron_ids"))}
		rRNA_idxs = np.array([rna_id_to_index[rRNA_id] for rRNA_id in cistron_rRNA_ids])
		RNAP_on_rRNA_counts = partial_rRNA_cistron_counts[:, rRNA_idxs]
		grouped_RNAP_on_rRNA_counts = np.array(
			rRNA_cistron_to_type_mapping_matrix @ RNAP_on_rRNA_counts.T)

		# Active RNAP fraction on rrn genes
		unique_molecules = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		unique_molecule_ids = unique_molecules.readAttribute("uniqueMoleculeIds")
		rnap_idx = unique_molecule_ids.index('active_RNAP')
		unique_molecule_counts = unique_molecules.readColumn("uniqueMoleculeCounts")
		active_RNAP_counts = unique_molecule_counts[:, rnap_idx]
		RNAP_active_fraction_on_rRNA = np.array(
			[RNAP_counts / active_RNAP_counts for RNAP_counts in RNAP_on_rRNA_counts.T]).T
		grouped_RNAP_active_fraction_on_rRNA = np.array(
			rRNA_cistron_to_type_mapping_matrix @ RNAP_active_fraction_on_rRNA.T)

		# Counts of accumulated rRNAs and ppGpp
		ribosome_30S_id = [sim_data.molecule_ids.s30_full_complex]
		ribosome_50S_id = [sim_data.molecule_ids.s50_full_complex]
		rRNA_5S_ids = sim_data.molecule_groups.s50_5s_rRNA
		rRNA_16S_ids = sim_data.molecule_groups.s30_16s_rRNA
		rRNA_23S_ids = sim_data.molecule_groups.s50_23s_rRNA
		ppGpp_id = [sim_data.molecule_ids.ppGpp]

		(ribosome_30S_counts, ribosome_50S_counts, rRNA_5S_counts,
		rRNA_16S_counts, rRNA_23S_counts, ppGpp_counts) = read_bulk_molecule_counts(
			simOutDir, (ribosome_30S_id, ribosome_50S_id, rRNA_5S_ids,
			rRNA_16S_ids, rRNA_23S_ids, ppGpp_id))
		free_rRNA_5S_counts = np.sum(rRNA_5S_counts, axis=1)
		free_rRNA_16S_counts = np.sum(rRNA_16S_counts, axis=1)
		free_rRNA_23S_counts = np.sum(rRNA_23S_counts, axis=1)

		ribosome_idx = unique_molecule_ids.index('active_ribosome')
		active_ribosome = unique_molecule_counts[:, ribosome_idx]
		total_30S_counts = ribosome_30S_counts + active_ribosome
		total_50S_counts = ribosome_50S_counts + active_ribosome
		total_rRNA_5S_counts = free_rRNA_5S_counts + total_50S_counts
		total_rRNA_16S_counts = free_rRNA_16S_counts + total_30S_counts
		total_rRNA_23S_counts = free_rRNA_23S_counts + total_50S_counts

		# Concentration of ppGpp
		counts_to_molar = TableReader(os.path.join(simOutDir, 'EnzymeKinetics')).readColumn('countsToMolar').squeeze()
		ppGpp_conc = ppGpp_counts * counts_to_molar * 1000  # uM

		## Make figure
		# ncols is number of rRNA operons + 1 for total
		ncols = len(rRNA_operon_names) + 1
		nrows = 8

		fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3.5 * ncols, 3.5 * nrows), sharex='all')

		# Make the plots
		axs[0, 0].set_ylabel("Actual initiations per second", labelpad=LABELPAD)
		axs[1, 0].set_ylabel("Expected initiations per second", labelpad=LABELPAD)
		axs[2, 0].set_ylabel("Synthesis probabilities", labelpad=LABELPAD)
		axs[3, 0].set_ylabel("Gene dosages", labelpad=LABELPAD),
		axs[4, 0].set_ylabel("Initiations per second per gene copy", labelpad=LABELPAD)
		axs[5, 0].set_ylabel("Active RNAP counts", labelpad=LABELPAD)
		axs[6, 0].set_ylabel("Fraction of active RNAPs", labelpad=LABELPAD)
		axs[7, 0].set_ylabel("Counts")

		for i in range(ncols - 1):
			this_operon_idxs = np.where(rRNA_TU_to_cistron_mapping_matrix[i, :] == 1)[0]

			for (n, idx) in enumerate(this_operon_idxs):
				axs[0, i].plot(
					time, rRNA_actual_rna_init_rate[:, idx],
					c=COLORS[n], label=cistron_rRNA_ids[idx])
				axs[0, i].legend(loc=1)
				axs[0, i].set_ylim(0, 10)

				axs[1, i].plot(
					time, rRNA_expected_rna_init_rate[:, idx],
					c=COLORS[n], label=cistron_rRNA_ids[idx])
				axs[1, i].legend(loc=1)
				axs[1, i].set_ylim(0, 3)

				axs[2, i].plot(time, rRNA_actual_synth_prob[:, idx],
							   c=COLORS[n],label=cistron_rRNA_ids[idx]+" actual")
				axs[2, i].plot(time, rRNA_target_synth_prob[:, idx],
							   linestyle='dashed', c=COLORS[n],
							   label=cistron_rRNA_ids[idx]+" target")
				axs[2, i].legend(loc=1)
				axs[2, i].set_ylim(0, 0.03)

				axs[3, i].plot(time, rRNA_gene_dosages[:, idx],
							   c=COLORS[n], label=cistron_rRNA_ids[idx])
				axs[3, i].legend(loc=1)
				axs[3, i].set_ylim(0, 10)

				axs[4, i].plot(time, rRNA_initiation_per_gene_dosage[:, idx],
							   c=COLORS[n], label=cistron_rRNA_ids[idx])
				axs[4, i].legend(loc=1)
				axs[4, i].set_ylim(0, 1)

				axs[5, i].plot(time, RNAP_on_rRNA_counts[:, idx],
							   c=COLORS[n], label=cistron_rRNA_ids[idx])
				axs[5, i].legend(loc=1)
				axs[5, i].set_ylim(0, 200)

				axs[6, i].plot(time, RNAP_active_fraction_on_rRNA[:, idx],
							   c=COLORS[n], label=cistron_rRNA_ids[idx])
				axs[6, i].legend(loc=1)
				axs[6, i].set_ylim(0, 0.15)

			axs[0, i].set_title(rRNA_operon_names[i])
			axs[nrows - 1, i].set_xlabel("Time (min)", labelpad=LABELPAD, fontsize=FONTSIZE)

		for idx in range(np.shape(rRNA_cistron_to_type_mapping_matrix)[0]):
			axs[0, ncols - 1].plot(
				time, grouped_rRNA_actual_rna_init_rate[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx])
			axs[0, ncols - 1].legend(loc=1)

			axs[1, ncols - 1].plot(
				time, grouped_rRNA_expected_rna_init_rate[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx])
			axs[1, ncols - 1].legend(loc=1)

			axs[2, ncols - 1].plot(
				time, grouped_rRNA_actual_synth_prob[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx] + " actual")
			axs[2, ncols - 1].plot(
				time, grouped_rRNA_target_synth_prob[idx],
				c=COLORS[idx], linestyle='dashed',
				label=rRNA_cistron_types[idx] + " target")
			axs[2, ncols - 1].legend(loc=1)

			axs[3, ncols - 1].plot(
				time, grouped_rRNA_gene_dosages[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx])
			axs[3, ncols - 1].legend(loc=1)

			axs[4, ncols - 1].plot(
				time, grouped_rRNA_initiation_per_gene_dosage[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx])
			axs[4, ncols - 1].legend(loc=1)

			axs[5, ncols - 1].plot(
				time, grouped_RNAP_on_rRNA_counts[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx])
			axs[5, ncols - 1].legend(loc=1)

			axs[6, ncols - 1].plot(
				time, grouped_RNAP_active_fraction_on_rRNA[idx],
				c=COLORS[idx], label=rRNA_cistron_types[idx])
			axs[6, ncols - 1].legend(loc=1)

		axs[0, ncols - 1].set_title("All operons")
		axs[nrows - 1, ncols - 1].set_xlabel(
			"Time (min)", labelpad=LABELPAD, fontsize=FONTSIZE)

		axs[7, 0].plot(time, total_rRNA_5S_counts, c=COLORS[0], label="5S rRNA")
		axs[7, 0].plot(time, total_rRNA_23S_counts, c=COLORS[1], label="23S rRNA")
		axs[7, 0].plot(time, total_rRNA_16S_counts, c=COLORS[2], label="16S rRNA")
		axs[7, 0].set_title("Total rRNA counts")
		axs[7, 0].legend(loc=1)

		axs[7, 1].plot(time, total_30S_counts, label="30S subunit")
		axs[7, 1].plot(time, total_50S_counts, label="50S subunit")
		axs[7, 1].set_title("Total ribosomal subunit counts")
		axs[7, 1].legend(loc=1)

		axs[7, 2].plot(time, free_rRNA_5S_counts, c=COLORS[0], label="5S rRNA")
		axs[7, 2].plot(time, free_rRNA_23S_counts, c=COLORS[1], label="23S rRNA")
		axs[7, 2].plot(time, free_rRNA_16S_counts, c=COLORS[2], label="16S rRNA")
		axs[7, 2].set_title("Free rRNA counts")
		axs[7, 2].legend(loc=1)

		axs[7, 3].plot(time, ppGpp_conc)
		axs[7, 3].set_title("ppGpp concentration")
		axs[7, 3].set_ylabel("Concentration (uM)")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

if __name__ == "__main__":
	Plot().cli()
