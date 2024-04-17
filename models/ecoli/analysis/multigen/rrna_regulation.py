"""
Plots timetraces of components in the model related to rRNA regulation.
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


COMPLEXATION_RXN_IDS = ['CPLX0-3953_RXN', 'CPLX0-3962_RXN']

class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Load from sim_data
		transcription = sim_data.process.transcription
		uncharged_trna_names = transcription.uncharged_trna_names
		charged_trna_names = transcription.charged_trna_names
		aa_from_trna = transcription.aa_from_trna.T
		aa_names = sim_data.molecule_groups.amino_acids
		trp_index = aa_names.index('TRP[c]')

		cell_paths = self.ap.get_cells()

		# Load data
		time = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True).flatten()
		gen_start_times = read_stacked_columns(
			cell_paths, 'Main', 'time', remove_first=True,
			fun=lambda x: x[0]).flatten()

		# Get rna_ids attribute from RnaSynthProb table in reference cell path
		reference_cell_path = cell_paths[0]
		sim_out_dir = os.path.join(reference_cell_path, 'simOut')
		rna_synth_prob_reader = TableReader(
			os.path.join(sim_out_dir, 'RnaSynthProb'))
		rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')

		# Get indexes of rRNAs in RnaSynthProb table
		rna_id_to_is_rRNA = {
			rna['id']: rna['is_rRNA'] for rna in transcription.rna_data
			}
		rrna_indexes_rna_synth_prob = np.array([
			i for (i, rna_id) in enumerate(rna_ids)
			if rna_id_to_is_rRNA[rna_id]
			])

		# Get reactionIDs attribute from ComplexationListener table
		complexation_reader = TableReader(
			os.path.join(sim_out_dir, 'ComplexationListener'))
		reaction_ids = complexation_reader.readAttribute('reactionIDs')

		# Get indexes of ribosomal subunit complexation reactions in
		# ComplexationListener table
		subunit_complexation_rxn_indexes = np.array([
			reaction_ids.index(rxn_id) for rxn_id in COMPLEXATION_RXN_IDS
			])

		# Get rna_ids attribute from RnapData table in reference cell path
		rna_synth_prob_reader = TableReader(
			os.path.join(sim_out_dir, 'RnapData'))
		rna_ids = rna_synth_prob_reader.readAttribute('rnaIds')

		# Get indexes of rRNAs in RnapData table
		rna_id_to_is_rRNA = {
			rna['id']: rna['is_rRNA'] for rna in transcription.rna_data
			}
		rrna_indexes_rnap_data = np.array([
			i for (i, rna_id) in enumerate(rna_ids)
			if rna_id_to_is_rRNA[rna_id]
			])

		# Get index of active ribosomes in the unique molecule counts reader
		unique_molecule_counts_reader = TableReader(
			os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_counts_reader.readAttribute(
			'uniqueMoleculeIds')
		active_ribosome_idx = unique_molecule_ids.index('active_ribosome')

		# Get summed rRNA synthesis probabilities
		rrna_synthesis_probs = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob',
			remove_first=True, fun=lambda x: x[:, rrna_indexes_rna_synth_prob])
		rrna_synthesis_probs = rrna_synthesis_probs.sum(axis=1)

		# Get summed rRNA copy numbers
		rrna_copy_numbers = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			remove_first=True, fun=lambda x: x[:, rrna_indexes_rna_synth_prob])
		rrna_copy_numbers = rrna_copy_numbers.sum(axis=1)

		# Get summed rRNA initiation event counts
		rrna_init_events = read_stacked_columns(
			cell_paths, 'RnapData', 'rnaInitEvent',
			remove_first=True, fun=lambda x: x[:, rrna_indexes_rnap_data])
		rrna_init_events = rrna_init_events.sum(axis=1)

		# Get number of ribosomal subunit complexation events
		ribosome_complexation_events = read_stacked_columns(
			cell_paths, 'ComplexationListener', 'complexationEvents',
			remove_first=True, fun=lambda x: x[:, subunit_complexation_rxn_indexes])
		ribosome_complexation_events = ribosome_complexation_events.sum(axis=1)

		# Get counts of active ribosomes
		unique_molecule_counts = read_stacked_columns(
			cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
			remove_first=True)
		active_ribosome_counts = unique_molecule_counts[:, active_ribosome_idx]

		# Get counts of imported amino acids
		imported_aas = read_stacked_columns(
			cell_paths, 'GrowthLimits', 'aa_import', remove_first=True)
		imported_aas = imported_aas.sum(axis=1)

		# Get counts of translated amino acids
		translated_aas = read_stacked_columns(
			cell_paths, 'RibosomeData', 'actualElongations',
			remove_first=True).squeeze()

		# Get fraction of charged leucine tRNAs
		(uncharged_trna_counts, charged_trna_counts, ) = read_stacked_bulk_molecules(
			cell_paths,	(uncharged_trna_names, charged_trna_names, ),
			remove_first=True)
		charged_trna_counts = charged_trna_counts @ aa_from_trna
		uncharged_trna_counts = uncharged_trna_counts @ aa_from_trna
		fraction_charged = charged_trna_counts / (
				uncharged_trna_counts + charged_trna_counts)
		trp_fraction_charged = fraction_charged[:, trp_index]

		# Get ppGpp concentrations
		ppgpp_conc = read_stacked_columns(
			cell_paths, 'GrowthLimits', 'ppgpp_conc',
			remove_first=True).squeeze()  # uM

		# Get instantaneous doubling times
		growth_rates = (1 / units.s) * read_stacked_columns(
			cell_paths, 'Mass', 'instantaneous_growth_rate',
			remove_first=True)
		doubling_times = (
				(1 / growth_rates) * np.log(2)).asNumber(units.min)

		def plot_ax(ax, y, ylabel, ylim, yticks, clip_on=False):
			ax.plot(time / 60, y, color='#8c8c8c', clip_on=clip_on)
			ax.set_ylabel(ylabel, fontsize=8)
			ax.spines["top"].set_visible(False)
			ax.spines["right"].set_visible(False)
			ax.spines["bottom"].set_position(("outward", 10))
			ax.spines["left"].set_position(("outward", 10))
			ax.spines["bottom"].set_visible(False)
			ax.get_xaxis().set_visible(False)
			ax.set_xlim([0, time[-1] / 60])
			ax.set_ylim(ylim)
			ax.set_yticks(yticks)

		plt.figure(figsize=(12, 10))

		# Plot per-copy transcription probability of rRNA operons
		ax1 = plt.subplot(11, 1, 1)
		plot_ax(
			ax1, rrna_synthesis_probs / rrna_copy_numbers,
			'p_trc per\nrRNA copy', [0, 0.015], [0, 0.015])

		# Plot copy numbers of rRNA operons
		ax2 = plt.subplot(11, 1, 2, sharex=ax1)
		plot_ax(
			ax2, rrna_copy_numbers, 'rRNA\ncopy #',
			[0, 60], [0, 60])

		# Plot total transcription probability of rRNA operons
		ax3 = plt.subplot(11, 1, 3)
		plot_ax(
			ax3, rrna_synthesis_probs,
			'total rRNA\np_trc', [0, 0.25], [0, 0.25])

		# Plot number of transcription initiation events
		ax4 = plt.subplot(11, 1, 4, sharex=ax1)
		plot_ax(
			ax4, rrna_init_events, 'rRNA\ninit #',
			[0, 80], [0, 40, 80])

		# Plot number of ribosome subunit complexation events
		ax5 = plt.subplot(11, 1, 5, sharex=ax1)
		plot_ax(
			ax5, ribosome_complexation_events, '# of comp.\nevents',
			[0, 150], [0, 150])

		# Plot counts of active ribosomes
		ax6 = plt.subplot(11, 1, 6, sharex=ax1)
		plot_ax(
			ax6, active_ribosome_counts, '# of active\nribosomes',
			[0, 80000], [0, 40000, 80000])

		# Plot number of imported AAs
		ax7 = plt.subplot(11, 1, 7, sharex=ax1)
		plot_ax(
			ax7, imported_aas, '# of aas\nimported',
			[0, 8e6], [0, 8e6])

		# Plot number of translated AAs
		ax8 = plt.subplot(11, 1, 8, sharex=ax1)
		plot_ax(
			ax8, translated_aas, '# of aas\ntranslated',
			[0, 2e6], [0, 1e6, 2e6])

		# Plot fraction of charged TRP tRNAs
		ax9 = plt.subplot(11, 1, 9, sharex=ax1)
		plot_ax(
			ax9, trp_fraction_charged, 'TRP charged\nfraction',
			[0, 1], [0, 0.5, 1])

		# Plot ppGpp concentrations
		ax10 = plt.subplot(11, 1, 10, sharex=ax1)
		plot_ax(
			ax10, ppgpp_conc, 'ppGpp\nconc (uM)',
			[0, 100], [0, 50, 100])

		# Plot doubling times
		ax11 = plt.subplot(11, 1, 11, sharex=ax1)
		plot_ax(
			ax11, doubling_times, 'DT\n(min)',
			[0, 80], [0, 40, 80], clip_on=True)
		ax11.get_xaxis().set_visible(True)
		ax11.spines["bottom"].set_visible(True)
		ax11.set_xlabel('Time (generations)')
		ax11.set_xticks((gen_start_times / 60).tolist() + [time[-1] / 60])
		ax11.set_xticklabels(np.arange(len(gen_start_times) + 1).tolist())

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')



if __name__ == '__main__':
	Plot().cli()
