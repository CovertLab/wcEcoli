'''
Plots fold change for molecules related to tRNA charging (charged/uncharged tRNA
and synthetases) along with ribosome elongation rate.

Useful for seeing if any molecules have deviated far from initial values and
are possible causes for changes to the ribosomal elongation rate if the rate is
specified by tRNA charging in the simulation. tRNA and synthetases are shown
on a per amino acid basis with a sum from all relevant species.

TODO: add amino acids and other metabolites involved in charging reactions
'''

import os
import pickle

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.analysis.plotting_tools import COLORS_SMALL
from wholecell.io.tablereader import TableReader
from wholecell.utils.sparkline import whitePadSparklineAxis


SECONDARY_COLOR = 'g'

def plot_ax(ax, x, y, secondary=False, dashed=False):
	'''
	Plots data and sets some common features for the axes

	Inputs:
		ax (matplotlib axes): axes to plot on
		x (numpy array of floats): x values to plot
		y (numpy array (1D or 2D) of floats): y values to plot
		secondary (bool): if True, plots with secondary color, otherwise
			defaults to a color cycle
		dashed (bool): if True, plots a dashed line instead of solid
	'''

	if secondary:
		ax.set_prop_cycle('color', SECONDARY_COLOR)
	else:
		ax.set_prop_cycle('color', COLORS_SMALL)  # Provides same color for each trace in different generations
	style = '--' if dashed else '-'
	ax.plot(x, y, style)

def post_plot_formatting(ax, division_times, y_label, draw_horizontal=None, y_lim=None, show_x_axis=False, secondary=False):
	'''
	Formats axes after all data has been plotted

	Inputs:
		ax (matplotlib axes): axes to plot on
		division_times (numpy array of floats): times of division
		y_label (str): label for the y axis
		draw_horizontal (float): if not None, draws a horizontal line at the given y position
		y_lim (numpy array of floats): if not None, specifies the lower and upper y limits
		show_x_axis (bool): if True, displays x axis
	'''

	if y_lim is not None:
		ax.set_ylim(y_lim)
	ax.set_xlim([0, division_times[-1]])
	whitePadSparklineAxis(ax, xAxis=show_x_axis, secondary=secondary)
	if show_x_axis:
		ax.set_xticks(division_times[:-1], minor=True)
		ax.set_xticks([0, division_times[-1]])

	if draw_horizontal is not None:
		ax.axhline(draw_horizontal, color='k', linestyle='--', linewidth=0.5, alpha=0.8)
		ax.set_yticks(np.hstack((ax.get_yticks(), draw_horizontal)))

	str_format = FormatStrFormatter('%.3g')
	ax.xaxis.set_major_formatter(str_format)
	ax.yaxis.set_major_formatter(str_format)

	if secondary:
		color = SECONDARY_COLOR
	else:
		color = 'k'
	ax.set_ylabel(y_label, fontsize=8, color=color)
	ax.tick_params(axis='both', labelsize=6)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		transcription = sim_data.process.transcription
		synthetase_names = transcription.synthetase_names
		uncharged_trna_names = transcription.uncharged_trna_names
		charged_trna_names = transcription.charged_trna_names
		aa_from_synthetase = transcription.aa_from_synthetase.T
		aa_from_trna = transcription.aa_from_trna.T

		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		mol_ids = sim_data.molecule_ids
		ppgpp_molecules = [mol_ids.RelA, mol_ids.SpoT, mol_ids.ppGpp]


		# Create plot and axes
		n_subplots = 13
		fig = plt.figure(figsize=(5, 20))
		growth_ax = plt.subplot(n_subplots, 1, 1)
		growth_ax2 = growth_ax.twinx()
		ppgpp_ax = plt.subplot(n_subplots, 1, 2)
		spot_ax = plt.subplot(n_subplots, 1, 3)
		rela_ax = spot_ax.twinx()
		synth_ax = plt.subplot(n_subplots, 1, 4)
		trna_ax = plt.subplot(n_subplots, 1, 5)
		expected_frac_ax = plt.subplot(n_subplots, 1, 6)
		frac_ax = plt.subplot(n_subplots, 1, 7)
		diff_ax = plt.subplot(n_subplots, 1, 8)
		uncharged_trna_ax = plt.subplot(n_subplots, 1, 9)
		charged_trna_ax = plt.subplot(n_subplots, 1, 10)
		ppgpp_syn_ax = plt.subplot(n_subplots, 1, 11)
		ppgpp_syn_ax2 = ppgpp_syn_ax.twinx()
		ppgpp_deg_ax = plt.subplot(n_subplots, 1, 12)
		ppgpp_deg_ax2 = ppgpp_deg_ax.twinx()
		legend_ax = plt.subplot(n_subplots, 1, 13)

		initial_synthetase_conc = None
		initial_uncharged_trna_conc = None
		initial_charged_trna_conc = None
		initial_ppgpp_protein_conc = None
		initial_total_trna_conc = None
		division_times = []
		total_elong = 0.
		total_growth = 0.
		total_ppgpp = 0.
		timesteps = 0.
		for sim_dir in self.ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))
			ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))
			mass_reader = TableReader(os.path.join(simOutDir, 'Mass'))
			enzyme_kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))
			growth_reader = TableReader(os.path.join(simOutDir, 'GrowthLimits'))

			# Load data
			time = main_reader.readColumn('time') / 3600
			division_times.append(time[-1])
			elong_rate = ribosome_reader.readColumn('effectiveElongationRate')
			growth_rate = mass_reader.readColumn('instantaneous_growth_rate') * 3600
			counts_to_molar = enzyme_kinetics_reader.readColumn('countsToMolar').reshape(-1, 1)
			(synthetase_counts, uncharged_trna_counts, charged_trna_counts, ppgpp_mol_counts
				) = read_bulk_molecule_counts(simOutDir,
				(synthetase_names, uncharged_trna_names, charged_trna_names, ppgpp_molecules))
			expected_fraction = growth_reader.readColumn('fraction_trna_charged') @ aa_from_trna / aa_from_trna.sum(0)
			rela_rxns = growth_reader.readColumn('rela_syn')
			spot_syn_rxns = growth_reader.readColumn('spot_syn')
			spot_deg_rxns = growth_reader.readColumn('spot_deg')
			spot_deg_inhibited = growth_reader.readColumn('spot_deg_inhibited')

			## Running totals for elongation and growth
			total_elong += elong_rate.sum()
			total_growth += np.nansum(growth_rate)

			## Synthetase counts
			synthetase_conc = counts_to_molar * np.dot(synthetase_counts, aa_from_synthetase)
			if initial_synthetase_conc is None:
				initial_synthetase_conc = synthetase_conc[1, :]
			normalized_synthetase_conc = synthetase_conc / initial_synthetase_conc

			## Uncharged tRNA counts
			uncharged_trna_conc = counts_to_molar * np.dot(uncharged_trna_counts, aa_from_trna)
			if initial_uncharged_trna_conc is None:
				initial_uncharged_trna_conc = uncharged_trna_conc[1, :]
			normalized_uncharged_trna_conc = uncharged_trna_conc / initial_uncharged_trna_conc

			## Charged tRNA counts
			charged_trna_conc = counts_to_molar * np.dot(charged_trna_counts, aa_from_trna)
			if initial_charged_trna_conc is None:
				initial_charged_trna_conc = charged_trna_conc[1, :]
			normalized_charged_trna_conc = charged_trna_conc / initial_charged_trna_conc

			## Total tRNA and fraction charged
			total_trna_conc = charged_trna_conc + uncharged_trna_conc
			if initial_total_trna_conc is None:
				initial_total_trna_conc = total_trna_conc[1, :]
			normalized_total_trna_conc = total_trna_conc / initial_total_trna_conc
			fraction_charged = charged_trna_conc / (charged_trna_conc + uncharged_trna_conc)

			## ppGpp related counts and concentration
			ppgpp_protein_conc = counts_to_molar * ppgpp_mol_counts[:, :2]
			if initial_ppgpp_protein_conc is None:
				initial_ppgpp_protein_conc = ppgpp_protein_conc[1, :]
			normalized_ppgpp_protein_conc = ppgpp_protein_conc / initial_ppgpp_protein_conc
			ppgpp_conc = ppgpp_mol_counts[:, 2] * counts_to_molar[:, 0] * 1000
			total_ppgpp += ppgpp_conc.sum()
			timesteps += len(ppgpp_conc)

			# Plot data
			plot_ax(growth_ax, time[1:], elong_rate[1:])  # [1:] to remove spike down
			plot_ax(growth_ax2, time[1:], growth_rate[1:], secondary=True)
			plot_ax(spot_ax, time, np.log2(normalized_ppgpp_protein_conc[:, 1]))
			plot_ax(rela_ax, time, np.log2(normalized_ppgpp_protein_conc[:, 0]), secondary=True)
			plot_ax(ppgpp_ax, time[1:], ppgpp_conc[1:])
			plot_ax(synth_ax, time, np.log2(normalized_synthetase_conc))
			plot_ax(trna_ax, time, np.log2(normalized_total_trna_conc))
			plot_ax(expected_frac_ax, time[1:], expected_fraction[1:, :])
			plot_ax(frac_ax, time, fraction_charged)
			plot_ax(diff_ax, time, fraction_charged - expected_fraction)
			plot_ax(uncharged_trna_ax, time, np.log2(normalized_uncharged_trna_conc))
			plot_ax(charged_trna_ax, time, np.log2(normalized_charged_trna_conc))
			plot_ax(ppgpp_syn_ax, time[1:], rela_rxns[1:, :])
			plot_ax(ppgpp_syn_ax2, time[1:], spot_syn_rxns[1:], secondary=True, dashed=True)
			plot_ax(ppgpp_deg_ax, time[1:], spot_deg_inhibited[1:, :])
			plot_ax(ppgpp_deg_ax2, time[1:], spot_deg_rxns[1:], secondary=True, dashed=True)

		elong_mean = total_elong / timesteps
		growth_mean = total_growth / timesteps
		ppgpp_mean = total_ppgpp / timesteps

		# Format plot axes
		post_plot_formatting(growth_ax, division_times, 'Ribosome\nElongation Rate', y_lim=[0, 22], draw_horizontal=elong_mean)
		post_plot_formatting(growth_ax2, division_times, 'Growth Rate\n(1/hr)', y_lim=0, draw_horizontal=growth_mean, secondary=True)
		post_plot_formatting(spot_ax, division_times, 'SpoT Conc\nFold Change', draw_horizontal=0)
		post_plot_formatting(rela_ax, division_times, 'RelA Conc\nFold Change', draw_horizontal=0, secondary=True)
		post_plot_formatting(ppgpp_ax, division_times, 'ppGpp Conc\n(uM)', y_lim=0, draw_horizontal=ppgpp_mean)
		post_plot_formatting(synth_ax, division_times, 'Synthetase Conc\nFold Change', draw_horizontal=0)
		post_plot_formatting(trna_ax, division_times, 'Total tRNA Conc\nFold Change', draw_horizontal=0)
		post_plot_formatting(expected_frac_ax, division_times, 'Expected fraction\ntRNA Charged', y_lim=[0, 1])
		post_plot_formatting(frac_ax, division_times, 'Actual fraction\ntRNA Charged', y_lim=[0, 1])
		post_plot_formatting(diff_ax, division_times, 'Difference in fraction\ntRNA Charged')
		post_plot_formatting(uncharged_trna_ax, division_times, 'Uncharged tRNA Conc\nFold Change', draw_horizontal=0)
		post_plot_formatting(charged_trna_ax, division_times, 'Charged tRNA Conc\nFold Change', draw_horizontal=0)
		post_plot_formatting(ppgpp_syn_ax, division_times, 'ppGpp Synthesis\nby RelA (uM/s)', y_lim=0)
		post_plot_formatting(ppgpp_syn_ax2, division_times, 'ppGpp Synthesis\nby SpoT (uM/s)', y_lim=0, secondary=True)
		post_plot_formatting(ppgpp_deg_ax, division_times, 'ppGpp degradation\ninhibited (uM/s)', y_lim=0, show_x_axis=True)
		post_plot_formatting(ppgpp_deg_ax2, division_times, 'ppGpp degradation\n(uM/s)', y_lim=0, show_x_axis=True, secondary=True)
		ppgpp_deg_ax.set_xlabel('Time (hr)', fontsize=8)

		# Format and display legend below all plots
		legend_ax.set_prop_cycle('color', COLORS_SMALL)
		legend_ax.axis('off')
		legend_ax.plot(0, np.zeros((1, n_aas)))
		legend_ax.legend(aa_ids, ncol=3, fontsize=6, loc=10)

		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
