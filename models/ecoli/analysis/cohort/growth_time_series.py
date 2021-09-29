"""
Show mean traces and confidence interval for growth related properties over
multiple initial seeds.
"""

import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def plot_time_series(self, gs, t_flat, y_flat, ylabel):
		ax = plt.subplot(gs)

		data = {}
		for t, y in zip(t_flat, y_flat):
			data.setdefault(t, []).append(y)

		mean = []
		std = []
		t = np.array(sorted(data.keys()))
		for _t in t:
			d = np.vstack(data[_t])
			mean.append(d.mean(0))
			std.append(d.std(0))
		mean = np.array(mean).squeeze()
		std = np.array(std).squeeze()

		# # Plot all single traces
		# new_cell = np.where(t_flat[:-1] > t_flat[1:])[0] + 1
		# splits = [0] + list(new_cell) + [None]
		# for start, end in zip(splits[:-1], splits[1:]):
		# 	ax.plot(t_flat[start:end], y_flat[start:end], 'k', alpha=0.05, linewidth=0.5)

		ax.plot(t, mean)
		if len(mean.shape) > 1:
			for m, s in zip(mean.T, std.T):
				ax.fill_between(t, m - s, m + s, alpha=0.1)
		else:
			ax.fill_between(t, mean - std, mean + std, alpha=0.1)

		ax.set_xlabel('Time (min)', fontsize=8)
		ax.set_ylabel(ylabel, fontsize=8)
		ax.tick_params(labelsize=6)
		self.remove_border(ax)

	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		transcription = sim_data.process.transcription
		ppgpp_id = sim_data.molecule_ids.ppGpp
		uncharged_trna_names = transcription.rna_data['id'][transcription.rna_data['is_tRNA']]
		charged_trna_names = transcription.charged_trna_names
		aa_from_trna = transcription.aa_from_trna.T

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		time = read_stacked_columns(cell_paths, 'Main', 'time',
			remove_first=True, ignore_exception=True).squeeze() / 60
		growth_rate = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
			remove_first=True, ignore_exception=True).squeeze() * 3600
		elong_rate = read_stacked_columns(cell_paths, 'RibosomeData', 'effectiveElongationRate',
			remove_first=True, ignore_exception=True).squeeze()
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True, ignore_exception=True).squeeze()
		(ppgpp_counts, uncharged_trna_counts, charged_trna_counts) = read_stacked_bulk_molecules(
			cell_paths, ([ppgpp_id], uncharged_trna_names, charged_trna_names),
			remove_first=True, ignore_exception=True)

		# Derived quantities
		ppgpp_conc = ppgpp_counts * counts_to_molar * 1000
		charged_trna_counts = charged_trna_counts @ aa_from_trna
		uncharged_trna_counts = uncharged_trna_counts @ aa_from_trna
		fraction_charged = charged_trna_counts / (uncharged_trna_counts + charged_trna_counts)

		plt.figure(figsize=(5, 15))
		gs = gridspec.GridSpec(4, 1)

		self.plot_time_series(gs[0, 0], time, growth_rate, 'Growth rate\n(1/hr)')
		self.plot_time_series(gs[1, 0], time, elong_rate, 'Elongation rate\n(AA/s)')
		self.plot_time_series(gs[2, 0], time, ppgpp_conc, 'ppGpp concentration\n(uM)')
		self.plot_time_series(gs[3, 0], time, fraction_charged, 'Fraction charged')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
