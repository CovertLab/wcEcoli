"""
Template for variant analysis plots

TODO:
	add labels for variants
	add labels for bar plot
	add x labels
	2D density plot instead of points
"""

import os
import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader


# +1 to make the window balanced on either side of the point of interest
MA_WINDOWS = np.array([10, 30, 60, 100, 300, 600, 1000, 3000]) + 1


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_bar(self, ax, cell_property, growth_rates, xlabel, tick_labels, growth_labels=None):
		n_subproperties = cell_property.shape[0]
		n_growth_rates = len(growth_rates)
		all_r = np.zeros((n_subproperties, n_growth_rates))
		for i, property in enumerate(cell_property):
			for j, growth in enumerate(growth_rates):
				all_r[i, j] = stats.pearsonr(property[:len(growth)], growth)[0]

		width = 0.8 / n_growth_rates
		offsets = np.arange(n_growth_rates) * width - 0.4 + width / 2
		x = np.arange(n_subproperties)
		for i, (r, offset) in enumerate(zip(all_r.T, offsets)):
			label = None if growth_labels is None else growth_labels[i]
			ax.bar(x + offset, r, width, label=label)

		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xticks(x)
		ax.set_xticklabels(tick_labels, fontsize=6, rotation=45, ha='right')
		ax.set_xlabel(xlabel, fontsize=8)

		return ax.get_legend_handles_labels()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids]

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		# TODO: add plots for ppGpp, ribosomes, RNAP etc
		scaling = 3
		n_cols = 9
		plt.figure(figsize=(scaling*n_cols, scaling*n_variants))
		gs = gridspec.GridSpec(nrows=n_variants, ncols=n_cols)

		for row, variant in enumerate(variants):
			cell_paths = ap.get_cells(variant=[variant])

			# Load attributes
			unique_molecule_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
			unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
			rnap_idx = unique_molecule_ids.index('active_RNAP')
			ribosome_idx = unique_molecule_ids.index('active_ribosome')

			# Load data
			time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec',
				remove_first=True).squeeze()
			sim_time = read_stacked_columns(cell_paths, 'Main', 'time',
				remove_first=True).squeeze()
			counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
				remove_first=True).squeeze()
			growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
				remove_first=True).squeeze() * 3600
			aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc',
				remove_first=True).T
			uncharged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'uncharged_trna_conc',
				remove_first=True).T
			charged_trna = read_stacked_columns(cell_paths, 'GrowthLimits', 'charged_trna_conc',
				remove_first=True).T
			synthetase_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'synthetase_conc',
				remove_first=True).T
			aa_supply_enzymes_fwd = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_fwd',
				remove_first=True).T * counts_to_molar
			aa_supply_enzymes_rev = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_supply_enzymes_rev',
				remove_first=True).T * counts_to_molar
			ppgpp_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'ppgpp_conc',
				remove_first=True).squeeze()
			unique_mol_counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
				remove_first=True)

			# Derived values
			fraction_charged = charged_trna / (uncharged_trna + charged_trna)
			active_ribosome_conc = unique_mol_counts[:, ribosome_idx] * counts_to_molar
			active_rnap_conc = unique_mol_counts[:, rnap_idx] * counts_to_molar

			# Apply moving average to growth to calculate correlation between current value and future growth
			all_growth = [growth] + [
				np.convolve(growth, np.ones(window) / window, mode='valid')
				for window in MA_WINDOWS
				]
			ma_labels = ['Current growth'] + [
				f'{np.convolve(time_step, np.ones(window-1), mode="valid").mean():.0f} sec window'
				for window in MA_WINDOWS
				]

			# Create subplots for this variant
			combined_data = np.vstack((ppgpp_conc, active_ribosome_conc, active_rnap_conc, sim_time))
			combined_labels = ['ppGpp', 'Ribosome', 'RNAP', 'Time']
			ax = plt.subplot(gs[row, 0])
			handles, labels = self.plot_bar(ax, combined_data, all_growth,
				'Other conc', combined_labels, growth_labels=ma_labels)
			ax.set_ylabel('Growth correlation', fontsize=8)

			ax = plt.subplot(gs[row, -1])
			ax.axis('off')
			ax.legend(handles, labels, loc='center', frameon=False, fontsize=8)

			plot_data = [
				(aa_conc, 'Amino acid conc'),
				(uncharged_trna, 'Uncharged tRNA conc'),
				(charged_trna, 'Charged tRNA conc'),
				(fraction_charged, 'Fraction charged'),
				(synthetase_conc, 'Synthetase conc'),
				(aa_supply_enzymes_fwd, 'Forward enzyme conc'),
				(aa_supply_enzymes_rev, 'Reverse enzyme conc'),
				]
			for col, (property, label) in enumerate(plot_data):
				self.plot_bar(plt.subplot(gs[row, col + 1]), property, all_growth, label, aa_ids)


		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
