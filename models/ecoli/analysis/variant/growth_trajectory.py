"""
Template for variant analysis plots
"""

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


def plot(ax, x, y, sim_time=None, xlabel=None, ylabel=None, label=None, background=False, markersize=3):
	if background:
		kwargs = {'alpha': 0.2, 'linewidth': 0.6, 'color': 'black'}
	else:
		kwargs = {'alpha': 0.7, 'label': label}
		ax.plot(x[0], y[0], 'og', markersize=markersize)
		ax.plot(x[-1], y[-1], 'or', markersize=markersize)

	trace, = ax.plot(x, y, **kwargs)
	if sim_time is not None:
		time_hr = np.floor(sim_time / 3600)
		hour_markers = np.where(np.diff(time_hr))[0] + 1
		ax.plot(x[hour_markers], y[hour_markers], 'o', alpha=kwargs['alpha'], markersize=markersize, color=trace.get_color())

	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel(ylabel, fontsize=8)
	ax.tick_params(labelsize=6)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		# Create plot
		_, axes = plt.subplots(3, 2, figsize=(8, 12))

		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		for variant in variants:
			all_mass_means = []
			all_growth_means = []
			all_ratio_means = []
			all_ratio_ma = []
			all_growth_ma = []
			all_protein_growth_ma = []
			all_rna_growth_ma = []
			all_small_mol_growth_ma = []
			all_times = []
			for seed in ap.get_seeds(variant):
				cell_paths = ap.get_cells(variant=[variant], seed=[seed])

				# Load data
				sim_time = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True).squeeze()
				time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec', remove_first=True).squeeze()
				growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True).squeeze()
				protein = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
				rna = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
				protein_growth = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', fun=growth_function).squeeze() / time_step
				rna_growth = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', fun=growth_function).squeeze() / time_step
				small_mol_growth = read_stacked_columns(cell_paths, 'Mass', 'smallMoleculeMass', fun=growth_function).squeeze() / time_step
				growth_means = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True, fun=np.mean).squeeze()
				protein_means = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True, fun=np.mean).squeeze()
				rna_means = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True, fun=np.mean).squeeze()
				mass_means = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True, fun=np.mean).squeeze()

				if len(np.unique(time_step)) > 1:
					raise ValueError('Check plot implementation to handle variable time step across sims.')

				moving_window = min(201, len(growth))
				convolution_array = np.ones(moving_window) / moving_window

				# Process data
				ratio = rna / protein
				ratio_means = rna_means / protein_means
				growth_ma = np.convolve(growth, convolution_array, mode='valid')
				ratio_ma = np.convolve(ratio, convolution_array, mode='valid')
				protein_growth_ma = np.convolve(protein_growth, convolution_array, mode='valid')
				rna_growth_ma = np.convolve(rna_growth, convolution_array, mode='valid')
				small_mol_growth_ma = np.convolve(small_mol_growth, convolution_array, mode='valid')
				time_ma = np.convolve(sim_time, convolution_array, mode='valid')

				plot(axes[0, 0], mass_means, growth_means, background=True)
				plot(axes[1, 0], ratio_means, growth_means, background=True)
				plot(axes[2, 0], ratio_ma, growth_ma, background=True)
				plot(axes[0, 1], ratio_ma, protein_growth_ma, background=True)
				plot(axes[1, 1], ratio_ma, rna_growth_ma, background=True)
				plot(axes[2, 1], ratio_ma, small_mol_growth_ma, background=True)

				all_mass_means.append(mass_means)
				all_growth_means.append(growth_means)
				all_ratio_means.append(ratio_means)
				all_ratio_ma.append(ratio_ma)
				all_growth_ma.append(growth_ma)
				all_protein_growth_ma.append(protein_growth_ma)
				all_rna_growth_ma.append(rna_growth_ma)
				all_small_mol_growth_ma.append(small_mol_growth_ma)
				all_times.append(time_ma)

			min_length = min([len(data) for data in all_growth_means])
			stacked_mass_means = np.vstack([data[:min_length] for data in all_mass_means]).mean(0)
			stacked_growth_means = np.vstack([data[:min_length] for data in all_growth_means]).mean(0)
			stacked_ratio_means = np.vstack([data[:min_length] for data in all_ratio_means]).mean(0)

			min_length_ma = min([len(data) for data in all_growth_ma])
			stacked_ratio_ma = np.vstack([data[:min_length_ma] for data in all_ratio_ma]).mean(0)
			stacked_growth_ma = np.vstack([data[:min_length_ma] for data in all_growth_ma]).mean(0)
			stacked_protein_growth_ma = np.vstack([data[:min_length_ma] for data in all_protein_growth_ma]).mean(0)
			stacked_rna_growth_ma = np.vstack([data[:min_length_ma] for data in all_rna_growth_ma]).mean(0)
			stacked_small_mol_growth_ma = np.vstack([data[:min_length_ma] for data in all_small_mol_growth_ma]).mean(0)
			stacked_times = np.vstack([data[:min_length_ma] for data in all_times]).mean(0)

			plot(axes[0, 0], stacked_mass_means, stacked_growth_means,
				xlabel='Cell cycle mass', ylabel='Cell cycle growth', label=variant)
			plot(axes[1, 0], stacked_ratio_means, stacked_growth_means,
				xlabel='Cell cycle RNA/protein', ylabel='Cell cycle growth', label=variant)
			plot(axes[2, 0], stacked_ratio_ma, stacked_growth_ma, sim_time=stacked_times,
				xlabel='RNA/protein', ylabel='Growth', label=variant)
			plot(axes[0, 1], stacked_ratio_ma, stacked_protein_growth_ma, sim_time=stacked_times,
				xlabel='RNA/protein', ylabel='Protein growth', label=variant)
			plot(axes[1, 1], stacked_ratio_ma, stacked_rna_growth_ma, sim_time=stacked_times,
				xlabel='RNA/protein', ylabel='RNA growth', label=variant)
			plot(axes[2, 1], stacked_ratio_ma, stacked_small_mol_growth_ma, sim_time=stacked_times,
				xlabel='RNA/protein', ylabel='Small molecule growth', label=variant)

		for ax in axes.reshape(-1):
			self.remove_border(ax)
			if len(variants) > 1:
				ax.legend(fontsize=6)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
