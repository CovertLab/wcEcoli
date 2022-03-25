"""
Template for variant analysis plots

TODO:
 - actual time
 - all aa for each variant
 - divide by control
 - highlight mutants
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate, signal

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants import remove_aa_inhibition
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


def butter_highpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	normal_cutoff = cutoff / nyq
	b, a = signal.butter(order, normal_cutoff, btype = "high", analog = False)
	return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
	b, a = butter_highpass(cutoff, fs, order=order)
	y = signal.filtfilt(b, a, data)
	return y


def switch(data, input_drop=5, input_ma=60, deriv_ma=300, splrep=False):
	if input_drop > 1:
		n_points = len(data) // input_drop
		n_drop = len(data) % input_drop
		if n_drop == 0:
			index = slice(None)
		else:
			index = slice(-n_drop)
		data = data[index].reshape(n_points, -1).mean(1)

	if input_ma:
		input_ma = min(len(data) - 1, input_ma)
		data = np.convolve(data, np.ones(input_ma) / input_ma, 'valid')

	t = np.arange(len(data)) / 3600 * input_drop

	if splrep:
		spline = interpolate.splrep(t, data, k=5)
		deriv = interpolate.splev(t, spline, der=2)
	else:
		spline = interpolate.CubicSpline(t, data)
		deriv = spline.derivative(2)(t)

	if deriv_ma:
		deriv_ma = min(len(deriv) - 1, deriv_ma)
		deriv = np.convolve(deriv, np.ones(deriv_ma) / deriv_ma, 'same')

	sign = np.sign(deriv)

	return np.sum(sign[:-1] != sign[1:])

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		variants = self.ap.get_variants()
		n_variants = len(variants)

		baseline = None
		all_corr = {}
		all_periods = {}
		all_periods_deriv = {}
		all_periods_mean_adjusted = {}
		all_conc_mean = {}
		all_conc_std = {}
		autocorrelate = lambda x: signal.correlate(x, x, method='fft')
		for variant in variants:
			var_corr = []
			Pxx_all = None
			# TODO: plot each variant? or make cohort?
			var_conc = []
			periods = []
			periods_deriv = []
			for seed in self.ap.get_seeds(variant):
				cell_paths = self.ap.get_cells(variant=[variant], seed=[seed])

				if not np.all([self.ap.get_successful(cell) for cell in cell_paths]):
					continue

				# Load data
				sim_time = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True).squeeze()
				aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc', remove_first=True).T

				corr = np.apply_along_axis(autocorrelate, 1, aa_conc)
				var_corr.append(corr[:, corr.shape[1]//2:] / corr.mean(1).reshape(-1, 1))

				# Other approaches (combined it should be pearson correlation)
				# corr = np.apply_along_axis(autocorrelate, 1, aa_conc - aa_conc.mean(1).reshape(-1, 1))
				# var_corr.append(corr[:, corr.shape[1]//2:] / corr[:, corr.shape[1]//2].reshape(-1, 1))

				# # Crude approach to oscillations for leu
				# print(f'\n{variant}')
				# mean = aa_conc.mean(1)[10]
				# std = aa_conc.std(1)[10]
				# print(std / mean)
				# above = aa_conc[10, :] > mean
				# switches = np.sum(above[:-1] != above[1:])
				# total_time = (sim_time[-1] - sim_time[0]) / 3600
				# period = total_time / ((switches - 1) / 2)
				# print(switches)
				# print(period)


				total_time = (sim_time[-1] - sim_time[0]) / 3600
				above = aa_conc > aa_conc.mean(1).reshape(-1, 1)
				switches = np.sum(above[:, :-1] != above[:, 1:], axis=1)
				period = total_time / (switches / 2)
				period[period == np.inf] = np.nan
				var_conc.append(aa_conc)
				periods.append(period)

				switches = np.apply_along_axis(switch, 1, aa_conc)
				period_deriv = total_time / (switches / 2)
				period_deriv[period_deriv == np.inf] = np.nan
				periods_deriv.append(period_deriv)

				#
				# np.save('leu.npy', aa_conc[10, :])
				#
				# # Normalize amplitude
				# mean = aa_conc.mean(1).reshape(-1, 1)
				# std = aa_conc.std(1).reshape(-1, 1)
				# aa_conc = (aa_conc - mean) / std
				#
				# # Low pass filter
				# moving_window = min(1001, aa_conc.shape[1])
				# convolution_array = np.ones(moving_window) / moving_window
				# # aa_conc = np.apply_along_axis(np.convolve, 1, aa_conc, convolution_array, mode='valid')
				#
				# n_samples = aa_conc.shape[1]
				# print(n_samples)
				# freq = (sim_time[-1] - sim_time[0]) / n_samples
				# n_points = min(n_samples, 7200)
				#
				# # High pass filter
				# # aa_conc = butter_highpass_filter(aa_conc, 0.0005, freq)
				#
				# # Periodogram
				# f, Pxx = signal.periodogram(aa_conc, fs=1/freq, nfft=n_points, window='flattop', scaling='spectrum')
				#
				# # fft
				# ## Taper window to reduce leakage: http://qingkaikong.blogspot.com/2016/10/signal-processing-why-do-we-need-taper.html
				# # aa_conc *= signal.cosine(n_samples)
				# # f = np.fft.fftfreq(n_points) / freq
				# # Pxx = np.abs(np.real(np.fft.fft(aa_conc, n_points)))
				# # mask = f > 0
				# # f = f[mask]
				# # Pxx = Pxx[:, mask] * f
				#
				#
				# if Pxx_all is None:
				# 	Pxx_all = Pxx
				# else:
				# 	Pxx_all += Pxx
				#
				# period = 1 / f / 3600
				#
				# if variant == 0:
				# 	baseline = Pxx_all
				#
				# # TODO: divide mutant variant Pxx by wt to normalize? but must make sure the same scale
				# # TODO: set the length of analysis to minimal length of seed and add from all seeds
				# # break

			# TODO: period based on average time series across seeds

			all_corr[variant] = var_corr

			stacked_conc = np.hstack(var_conc)
			all_conc_mean[variant] = stacked_conc.mean(1)
			all_conc_std[variant] = stacked_conc.std(1)
			all_periods[variant] = np.nanmean(np.vstack(periods), 0)
			all_periods_deriv[variant] = np.nanmean(np.vstack(periods_deriv), 0)

			mean_periods = []
			for aa_conc in var_conc:
				total_time = aa_conc.shape[1] / 3600  # TODO: actual time?
				above = aa_conc > all_conc_mean[variant].reshape(-1, 1)
				switches = np.sum(above[:, :-1] != above[:, 1:], axis=1)
				period = total_time / (switches / 2)
				period[period == np.inf] = np.nan
				mean_periods.append(period)
			all_periods_mean_adjusted[variant] = np.nanmean(np.vstack(mean_periods), 0)

		# Save average data for comparison across runs
		with open(f'{os.path.join(plotOutDir, plotOutFileName)}.tsv', 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			headers = ['Variant']
			for aa in aa_ids:
				headers += [f'{aa} mean', f'{aa} std', f'{aa} CV', f'{aa} period']
			writer.writerow(headers)

			for variant in variants:
				cols = [variant]
				for mean, std, period in zip(all_conc_mean[variant], all_conc_std[variant], all_periods[variant]):
					cols += [mean, std, std / mean, period]

				writer.writerow(cols)

		control_idx = 0
		all_conc_cv = {
			variant: all_conc_std[variant] / all_conc_mean[variant]
			for variant in variants
			}
		width = 0.4

		# Plot for each amino acid for all variants
		def plot_bar(data, ylabel, file_label=''):
			plt.figure(figsize=(1.5, 1.5))

			controls = [
				data.get(control_idx, np.zeros(n_aas))[aa_ids.index(remove_aa_inhibition.get_aa_and_ki_factor(variant)[0])]
				for variant in variants
				if variant != control_idx
				]
			mutants = [
				data[variant][aa_ids.index(remove_aa_inhibition.get_aa_and_ki_factor(variant)[0])]
				for variant in variants
				if variant != control_idx
				]
			xlabels = [
				remove_aa_inhibition.get_aa_and_ki_factor(variant)[0][:-3]
				for variant in variants
				if variant != control_idx
				]
			x = range(len(mutants))
			# TODO: plot value for all mutants not just specific to aa
			plt.bar(x, controls, width=-width, align='edge')
			plt.bar(x, mutants, width=width, align='edge')

			plt.xticks(x, xlabels, fontsize=8, rotation=45)
			plt.ylabel(ylabel, fontsize=8)

			plt.gca().tick_params(labelsize=8)

			self.remove_border()
			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + file_label, metadata)
			plt.close('all')

		plot_bar(all_periods, 'Period')
		plot_bar(all_conc_cv, 'Coefficient of variation', file_label='_cv')
		plot_bar(all_periods_deriv, 'Period', file_label='_deriv')
		plot_bar(all_periods_mean_adjusted, 'Period', file_label='_mean')

		mean_corr = {}
		for variant, var_corr in all_corr.items():
			max_n = max([corr.shape[1] for corr in var_corr])
			n = np.zeros(max_n)
			total_corr = np.zeros((n_aas, max_n))
			for corr in var_corr:
				n_timepoints = corr.shape[1]
				n[:n_timepoints] += 1
				total_corr[:, :n_timepoints] += corr
			mean_corr[variant] = total_corr / n

		control_corr = mean_corr.get(0)

		def plot(label='', normalized=False):
			_, axes = plt.subplots(nrows=n_variants, ncols=n_aas, figsize=(30, 2*n_variants))

			for row, all_corr, in enumerate(mean_corr.values()):
				for col, corr in enumerate(all_corr):
					ax = axes[row, col]

					if normalized and control_corr is not None:
						control_pad = 0.95  # approaches 0 at the longest lag which can cause ratio to explode
						n_overlap = int(min(len(corr), control_corr.shape[1]*control_pad))
						ax.plot(np.arange(n_overlap) / 3600, corr[:n_overlap] / control_corr[col, :n_overlap])
						ax.axhline(1, color='k', alpha=0.5, linestyle='--', linewidth=1)
					else:
						ax.plot(np.arange(len(corr)) / 3600, corr)
						if control_corr is not None:
							ax.plot(np.arange(control_corr.shape[1]) / 3600, control_corr[col, :])

					ax.tick_params(labelsize=6)

					if row == 0:
						ax.set_title(aa_ids[col], fontsize=6)

					if col == 0:
						ax.set_ylabel(f'Variant {variants[row]}\nautocorrelation', fontsize=6)

					if row == n_variants - 1:
						ax.set_xlabel('Lag (hr)', fontsize=6)

			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + label, metadata)
			plt.close('all')

		plot(label='_lag')
		plot(label='_normalized', normalized=True)

		# # Not including baseline can provide some insights but need to figure out how to normalize
		# # over the range of periods with a gradual upslope
		# # Using baseline gives a depression of short periods in cases where inhibition is removed
		# y = Pxx_all if baseline is None else Pxx_all / baseline
		# for ax, P in zip(axes.flatten(), y):
		# 	# ax.plot(period[:-1], P[:-1] / np.sqrt(-np.diff(period*3600)))
		# 	# ax.plot(period, np.sqrt(P))
		# 	# ax.plot(period, P / f)
		# 	# ax.plot(f, P / f)
		# 	# ax.plot(f, P)
		# 	ax.plot(period, P)
		# 	# ax.plot(period, np.sqrt(P) * f)
		# 	ax.set_xscale('log')
		# 	ax.set_yscale('log')
		# 	ax.set_xticks(np.logspace(-3, 0, 4))
		# 	ax.axhline(1, color='k', linestyle='--', alpha=0.5)  # useful for comparison to baseline



if __name__ == "__main__":
	Plot().cli()
