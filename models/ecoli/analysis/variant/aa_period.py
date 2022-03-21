"""
Template for variant analysis plots
"""

from matplotlib import pyplot as plt
import numpy as np
from scipy import signal

from models.ecoli.analysis import variantAnalysisPlot
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


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		variants = self.ap.get_variants()

		baseline = None
		for variant in variants:
			Pxx_all = None
			# TODO: plot each variant? or make cohort?
			for seed in self.ap.get_seeds(variant):
				cell_paths = self.ap.get_cells(variant=[variant], seed=[seed], only_successful=True)

				# Load data
				sim_time = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True).squeeze()
				aa_conc = read_stacked_columns(cell_paths, 'GrowthLimits', 'aa_conc', remove_first=True).T

				# Crude approach to oscillations
				print(f'\n{variant}')
				mean = aa_conc.mean(1)[10]
				std = aa_conc.std(1)[10]
				print(std / mean)
				above = aa_conc[10, :] > mean
				switches = np.sum(above[:-1] != above[1:])
				print(switches)

				np.save('leu.npy', aa_conc[10, :])

				# Normalize amplitude
				mean = aa_conc.mean(1).reshape(-1, 1)
				std = aa_conc.std(1).reshape(-1, 1)
				aa_conc = (aa_conc - mean) / std

				# Low pass filter
				moving_window = min(1001, aa_conc.shape[1])
				convolution_array = np.ones(moving_window) / moving_window
				# aa_conc = np.apply_along_axis(np.convolve, 1, aa_conc, convolution_array, mode='valid')

				n_samples = aa_conc.shape[1]
				print(n_samples)
				freq = (sim_time[-1] - sim_time[0]) / n_samples
				n_points = min(n_samples, 7200)

				# High pass filter
				# aa_conc = butter_highpass_filter(aa_conc, 0.0005, freq)

				# Periodogram
				f, Pxx = signal.periodogram(aa_conc, fs=1/freq, nfft=n_points, window='flattop', scaling='spectrum')

				# fft
				## Taper window to reduce leakage: http://qingkaikong.blogspot.com/2016/10/signal-processing-why-do-we-need-taper.html
				# aa_conc *= signal.cosine(n_samples)
				# f = np.fft.fftfreq(n_points) / freq
				# Pxx = np.abs(np.real(np.fft.fft(aa_conc, n_points)))
				# mask = f > 0
				# f = f[mask]
				# Pxx = Pxx[:, mask] * f

				if Pxx_all is None:
					Pxx_all = Pxx
				else:
					Pxx_all += Pxx

				period = 1 / f / 3600

				if variant == 0:
					baseline = Pxx_all

				# TODO: divide mutant variant Pxx by wt to normalize? but must make sure the same scale
				# TODO: set the length of analysis to minimal length of seed and add from all seeds
				# break

		_, axes = plt.subplots(5, 5, figsize=(10, 10))

		# Not including baseline can provide some insights but need to figure out how to normalize
		# over the range of periods with a gradual upslope
		# Using baseline gives a depression of short periods in cases where inhibition is removed
		y = Pxx_all if baseline is None else Pxx_all / baseline
		for ax, P in zip(axes.flatten(), y):
			# ax.plot(period[:-1], P[:-1] / np.sqrt(-np.diff(period*3600)))
			# ax.plot(period, np.sqrt(P))
			# ax.plot(period, P / f)
			# ax.plot(f, P / f)
			# ax.plot(f, P)
			ax.plot(period, P)
			# ax.plot(period, np.sqrt(P) * f)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_xticks(np.logspace(-3, 0, 4))
			ax.axhline(1, color='k', linestyle='--', alpha=0.5)  # useful for comparison to baseline

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
