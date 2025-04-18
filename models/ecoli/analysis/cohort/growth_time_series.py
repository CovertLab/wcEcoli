"""
Show mean traces and confidence interval for growth related properties over
multiple initial seeds.
"""

import os
import pickle

from matplotlib import pyplot as plt
import matplotlib
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.single.ribosome_limitation import calculate_ribosome_excesses
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def plot_time_series(self, ax, t_flat, y_flat, ylabel, timeline, filtered_t,
			downsample=5, log_scale=False, single_cells=False, distinct=False,
			show_x_labels=True, plot_options=None, single_scale=1, shade=True, use_hr=False):
		if plot_options is None:
			plot_options = {}
		t_scaling = 60 if use_hr else 1

		# TODO: add trace labels as arg
		# Extract y data for each time point (assumes time step lines up across samples)
		data = {}
		for t, y in zip(t_flat, y_flat):
			if t not in filtered_t:
				data.setdefault(t, []).append(y)

		# Calculate mean and standard deviation for each group of time points
		# based on downsampling times
		t = []
		mean = []
		std = []
		all_times = np.array(sorted(data.keys()))
		n_dropped = len(all_times) % downsample
		drop_slice = -n_dropped if n_dropped else None
		times = all_times[:drop_slice].reshape(-1, downsample)
		if len(y_flat.shape) > 1:
			stack = np.vstack
		else:
			stack = np.hstack
		for ts in times:
			d = stack([data[_t] for _t in ts])
			t.append(ts.mean())
			mean.append(d.mean(0))
			std.append(d.std(0))
		t = np.array(t) / t_scaling
		mean = np.array(mean).squeeze()
		std = np.array(std).squeeze()

		# Plot all single traces
		if single_cells:
			new_cell = np.where(t_flat[:-1] > t_flat[1:])[0] + 1
			n_cells = len(new_cell) + 1
			splits = [0] + list(new_cell) + [None]

			# Scale transparency and width with number of cells to make more
			# readable with many traces
			linewidth = 1 / (1 + np.log(n_cells)) * single_scale
			if distinct:
				alphas = np.linspace(0.2, 0.8, n_cells)
				color = matplotlib.colors.TABLEAU_COLORS['tab:blue']
			else:
				alphas = 1 / (1 + n_cells) * np.ones(n_cells)
				color = 'k'

			# Plot each cell trace
			for start, end, alpha in zip(splits[:-1], splits[1:], alphas):
				# Select cell specific trace
				t_single = t_flat[start:end] / t_scaling
				y_single = y_flat[start:end]

				# Downsample data
				n_dropped = len(t_single) % downsample
				drop_slice = -n_dropped if n_dropped else None
				t_single = t_single[:drop_slice].reshape(-1, downsample).mean(axis=1)
				if len(dim := y_flat.shape) > 1:
					y_single = y_single[:drop_slice].reshape(-1, downsample, dim[1]).mean(axis=1)
				else:
					y_single = y_single[:drop_slice].reshape(-1, downsample).mean(axis=1)

				# Plot single trace
				ax.plot(t_single, y_single, color,
					alpha=alpha, linewidth=linewidth)

		# Plot mean as a trace and standard deviation as a shaded area
		ax.plot(t, mean, **plot_options)
		if shade:
			if len(mean.shape) > 1:
				for m, s in zip(mean.T, std.T):
					ax.fill_between(t, m - s, m + s, alpha=0.1)
			else:
				ax.axhline(mean.mean(), linestyle='--', color='k', linewidth=0.5)
				ax.fill_between(t, mean - std, mean + std, alpha=0.1)

		# Format axes
		if log_scale:
			ax.set_yscale('log')
		if show_x_labels:
			t_units = 'hr' if use_hr else 'min'
			ax.set_xlabel(f'Time ({t_units})', fontsize=8)
		else:
			ax.xaxis.set_ticklabels([])
		ax.set_ylabel(ylabel, fontsize=8)
		ax.tick_params(labelsize=8)
		self.remove_border(ax)

		# Show any media changes
		t_max = t.max()
		y_min, y_max = ax.get_ylim()
		for i, (t_media, media) in enumerate(timeline):
			t_media /= 60 * t_scaling
			if t_media < t_max:
				ax.axvline(t_media, color='k', linestyle='--', linewidth=0.5, alpha=0.3)
				if log_scale:
					y_pos = 10**(np.log10(y_max) - (np.log10(y_max) - np.log10(y_min)) * i * 0.03)
				else:
					y_pos = y_max - (y_max - y_min) * i * 0.03
				ax.text(t_media, y_pos, media, fontsize=6)

	def plot_hist(self, ax, data, min_val, max_val, label, n_bins=40, sf=2):
		def plot(d):
			patch = ax.hist(d, bins=n_bins, range=(min_val, max_val), alpha=0.7, histtype='step')[-1][0]
			color = patch.get_edgecolor()

			# Add mean +/- std text
			mean = d.mean()
			ax.text(mean, ax.get_ylim()[1], f'{mean:.{sf}f} +/- {d.std():.{sf+1}f}', fontsize=6, color=color)

		if len(data.shape) > 1:
			for d in data.T:
				plot(d)
		else:
			plot(data)

		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xlabel(label, fontsize=8)
		ax.set_ylabel('Count', fontsize=8)


	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		transcription = sim_data.process.transcription
		ppgpp_id = sim_data.molecule_ids.ppGpp
		rnap_id = sim_data.molecule_ids.full_RNAP
		aa_ids = sim_data.molecule_groups.amino_acids
		ribosome_subunit_ids = [sim_data.molecule_ids.s30_full_complex, sim_data.molecule_ids.s50_full_complex]
		uncharged_trna_names = transcription.uncharged_trna_names
		charged_trna_names = transcription.charged_trna_names
		aa_from_trna = transcription.aa_from_trna.T
		cistron_data = transcription.cistron_data
		if sim_data.external_state.current_timeline_id:
			timeline = sim_data.external_state.saved_timelines[
				sim_data.external_state.current_timeline_id]
		else:
			timeline = []
		rna_fractions = ['is_mRNA', 'is_rRNA', 'is_tRNA']
		convert_to_fraction = lambda x: np.vstack([
			x[:, cistron_data[fraction]].sum(1)
			for fraction in rna_fractions
			]).T
		aa_mw = np.array([sim_data.getter.get_mass(aa[:-3]).asNumber(units.fg / units.count) for aa in aa_ids])
		n_transcribed_rnas = len(transcription.rna_data)
		rna_mw = np.concatenate((
			transcription.rna_data['mw'].asNumber(units.fg / units.count),
			transcription.mature_rna_data['mw'].asNumber(units.fg / units.count)
			))
		rna_mw_synth = rna_mw[:n_transcribed_rnas]
		is_mrna = np.concatenate((
			transcription.rna_data['is_mRNA'],
			np.zeros(len(transcription.mature_rna_data), dtype=bool)
			))
		is_mrna_synth = is_mrna[:n_transcribed_rnas]
		is_rrna = np.concatenate((
			transcription.rna_data['is_rRNA'],
			transcription.mature_rna_data['is_rRNA']
			))
		is_rrna_synth = is_rrna[:n_transcribed_rnas]
		is_trna = np.concatenate((
			transcription.rna_data['is_tRNA'],
			transcription.mature_rna_data['is_tRNA']
			))
		is_trna_synth = is_trna[:n_transcribed_rnas]

		cell_paths = self.ap.get_cells(only_successful=True)

		# Load attributes
		unique_molecule_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
		rnap_idx = unique_molecule_ids.index('active_RNAP')
		ribosome_idx = unique_molecule_ids.index('active_ribosome')

		# Load data
		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		def reduce_rna_degradation(mask):
			return lambda x: (x[:, mask] @ rna_mw[mask]).reshape(-1, 1)
		def reduce_rna_synthesis(mask):
			return lambda x: (x[:, mask] @ rna_mw_synth[mask]).reshape(-1, 1)
		axis_sum = lambda x: x.sum(1).reshape(-1, 1)
		time = read_stacked_columns(cell_paths, 'Main', 'time',
			remove_first=True).squeeze() / 60
		time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec',
			remove_first=True).squeeze()
		growth_rate = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
			remove_first=True).squeeze() * 3600
		protein_growth = read_stacked_columns(cell_paths, 'Mass', 'proteinMass',
			fun=growth_function).squeeze() / time_step * 3600
		rna_growth = read_stacked_columns(cell_paths, 'Mass', 'rnaMass',
			fun=growth_function).squeeze() / time_step * 3600
		small_mol_growth = read_stacked_columns(cell_paths, 'Mass', 'smallMoleculeMass',
			fun=growth_function).squeeze() / time_step * 3600
		protein_mass = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
		rna_mass = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
		mrna_mass = read_stacked_columns(cell_paths, 'Mass', 'mRnaMass', remove_first=True).squeeze()
		rrna_mass = read_stacked_columns(cell_paths, 'Mass', 'rRnaMass', remove_first=True).squeeze()
		trna_mass = read_stacked_columns(cell_paths, 'Mass', 'tRnaMass', remove_first=True).squeeze()
		cell_mass = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True).squeeze()
		rna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'count_RNA_degraded',
			remove_first=True, fun=reduce_rna_degradation(slice(None))).squeeze()
		mrna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'count_RNA_degraded',
			remove_first=True, fun=reduce_rna_degradation(is_mrna)).squeeze()
		rrna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'count_RNA_degraded',
			remove_first=True, fun=reduce_rna_degradation(is_rrna)).squeeze()
		trna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'count_RNA_degraded',
			remove_first=True, fun=reduce_rna_degradation(is_trna)).squeeze()
		ribosome_elong_rate = read_stacked_columns(cell_paths, 'RibosomeData', 'effectiveElongationRate',
			remove_first=True).squeeze()
		rnap_elongations = read_stacked_columns(cell_paths, 'RnapData', 'actualElongations',
			remove_first=True).squeeze()
		aas_elongated = read_stacked_columns(cell_paths, 'GrowthLimits', 'aasUsed', remove_first=True)
		ntps_elongated = read_stacked_columns(cell_paths, 'GrowthLimits', 'ntpUsed', remove_first=True)
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True)
		unique_mol_counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
			remove_first=True)
		rna_fraction_prob = read_stacked_columns(cell_paths, 'RnaSynthProb', 'actual_rna_synth_prob_per_cistron',
			remove_first=True, fun=convert_to_fraction)
		(ppgpp_counts, uncharged_trna_counts, charged_trna_counts, aa_counts,
			inactive_rnap_counts, ribosome_subunit_counts) = read_stacked_bulk_molecules(cell_paths,
				([ppgpp_id], uncharged_trna_names, charged_trna_names, aa_ids, [rnap_id], ribosome_subunit_ids),
				remove_first=True)
		excess, synth_fractions, _ = calculate_ribosome_excesses(sim_data, cell_paths)
		mrna_counts = read_stacked_columns(cell_paths, 'RNACounts', 'mRNA_counts', remove_first=True, fun=axis_sum).squeeze()
		rna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna_synthesis(slice(None)))
		mrna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna_synthesis(is_mrna_synth))
		rrna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna_synthesis(is_rrna_synth))
		trna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna_synthesis(is_trna_synth))

		# Derived quantities
		counts_to_molar_squeezed = counts_to_molar.squeeze()
		ppgpp_conc = ppgpp_counts * counts_to_molar_squeezed * 1000
		charged_trna_counts = charged_trna_counts @ aa_from_trna
		uncharged_trna_counts = uncharged_trna_counts @ aa_from_trna
		fraction_charged = charged_trna_counts / (uncharged_trna_counts + charged_trna_counts)
		aa_conc = counts_to_molar * aa_counts
		active_rnap_counts = unique_mol_counts[:, rnap_idx]
		active_ribosome_counts = unique_mol_counts[:, ribosome_idx]
		rnap_elong_rate = rnap_elongations / time_step / active_rnap_counts
		rnap_fraction_active = active_rnap_counts / (active_rnap_counts + inactive_rnap_counts)
		inactive_ribosome_counts = ribosome_subunit_counts.min(1)
		ribosome_fraction_active = active_ribosome_counts / (active_ribosome_counts + inactive_ribosome_counts)
		rnap_conc = counts_to_molar_squeezed * (active_rnap_counts + inactive_rnap_counts) * 1000
		ribosome_conc = counts_to_molar_squeezed * (active_ribosome_counts + inactive_ribosome_counts) * 1000
		rnap_output = counts_to_molar_squeezed * ntps_elongated.sum(1) / time_step
		ribosome_output = counts_to_molar_squeezed * aas_elongated.sum(1) / time_step
		rp_ratio = rna_mass / protein_mass
		aa_mass = aa_counts @ aa_mw
		rpa_ratio = rna_mass / (protein_mass + aa_mass)
		protein_fraction = protein_mass / cell_mass
		rna_fraction = rna_mass / cell_mass
		aa_fraction = (protein_mass + aa_mass) / cell_mass
		mrna_rrna_ratio = mrna_mass / rrna_mass
		mrna_fraction = mrna_mass / rna_mass
		rrna_fraction = rrna_mass / rna_mass
		trna_fraction = trna_mass / rna_mass
		mrna_conc = mrna_counts / cell_mass
		rna_deg_rate = rna_deg_mass / rna_mass / time_step * 3600
		mrna_deg_ratio = mrna_deg_mass / rna_deg_mass
		rrna_deg_ratio = rrna_deg_mass / rna_deg_mass
		trna_deg_ratio = trna_deg_mass / rna_deg_mass
		mrna_fraction_produced = mrna_produced_mass / rna_produced_mass
		rrna_fraction_produced = rrna_produced_mass / rna_produced_mass
		trna_fraction_produced = trna_produced_mass / rna_produced_mass
		unique_time, cell_count = np.unique(time, return_counts=True)
		filtered = set(unique_time[cell_count < cell_count.max() / 2.])

		# Aggregate data for easy subplotting and selection of data
		data = {
			'Growth rate\n(1/hr)': {'y': [growth_rate], 'lim': [0, 2]},
			'RNA growth rate\n(1/hr)': {'y': [rna_growth], 'lim': [0, 2]},
			'Protein growth rate\n(1/hr)': {'y': [protein_growth], 'lim': [0, 2]},
			'Small mol growth rate\n(1/hr)': {'y': [small_mol_growth], 'lim': [0, 2]},
			'RNAP elongation rate\n(nt/s)': {'y': [rnap_elong_rate], 'lim': [50, 70]},
			'RNAP active fraction': {'y': [rnap_fraction_active], 'lim': [0.1, 0.4]},
			'Ribosome elongation rate\n(AA/s)': {'y': [ribosome_elong_rate], 'lim': [0, 25]},
			'Ribosome active fraction': {'y': [ribosome_fraction_active], 'lim': [0.75, 0.9]},
			f'Fraction charged\n{aa_ids[0][:-3]} tRNA': {'y': [fraction_charged[:, 0]], 'lim': [0, 1]},
			f'Fraction charged\n{aa_ids[10][:-3]} tRNA': {'y': [fraction_charged[:, 10]], 'lim': [0, 1]},
			f'{aa_ids[0][:-3]} concentration\n(mM)': {'y': [aa_conc[:, 0]], 'lim': [0, 10]},
			f'{aa_ids[10][:-3]} concentration\n(mM)': {'y': [aa_conc[:, 10]], 'lim': [0, 2]},
			'ppGpp concentration\n(uM)': {'y': [ppgpp_conc], 'lim': [0, 300]},
			'Fraction charged': {'y': [fraction_charged], 'lim': [0, 1.2]},
			'Amino acid concentrations\n(mM)': {'y': [aa_conc], 'lim': [1e-4, 500], 'log': True},
			'RNA fraction\nsynthesis probability': {'y': [rna_fraction_prob], 'lim': [0, 1]},
			'RNA/protein mass fraction\n(with and without free AA)': {'y': [rp_ratio, rpa_ratio], 'lim': [0, 1]},
			'RNA mass fraction': {'y': [rna_fraction], 'lim': [0, 0.15]},
			'Protein mass fraction\n(with and without free AA)': {'y': [protein_fraction, aa_fraction], 'lim': [0, 0.3]},
			'# cells': {'x': unique_time, 'y': [cell_count]},
			'RNAP conc\n(uM)': {'y': [rnap_conc], 'lim': [0, 10]},
			'RNAP output\n(mM NTPs/s)': {'y': [rnap_output], 'lim': [0.03, 0.14]},
			'Ribosome conc\n(uM)': {'y': [ribosome_conc], 'lim': [0, 40]},
			'Ribosome output\n(mM AA/s)': {'y': [ribosome_output], 'lim': [0, 1]},
			'mRNA:rRNA ratio': {'y': [mrna_rrna_ratio], 'lim': [0.03, 0.09]},
			'RNA mass fractions': {'y': [mrna_fraction, rrna_fraction, trna_fraction], 'lim': [0, 1]},
			'RNA deg rate': {'y': [rna_deg_rate], 'lim': [0.2, 0.7]},
			'RNA deg ratio': {'y': [mrna_deg_ratio, rrna_deg_ratio, trna_deg_ratio], 'lim': [0, 1]},
			'Excess ribosome RNA/protein': {'y': [excess], 'lim': [0, 1]},
			'Synthesis fraction rRNA/rProtein/enzymes': {'y': [synth_fractions], 'lim': [0, 1]},
			'mRNA conc\n(count/fg)': {'y': [mrna_conc], 'lim': [0, 15]},
			'RNA mass fraction produced': {'y': [mrna_fraction_produced, rrna_fraction_produced, trna_fraction_produced], 'lim': [0, 1]},
			}
		data_4 = {
			f'{aa_ids[i][:-3]} concentration\n(mM)': {'y': [aa_conc[:, i]]}
			for i in range(len(aa_ids))
			}
		data_5 = {
			'Growth rates (1/hr)': {'y': [growth_rate, rna_growth, protein_growth], 'lim': [0, 2]},
			}
		# Subset of the data to plot for the paper
		paper_2_keys = [
			'Growth rate\n(1/hr)',
			f'{aa_ids[10][:-3]} concentration\n(mM)',
			f'Fraction charged\n{aa_ids[10][:-3]} tRNA',
			'ppGpp concentration\n(uM)',
			'RNA fraction\nsynthesis probability',
			'RNAP elongation rate\n(nt/s)',
			'RNAP active fraction',
			'Ribosome elongation rate\n(AA/s)',
			]
		paper_4_keys = [f'{aa_id[:-3]} concentration\n(mM)' for aa_id in aa_ids]
		paper_6_keys = [
			'RNA growth rate\n(1/hr)',
			'RNA fraction\nsynthesis probability',
			'RNAP elongation rate\n(nt/s)',
			'RNA deg rate',
			'RNAP output\n(mM NTPs/s)',
			'mRNA:rRNA ratio',
			'ppGpp concentration\n(uM)',
			]

		def subplots(filename, data, keys, filtered, rows=None, cols=None,
				trim=False, row_scale=3., col_scale=3., all_x_labeled=True,
				xlim=None, **kwargs):
			# Determine layout
			if rows:
				cols = int(np.ceil(len(keys) / rows))
			elif cols:
				rows = int(np.ceil(len(keys) / cols))
			else:
				rows = int(np.ceil(len(keys) / np.sqrt(len(keys))))
				cols = int(np.ceil(len(keys) / rows))

			# Plot data on subplots
			_, axes = plt.subplots(rows, cols, figsize=(col_scale*cols, row_scale*rows))
			for i, key in enumerate(keys):
				row = i % rows
				col = i // rows
				if rows == 1 and cols == 1:
					axes = np.array([axes])
				if cols == 1:
					ax = axes[row]
				elif rows == 1:
					ax = axes[col]
				else:
					ax = axes[row, col]

				show_x_labels = (row == rows - 1) or all_x_labeled

				entry = data[key]
				x = entry.get('x', time)
				for j, y in enumerate(entry['y']):
					self.plot_time_series(ax, x, y, key,
						timeline if j == 0 else [], filtered,
						log_scale=entry.get('log', False),
						show_x_labels=show_x_labels,
						**kwargs)

				if trim and (lim := entry.get('lim')):
					ax.set_ylim(lim)
				if xlim is not None:
					ax.set_xlim(xlim)

			# Cleanup text if trimming
			if trim:
				for ax in axes.flatten():
					for text in ax.texts:
						text.set_visible(False)

			plt.tight_layout()
			exportFigure(plt, plotOutDir, filename, metadata)
			plt.close('all')

		# Plot all time series data
		subplots(plotOutFileName, data, data.keys(), set(), downsample=1, rows=4)

		# Downsample for less data and better illustrator load
		subplots(f'{plotOutFileName}_downsampled', data, data.keys(), set(), rows=4)

		# Trim axes from all data for easier comparison across runs
		subplots(f'{plotOutFileName}_trimmed', data, data.keys(), set(), trim=True, rows=4)

		# Set axes limits for easier comparison across runs and filter time
		# points without all cells for smoother traces
		subplots(f'{plotOutFileName}_filtered', data, data.keys(), filtered, trim=True, rows=4)

		# Plots specific for figure 2 in paper
		subplots(f'{plotOutFileName}_fig2', data, paper_2_keys, filtered,
			downsample=10, trim=True, cols=1, row_scale=1.5, all_x_labeled=False)
		subplots(f'{plotOutFileName}_fig2_single', data, paper_2_keys, filtered,
			downsample=10, trim=True, cols=1, single_cells=True)

		# Plots specific for figure 4 in paper
		subplots(f'{plotOutFileName}_fig4_single', data_4, paper_4_keys, filtered,
			downsample=10, trim=True, cols=1, single_cells=True, distinct=True,
			plot_options=dict(color='k', alpha=0.6), single_scale=2,
			shade=False, use_hr=True, xlim=(0, 5))

		# Plots specific for figure 5 in paper
		subplots(f'{plotOutFileName}_fig5', data_5, data_5.keys(), filtered,
			downsample=10, trim=True, cols=1, row_scale=2.7, col_scale=2.7)

		# Plots specific for figure 6 in paper
		subplots(f'{plotOutFileName}_fig6', data, paper_6_keys, filtered,
			downsample=10, trim=True, cols=1, row_scale=1.5, all_x_labeled=False)
		subplots(f'{plotOutFileName}_fig6_xlim', data, paper_6_keys, filtered,
			downsample=10, trim=True, cols=1, row_scale=1.5, all_x_labeled=False,
			xlim=(-40, 950))

		# Plot histograms of data
		# TODO: use data dict from above to generalize this to match any changes in time series traces
		_, axes = plt.subplots(4, 6, figsize=(20, 15))
		self.plot_hist(axes[0, 0], growth_rate, 0, 2, 'Growth rate\n(1/hr)')
		self.plot_hist(axes[1, 0], rna_growth, 0, 2, 'RNA growth rate\n(1/hr)')
		self.plot_hist(axes[2, 0], protein_growth, 0, 2, 'Protein growth rate\n(1/hr)')
		self.plot_hist(axes[3, 0], small_mol_growth, 0, 2, 'Small mol growth rate\n(1/hr)')
		self.plot_hist(axes[0, 1], rnap_elong_rate, 0, 100, 'RNAP elongation rate\n(nt/s)', sf=1)
		self.plot_hist(axes[1, 1], rnap_fraction_active, 0, 1, 'RNAP active fraction')
		self.plot_hist(axes[2, 1], ribosome_elong_rate, 0, 25, 'Ribosome elongation rate\n(AA/s)', sf=1)
		self.plot_hist(axes[3, 1], ribosome_fraction_active, 0, 1, 'Ribosome active fraction')
		self.plot_hist(axes[0, 2], fraction_charged[:, 0], 0, 1, f'Fraction charged\n{aa_ids[0][:-3]} tRNA')
		self.plot_hist(axes[1, 2], fraction_charged[:, 10], 0, 1, f'Fraction charged\n{aa_ids[10][:-3]} tRNA')
		self.plot_hist(axes[2, 2], aa_conc[:, 0], 0, 10, f'{aa_ids[0][:-3]} concentration\n(mM)')
		self.plot_hist(axes[3, 2], aa_conc[:, 10], 0, 1, f'{aa_ids[10][:-3]} concentration\n(mM)')
		self.plot_hist(axes[0, 3], ppgpp_conc, 0, 300, 'ppGpp concentration\n(uM)', n_bins=100, sf=1)
		self.plot_hist(axes[1, 3], fraction_charged, 0, 1, 'Fraction charged')
		self.plot_hist(axes[2, 3], aa_conc, 0, 100, 'Amino acid concentrations\n(mM)')  # TODO: handle log scale and variable ranges
		self.plot_hist(axes[3, 3], rna_fraction_prob, 0, 1, 'RNA fraction\nsynthesis probability')
		self.plot_hist(axes[0, 4], rp_ratio, 0, 1, 'RNA/protein\nmass fraction')
		self.plot_hist(axes[1, 4], rna_fraction, 0, 0.15, 'RNA mass fraction')
		self.plot_hist(axes[2, 4], protein_fraction, 0, 0.3, 'Protein mass fraction')
		self.plot_hist(axes[0, 5], rnap_conc, 0, 10, 'RNAP conc\n(uM)')
		self.plot_hist(axes[1, 5], rnap_output, 0, 0.2, 'RNAP output\n(mM NTPs/s)')
		self.plot_hist(axes[2, 5], ribosome_conc, 0, 40, 'Ribosome conc\n(uM)')
		self.plot_hist(axes[3, 5], ribosome_output, 0, 1, 'Ribosome output\n(mM AA/s)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, f'{plotOutFileName}_hist', metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
