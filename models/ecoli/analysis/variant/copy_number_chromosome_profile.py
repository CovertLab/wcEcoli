"""
Figure A (PLAN.md §6.1) — the chromosome-wide gene dosage profile.

Mean gene copy number against replication coordinate, for the control variant
and the highest-burden variant, with the Cooper-Helmstetter prediction overlaid
at each variant's own measured doubling time. The rRNA operons and the new gene
are marked, because the whole position argument is about where those sit
relative to each other.

This is a single-batch analysis: it needs one output directory and says nothing
about insertion position on its own. The cross-position statement is §6.6.

Memory note: `gene_copy_number` has ~4,500 subcolumns. Stacking every time point
across ~300 cells would be tens of GB, so every read reduces each cell to its
time-average via `fun=` before stacking (`read_stacked_columns` applies `fun`
per cell). `new_gene_counts_save_dt.py` OOMed at 96 GB for want of this.

Measured on Batch 1 (P4): origin-proximal genes lose 0.69 copies under burden
against 0.21 at the terminus, a 3.2x asymmetry, and the simulation tracks
Cooper-Helmstetter to a median 4.6% (control) / 6.0% (burdened).
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.variant.dosage_channel_decomposition import _window
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Number of bins across the replichore for the plotted curve. ~4,500 raw points
# is unreadable; the CSV keeps every gene regardless.
N_BINS = 30


def _cell_mean(x):
	"""Reduce one cell's (time, gene) block to a single row of time-averages."""
	return x.mean(axis=0)[None, :]


def _tau_minutes(x):
	"""Reduce one cell's time column to its doubling time in minutes."""
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if len(variants) < 2:
			print('Need at least two variants (control and burden). Found %d.'
				% len(variants))
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		cistron_data = sim_data.process.transcription.cistron_data.struct_array
		replication = sim_data.process.replication

		coord_of = dict(zip(
			cistron_data['id'], cistron_data['replication_coordinate']))
		rrna_of = dict(zip(cistron_data['id'], cistron_data['is_rRNA']))
		new_of = dict(zip(cistron_data['id'], cistron_data['is_new_gene']))

		per_variant = {}
		gene_ids = cistron_ids = None

		for variant in variants:
			cell_paths = self.ap.get_cells(
				variant=[variant], generation=generations)
			if len(cell_paths) == 0:
				print('No cells for variant %d; skipping.' % variant)
				continue

			if gene_ids is None:
				readable = first_cell_with_table(cell_paths, 'RnaSynthProb')
				if readable is None:
					print('No readable RnaSynthProb for variant %d; skipping.'
						% variant)
					continue
				rsp = TableReader(os.path.join(
					readable, 'simOut', 'RnaSynthProb'))
				gene_ids = rsp.readAttribute('gene_ids')
				cistron_ids = rsp.readAttribute('cistron_ids')

			copy_number = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'gene_copy_number',
				ignore_exception=True, fun=_cell_mean)
			taus = read_stacked_columns(
				cell_paths, 'Main', 'time',
				ignore_exception=True, fun=_tau_minutes)
			if copy_number.size == 0 or taus.size == 0:
				print('No readable cells for variant %d; skipping.' % variant)
				continue

			per_variant[variant] = dict(
				copy_number=copy_number.mean(axis=0),
				tau=float(np.mean(taus)),
				n_cells=copy_number.shape[0],
				)
			print('variant %d: %d cells, mean tau %.1f min'
				% (variant, copy_number.shape[0], per_variant[variant]['tau']))

		if len(per_variant) < 2:
			print('Fewer than two usable variants; nothing to compare.')
			return

		coords = np.array(
			[coord_of.get(c, 0) for c in cistron_ids], dtype=float)
		arm = np.where(
			coords >= 0,
			replication.replichore_lengths[0],
			replication.replichore_lengths[1])
		frac = np.abs(coords) / arm
		is_rrna = np.array(
			[bool(rrna_of.get(c, False)) for c in cistron_ids])
		is_new = np.array([bool(new_of.get(c, False)) for c in cistron_ids])

		# Analytical prediction per variant, at that variant's own tau.
		for variant, rec in per_variant.items():
			rec['analytical'] = replication.get_average_copy_number(
				rec['tau'], coords)

		self._write_csv(plotOutDir, plotOutFileName, gene_ids, cistron_ids,
			coords, frac, is_rrna, is_new, per_variant)
		self._report(frac, is_rrna, is_new, per_variant, variants)
		self._plot(plotOutDir, plotOutFileName, metadata, frac, is_rrna,
			is_new, per_variant, variants)

	def _write_csv(self, plotOutDir, plotOutFileName, gene_ids, cistron_ids,
			coords, frac, is_rrna, is_new, per_variant):
		path = os.path.join(plotOutDir, plotOutFileName + '.csv')
		header = ['gene_id', 'cistron_id', 'replication_coordinate',
			'replichore_fraction', 'is_rRNA', 'is_new_gene']
		for variant in sorted(per_variant):
			header += ['cn_sim_v%02d' % variant, 'cn_ch_v%02d' % variant]

		with open(path, 'w', newline='') as handle:
			writer = csv.writer(handle)
			writer.writerow(header)
			for i, gene_id in enumerate(gene_ids):
				row = [gene_id, cistron_ids[i], int(coords[i]),
					'%.5f' % frac[i], int(is_rrna[i]), int(is_new[i])]
				for variant in sorted(per_variant):
					rec = per_variant[variant]
					row += ['%.5f' % rec['copy_number'][i],
						'%.5f' % rec['analytical'][i]]
				writer.writerow(row)
		print('Wrote %s' % path)

	def _report(self, frac, is_rrna, is_new, per_variant, variants):
		control = variants[0]
		burden = variants[-1]
		if control not in per_variant or burden not in per_variant:
			return
		c0 = per_variant[control]['copy_number']
		c1 = per_variant[burden]['copy_number']
		lost = c0 - c1

		near = frac < 0.10
		far = frac > 0.90

		print('')
		print('=' * 62)
		print('CHROMOSOME DOSAGE PROFILE (PLAN.md §6.1)')
		print('=' * 62)
		print('control variant %d (tau %.1f min) vs burden variant %d '
			'(tau %.1f min)' % (control, per_variant[control]['tau'],
				burden, per_variant[burden]['tau']))
		if near.any():
			print('  origin-proximal   (f<0.10, n=%4d): %.3f -> %.3f, '
				'lost %.3f' % (near.sum(), c0[near].mean(), c1[near].mean(),
					lost[near].mean()))
		if far.any():
			print('  terminus-proximal (f>0.90, n=%4d): %.3f -> %.3f, '
				'lost %.3f' % (far.sum(), c0[far].mean(), c1[far].mean(),
					lost[far].mean()))
		if near.any() and far.any() and lost[far].mean() != 0:
			print('  asymmetry: %.2fx more copies lost near the origin'
				% (lost[near].mean() / lost[far].mean()))
		if is_rrna.any():
			print('  rRNA operons (n=%d, f=%.3f): %.3f -> %.3f, lost %.3f'
				% (is_rrna.sum(), frac[is_rrna].mean(), c0[is_rrna].mean(),
					c1[is_rrna].mean(), lost[is_rrna].mean()))
		if is_new.any():
			print('  new gene     (n=%d, f=%.3f): %.3f -> %.3f, lost %.3f'
				% (is_new.sum(), frac[is_new].mean(), c0[is_new].mean(),
					c1[is_new].mean(), lost[is_new].mean()))
		if is_rrna.any() and is_new.any() and lost[is_new].mean() != 0:
			print('  rRNA vs construct asymmetry: %.2fx'
				% (lost[is_rrna].mean() / lost[is_new].mean()))

		# How closely does the simulation follow the analytical relation?
		for variant in sorted(per_variant):
			rec = per_variant[variant]
			resid = np.abs(rec['copy_number'] - rec['analytical'])
			resid = resid / np.maximum(rec['analytical'], 1e-12)
			print('  variant %d vs Cooper-Helmstetter: median %.2f%%, '
				'90th pct %.2f%%' % (variant, 100 * np.median(resid),
					100 * np.percentile(resid, 90)))
		print('=' * 62)

	def _binned(self, frac, values):
		edges = np.linspace(0, 1, N_BINS + 1)
		centres, means = [], []
		for i in range(N_BINS):
			sel = (frac >= edges[i]) & (frac < edges[i + 1])
			if sel.any():
				centres.append((edges[i] + edges[i + 1]) / 2)
				means.append(values[sel].mean())
		return np.array(centres), np.array(means)

	def _plot(self, plotOutDir, plotOutFileName, metadata, frac, is_rrna,
			is_new, per_variant, variants):
		control = variants[0]
		burden = variants[-1]
		if control not in per_variant or burden not in per_variant:
			return

		fig, ax = plt.subplots(figsize=(7.5, 4.8))

		for variant, colour, label in (
				(control, 'C0', 'control'), (burden, 'C1', 'max burden')):
			rec = per_variant[variant]
			x, y = self._binned(frac, rec['copy_number'])
			ax.plot(x, y, color=colour, lw=2,
				label='%s, simulated (tau %.0f min)' % (label, rec['tau']))
			xa, ya = self._binned(frac, rec['analytical'])
			ax.plot(xa, ya, color=colour, lw=1.2, ls='--', alpha=0.65,
				label='%s, Cooper-Helmstetter' % label)

		for sel, marker, name in (
				(is_rrna, 'o', 'rRNA operons'), (is_new, 's', 'new gene')):
			if not sel.any():
				continue
			ax.plot([frac[sel].mean()],
				[per_variant[control]['copy_number'][sel].mean()],
				marker=marker, color='C2', ms=8, ls='none', label=name)
			ax.plot([frac[sel].mean()],
				[per_variant[burden]['copy_number'][sel].mean()],
				marker=marker, color='C2', ms=8, ls='none',
				markerfacecolor='none')

		ax.set_xlabel('position along replichore (0 = oriC, 1 = terC)')
		ax.set_ylabel('mean gene copy number')
		ax.set_title('Gene dosage across the chromosome, control vs burden',
			fontsize='medium')
		ax.legend(fontsize='x-small', frameon=False)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
