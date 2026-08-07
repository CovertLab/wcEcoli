"""
Figure C (PLAN.md §6.3) — who pays.

Bins every cistron by |replichore fraction| into deciles and reports, per bin,
the mean copy number in the control and in the highest-burden variant, the
absolute drop and the fractional drop. The rRNA operon bin and the new gene's
bin are called out.

This is the quantitative statement behind "origin-proximal genes lose most".
Figure A shows the shape; this table is what you quote. On Batch 1 the first
decile loses ~0.69 copies against ~0.21 in the last, a 3.2x asymmetry, and the
*fractional* drop is far flatter than the absolute one -- which is the honest
caveat, since a gene near terC starts from a lower base.

Single-batch analysis.
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

N_DECILES = 10


def _cell_mean(x):
	return x.mean(axis=0)[None, :]


def _tau_minutes(x):
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if len(variants) < 2:
			print('Need at least two variants. Found %d.' % len(variants))
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

		control, burden = variants[0], variants[-1]
		measured = {}
		cistron_ids = None

		for variant in (control, burden):
			cell_paths = self.ap.get_cells(
				variant=[variant], generation=generations)
			if len(cell_paths) == 0:
				print('No cells for variant %d.' % variant)
				return
			if cistron_ids is None:
				readable = first_cell_with_table(cell_paths, 'RnaSynthProb')
				if readable is None:
					print('No readable RnaSynthProb for variant %d.' % variant)
					return
				cistron_ids = TableReader(os.path.join(
					readable, 'simOut', 'RnaSynthProb')
					).readAttribute('cistron_ids')
			copy_number = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'gene_copy_number',
				ignore_exception=True, fun=_cell_mean)
			taus = read_stacked_columns(
				cell_paths, 'Main', 'time',
				ignore_exception=True, fun=_tau_minutes)
			if copy_number.size == 0:
				print('No readable cells for variant %d.' % variant)
				return
			measured[variant] = dict(
				cn=copy_number.mean(axis=0),
				tau=float(np.mean(taus)) if taus.size else float('nan'),
				n_cells=copy_number.shape[0])

		coords = np.array(
			[coord_of.get(c, 0) for c in cistron_ids], dtype=float)
		arm = np.where(
			coords >= 0,
			replication.replichore_lengths[0],
			replication.replichore_lengths[1])
		frac = np.abs(coords) / arm
		is_rrna = np.array([bool(rrna_of.get(c, False)) for c in cistron_ids])
		is_new = np.array([bool(new_of.get(c, False)) for c in cistron_ids])

		c0 = measured[control]['cn']
		c1 = measured[burden]['cn']

		edges = np.linspace(0, 1, N_DECILES + 1)
		rows = []
		for i in range(N_DECILES):
			lo, hi = edges[i], edges[i + 1]
			sel = (frac >= lo) & (frac < hi if i < N_DECILES - 1
				else frac <= hi)
			if not sel.any():
				continue
			cn0 = float(c0[sel].mean())
			cn1 = float(c1[sel].mean())
			rows.append(dict(
				bin_low='%.1f' % lo, bin_high='%.1f' % hi,
				n_genes=int(sel.sum()),
				cn_control='%.4f' % cn0,
				cn_burdened='%.4f' % cn1,
				delta='%.4f' % (cn0 - cn1),
				frac_delta='%.4f' % ((cn0 - cn1) / cn0 if cn0 else float('nan')),
				n_rrna=int((sel & is_rrna).sum()),
				has_new_gene=int((sel & is_new).any()),
				))

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows, frac, is_rrna, is_new, c0, c1, control, burden,
			measured)
		self._plot(plotOutDir, plotOutFileName, metadata, rows, frac, is_rrna,
			is_new, c0, c1)

	def _write_csv(self, plotOutDir, plotOutFileName, rows):
		fields = ['bin_low', 'bin_high', 'n_genes', 'cn_control',
			'cn_burdened', 'delta', 'frac_delta', 'n_rrna', 'has_new_gene']
		path = os.path.join(plotOutDir, plotOutFileName + '.csv')
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow(row)
		print('Wrote %s' % path)

	def _report(self, rows, frac, is_rrna, is_new, c0, c1, control, burden,
			measured):
		print('')
		print('=' * 72)
		print('SECTOR SUMMARY — who pays (PLAN.md §6.3)')
		print('=' * 72)
		print('control variant %d (tau %.1f min, %d cells) vs burden variant '
			'%d (tau %.1f min, %d cells)'
			% (control, measured[control]['tau'], measured[control]['n_cells'],
				burden, measured[burden]['tau'], measured[burden]['n_cells']))
		print('')
		print('%-12s %7s %10s %10s %9s %9s' % (
			'decile', 'genes', 'control', 'burdened', 'lost', 'lost %'))
		for row in rows:
			flag = ''
			if int(row['n_rrna']):
				flag += '  <- %s rRNA' % row['n_rrna']
			if int(row['has_new_gene']):
				flag += '  <- new gene'
			print('%-12s %7s %10s %10s %9s %8.1f%%%s' % (
				'%s-%s' % (row['bin_low'], row['bin_high']),
				row['n_genes'], row['cn_control'], row['cn_burdened'],
				row['delta'], 100 * float(row['frac_delta']), flag))
		if rows:
			first, last = rows[0], rows[-1]
			d0, d1 = float(first['delta']), float(last['delta'])
			print('')
			print('  absolute asymmetry first vs last decile: %.2fx'
				% (d0 / d1 if d1 else float('nan')))
			f0 = float(first['frac_delta'])
			f1 = float(last['frac_delta'])
			print('  fractional asymmetry:                    %.2fx'
				% (f0 / f1 if f1 else float('nan')))
			print('  (the fractional figure is the conservative one — genes '
				'near terC start from a lower base)')
		print('=' * 72)

	def _plot(self, plotOutDir, plotOutFileName, metadata, rows, frac,
			is_rrna, is_new, c0, c1):
		if not rows:
			return
		centres = [(float(r['bin_low']) + float(r['bin_high'])) / 2
			for r in rows]
		delta = [float(r['delta']) for r in rows]
		frac_delta = [100 * float(r['frac_delta']) for r in rows]

		fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

		axes[0].bar(centres, delta, width=0.085, color='C0', alpha=0.85)
		axes[0].set_xlabel('replichore decile (0 = oriC)')
		axes[0].set_ylabel('copies lost under burden')
		axes[0].set_title('Absolute dosage loss by sector', fontsize='medium')

		axes[1].bar(centres, frac_delta, width=0.085, color='C1', alpha=0.85)
		axes[1].set_xlabel('replichore decile (0 = oriC)')
		axes[1].set_ylabel('fractional loss (%)')
		axes[1].set_title('Fractional dosage loss by sector', fontsize='medium')

		for ax in axes:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			if is_rrna.any():
				ax.axvline(frac[is_rrna].mean(), color='C2', ls=':', lw=1.5,
					label='rRNA operons')
			if is_new.any():
				ax.axvline(frac[is_new].mean(), color='C3', ls='--', lw=1.5,
					label='new gene')
			ax.legend(fontsize='x-small', frameon=False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
