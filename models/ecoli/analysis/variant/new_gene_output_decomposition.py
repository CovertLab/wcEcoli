"""
Figure D (PLAN.md §6.4) — how much of the construct's production ceiling is
dosage.

Per variant this computes three quantities and plots them against translation
efficiency:

	copy number      n  -- new gene template copies per cell
	output           O  -- new gene transcript initiation events per unit time
	per-copy output  r  -- O / n

Because `O = n * r` is an identity rather than a fit, the panel shows directly
how much of the shortfall at high drive is the cell losing gene copies and how
much is each template being throttled. §6.0 answers that question for the rRNA
operons, which is the decisive one for the biology; this answers it for the
construct itself, which is the one a synthetic biologist asks.

The same TU-level reads as §6.0 are used, so the numbers here and there are
directly comparable. Never mix gene-level and TU-level indexing in a ratio --
`gene_copy_number` is indexed by `gene_ids` and `promoter_copy_number` by
`rnaIds`, and silently mixing them gives a plausible-looking wrong answer.

Single-batch analysis.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.variant.dosage_channel_decomposition import (
	_decompose, _window)
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_columns)
from wholecell.io.tablereader import TableReader

try:
	from models.ecoli.sim.variants.new_gene_burden_ladder import (
		TRL_EFF_VALUES)
except ImportError:  # analysis should still run against other variant types
	TRL_EFF_VALUES = None


def _tau_minutes(x):
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


def _trl_eff_for(variant):
	"""Translation efficiency for a burden-ladder variant index, or None."""
	if TRL_EFF_VALUES is None or variant == 0:
		return 0.0 if variant == 0 else None
	idx = variant - 1
	if 0 <= idx < len(TRL_EFF_VALUES):
		return TRL_EFF_VALUES[idx]
	return None


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		rows = []
		for variant in variants:
			cell_paths = self.ap.get_cells(
				variant=[variant], generation=generations)
			if len(cell_paths) == 0:
				continue

			rnap_cell = first_cell_with_table(cell_paths, 'RnapData')
			synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
			if rnap_cell is None or synth_cell is None:
				print('Variant %d: no cell with readable listeners; skipping.'
					% variant)
				continue
			rnap_ids = TableReader(os.path.join(
				rnap_cell, 'simOut', 'RnapData')).readAttribute('rnaIds')
			synth_ids = TableReader(os.path.join(
				synth_cell, 'simOut', 'RnaSynthProb')).readAttribute('rnaIds')

			tu_ids = self._new_gene_tu_ids(sim_data, rnap_ids, synth_ids)
			if not tu_ids:
				print('Variant %d: construct not TU-indexed; skipping.'
					% variant)
				continue

			rnap_idx = np.array([rnap_ids.index(t) for t in tu_ids])
			synth_idx = np.array([synth_ids.index(t) for t in tu_ids])

			init_events = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True, fun=lambda x: x[:, rnap_idx])
			copy_numbers = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'promoter_copy_number',
				ignore_exception=True, fun=lambda x: x[:, synth_idx])
			taus = read_stacked_columns(
				cell_paths, 'Main', 'time',
				ignore_exception=True, fun=_tau_minutes)
			if init_events.size == 0 or copy_numbers.size == 0:
				continue

			total_init = init_events.sum(axis=1)
			total_copies = copy_numbers.sum(axis=1)
			n = float(np.mean(total_copies))
			valid = total_copies > 0
			r = float(np.mean(total_init[valid] / total_copies[valid])) \
				if valid.any() else 0.0

			rows.append(dict(
				variant=variant,
				trl_eff=_trl_eff_for(variant),
				tau_min=float(np.mean(taus)) if taus.size else float('nan'),
				copy_number=n,
				per_copy_rate=r,
				synth_rate=n * r,
				n_cells=int(init_events.shape[0]),
				))

		if len(rows) < 2:
			print('Fewer than two usable variants; nothing to decompose.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, metadata, rows)

	def _new_gene_tu_ids(self, sim_data, rnap_ids, synth_ids):
		"""Return the construct's TU ids, or [] if it cannot be located."""
		try:
			cistron_data = \
				sim_data.process.transcription.cistron_data.struct_array
			new_cistrons = \
				cistron_data[cistron_data['is_new_gene']]['id'].tolist()
		except Exception as exc:  # noqa: BLE001 - best effort
			print('Could not resolve new gene ids (%s).' % exc)
			return []
		found = []
		for cistron in new_cistrons:
			for candidate in (cistron, '%s[c]' % cistron):
				if candidate in rnap_ids and candidate in synth_ids:
					found.append(candidate)
					break
		return found

	def _write_csv(self, plotOutDir, plotOutFileName, rows):
		fields = ['variant', 'trl_eff', 'tau_min', 'copy_number',
			'synth_rate', 'per_copy_rate', 'n_cells']
		path = os.path.join(plotOutDir, plotOutFileName + '.csv')
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow({k: row[k] for k in fields})
		print('Wrote %s' % path)

	def _report(self, rows):
		print('')
		print('=' * 72)
		print('NEW GENE OUTPUT DECOMPOSITION (PLAN.md §6.4)')
		print('=' * 72)
		print('%-8s %9s %9s %11s %13s %11s' % (
			'variant', 'trl_eff', 'tau', 'copies n', 'per-copy r', 'O = n*r'))
		for row in rows:
			te = '-' if row['trl_eff'] is None else '%.1f' % row['trl_eff']
			print('%-8d %9s %9.1f %11.4f %13.6g %11.6g' % (
				row['variant'], te, row['tau_min'], row['copy_number'],
				row['per_copy_rate'], row['synth_rate']))

		# Decompose the full span, control-with-expression to max burden. The
		# knockout control has no construct output, so the meaningful span
		# starts at the lowest induced rung.
		induced = [r for r in rows if r['synth_rate'] > 0]
		if len(induced) >= 2:
			lo, hi = induced[0], induced[-1]
			result = _decompose(lo['copy_number'], lo['per_copy_rate'],
				hi['copy_number'], hi['per_copy_rate'])
			print('')
			print('span variant %d -> %d (trl_eff %s -> %s):'
				% (lo['variant'], hi['variant'], lo['trl_eff'], hi['trl_eff']))
			print('  output %.6g -> %.6g   (dO = %.6g)'
				% (result['output_control'], result['output_burden'],
					result['d_output']))
			print('  dosage term   %.6g' % result['dosage_term'])
			print('  per-copy term %.6g' % result['per_copy_term'])
			print('  identity residual %.3g (must be ~0)'
				% result['identity_residual'])
			if result['d_output'] != 0:
				print('  dosage share of the construct\'s output loss: %.1f%%'
					% result['dosage_share_pct'])
			print('  NOTE: the ladder varies *translation* efficiency, while '
				'n and r here are both transcription-level, so this span is '
				'a genuine burden decomposition for the construct rather '
				'than an artefact of the sweep. Compare it against the rRNA '
				'figure from §6.0: the construct\'s promoter is not under '
				'stringent control, so its per-copy rate barely moves and '
				'almost all of its loss is dosage, whereas rRNA promoters '
				'are downregulated and theirs is mostly per-copy.')
		print('=' * 72)

	def _plot(self, plotOutDir, plotOutFileName, metadata, rows):
		# x axis: translation efficiency where known, else variant index.
		have_te = all(r['trl_eff'] is not None for r in rows)
		x = [r['trl_eff'] if have_te else r['variant'] for r in rows]
		xlabel = 'translation efficiency' if have_te else 'variant index'

		fig, axes = plt.subplots(1, 3, figsize=(12, 3.9))

		axes[0].plot(x, [r['copy_number'] for r in rows], 'o-', color='C0',
			lw=1.8, ms=5)
		axes[0].set_ylabel('template copies  n')
		axes[0].set_title('Dosage', fontsize='medium')

		axes[1].plot(x, [r['per_copy_rate'] for r in rows], 'o-', color='C1',
			lw=1.8, ms=5)
		axes[1].set_ylabel('initiations per copy  r')
		axes[1].set_title('Per-copy activity', fontsize='medium')

		axes[2].plot(x, [r['synth_rate'] for r in rows], 'o-', color='C2',
			lw=1.8, ms=5)
		axes[2].set_ylabel('output  O = n x r')
		axes[2].set_title('Total output', fontsize='medium')

		for ax in axes:
			ax.set_xlabel(xlabel)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
