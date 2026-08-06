"""
Measure how much of the new gene's per-copy transcription rate is set by the
model's built-in dosage normalisation, rather than emerging from the run.

Why this exists
---------------
`transcript_initiation` does not use a free per-promoter rate. When ppGpp
regulation is enabled (the default), it calls
`transcription.synth_prob_from_ppgpp`, which every timestep computes

	tau        = ln(2) / growth_inferred_from_ppGpp / 60
	n_avg_copy = replication.get_average_copy_number(tau, wt_coordinate)
	factor     = (growth + deg_rate) / n_avg_copy
	basal_prob = normalize(expression * factor)

`n_avg_copy` is the Cooper-Helmstetter mean copy number at the transcription
unit's own chromosomal coordinate. So the per-copy initiation probability the
model hands a gene is inversely proportional to how many copies it *expects*
that gene to have. Realised copies then multiply back in, and the cycle-average
total lands on target.

That is correct bookkeeping for a wild-type gene fitted to a measured mRNA
abundance. It has two consequences for the copy-number study:

  Burden axis. If the ppGpp-inferred tau tracks the real slowdown, n_avg_copy
  falls under burden, `basal_prob` rises to compensate, and the per-copy
  channel of the dosage/per-copy decomposition is partly an artefact of the
  normaliser rather than a measurement of promoter behaviour.

  Position axis. The construct is inserted pre-ParCa, so its
  `wt_replication_coordinate` IS its insertion coordinate and differs between
  position batches. Moving it toward the origin raises n_avg_copy, lowers
  basal_prob by the same factor, and can cancel the dosage advantage in the
  cycle-average.

What this script measures
-------------------------
Per variant, for the new gene's transcription units:

  basal_prob     `RnaSynthProb/basal_prob_ppgpp_synth_prob`, the normaliser's
                 own output, summed over the construct's TUs and time-averaged
  n_copies       realised `RnaSynthProb/promoter_copy_number`
  tau            measured doubling time, birth to division
  n_ch           Cooper-Helmstetter at the construct's wt coordinate and tau
  product        basal_prob * n_ch

`product` is the diagnostic. If the normaliser dominates the per-copy rate,
`product` is flat across variants while `basal_prob` rises. If `basal_prob` is
itself flat, the ppGpp-inferred tau is not tracking the real slowdown and the
burden-axis concern does not apply.

The construct's share of total RNAP initiation is also reported, because the
cancellation argument assumes the construct is a small perturbation on the
global `normalize()` denominator. At high expression factors it may not be.

Cross-batch use
---------------
The position-axis test needs two batches. Run this on each insertion-position
batch and compare `basal_prob` for the construct at the same variant index. If
the ratio matches `n_ch` inverted between the two positions, the normaliser
moved with the gene.

Index-matching warning: `basal_prob_ppgpp_synth_prob` and
`promoter_copy_number` are TU-indexed against the `rnaIds` attribute, while
`rnaInitEvent` is TU-indexed against `RnapData`'s own `rnaIds`. Always match by
ID, never by position.

Emits a CSV alongside the PDF.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Discard early generations, matching dosage_channel_decomposition.
IGNORE_FIRST_N_GENS = 16


def _window(n_generation, half):
	"""Generation indices for the primary (half=0) or late (half=1) window."""
	if n_generation <= IGNORE_FIRST_N_GENS:
		return None
	start = IGNORE_FIRST_N_GENS
	if half:
		start = IGNORE_FIRST_N_GENS + (n_generation - IGNORE_FIRST_N_GENS) // 2
	if start >= n_generation:
		return None
	return np.arange(start, n_generation)


def _tau_minutes(x):
	"""Reduce one cell's time column to its doubling time in minutes."""
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


def _sum_of(idx):
	"""Per cell: sum the given subcolumns, then average over time."""
	def fn(x):
		return np.array([[float(np.mean(x[:, idx].sum(axis=1)))]])
	return fn


def _grand_total(x):
	"""Per cell: sum every subcolumn, then average over time."""
	return np.array([[float(np.mean(x.sum(axis=1)))]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		generations = _window(self.ap.n_generation, 0)
		if generations is None:
			print('Run has only %d generations, fewer than the standing '
				'window start of %d. Using all generations.'
				% (self.ap.n_generation, IGNORE_FIRST_N_GENS))
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		replication = sim_data.process.replication
		rna_data = sim_data.process.transcription.rna_data.struct_array
		new_gene_tu_ids = self._new_gene_tu_ids(sim_data)
		if not new_gene_tu_ids:
			print('No new gene transcription units found in sim_data; '
				'nothing to diagnose.')
			return

		mask = np.isin(rna_data['id'], new_gene_tu_ids)
		wt_coords = rna_data['wt_replication_coordinate'][mask]
		coords = rna_data['replication_coordinate'][mask]
		if not np.array_equal(wt_coords, coords):
			print('Note: wt_replication_coordinate differs from '
				'replication_coordinate for the construct (%s vs %s). The '
				'normaliser uses the wt value.' % (wt_coords, coords))
		# One construct locus, possibly several tandem TUs; they share a
		# coordinate to within a few kb, so the mean is the right summary.
		wt_coord = float(np.mean(wt_coords))

		rows = []
		for variant in variants:
			measured = self._measure(variant, generations, new_gene_tu_ids)
			if measured is None:
				print('No usable cells for variant %d; skipping.' % variant)
				continue
			tau = measured['tau']
			n_ch = float(replication.get_average_copy_number(
				tau, np.array([wt_coord])))
			measured.update(
				variant=variant,
				wt_coordinate=wt_coord,
				n_ch=n_ch,
				product=measured['basal_prob'] * n_ch,
				n_over_n_ch=(measured['n_copies'] / n_ch
					if n_ch else float('nan')),
				)
			rows.append(measured)

		if not rows:
			print('Nothing measurable.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, rows, metadata)

	def _new_gene_tu_ids(self, sim_data):
		"""Return the construct's TU ids, or [] if it cannot be located."""
		try:
			cistron_data = \
				sim_data.process.transcription.cistron_data.struct_array
			new_cistrons = set(
				cistron_data[cistron_data['is_new_gene']]['id'].tolist())
		except Exception as exc:  # noqa: BLE001 - best effort, never fatal
			print('Could not read new gene cistrons (%s).' % exc)
			return []
		if not new_cistrons:
			return []

		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data.struct_array
		try:
			matrix = transcription.cistron_tu_mapping_matrix.toarray()
			cistron_ids = list(cistron_data['id'])
			rows = [cistron_ids.index(c) for c in sorted(new_cistrons)]
			tu_mask = matrix[rows, :].sum(axis=0) > 0
			return [str(i) for i in rna_data['id'][tu_mask]]
		except Exception:  # noqa: BLE001 - fall back to an id-prefix match
			return [str(i) for i in rna_data['id']
				if any(str(i).startswith(c) for c in new_cistrons)]

	def _measure(self, variant, generations, tu_ids):
		"""Time- and seed-averaged quantities for one variant."""
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		sim_out_dir = os.path.join(cell_paths[0], 'simOut')
		synth_ids = TableReader(
			os.path.join(sim_out_dir, 'RnaSynthProb')).readAttribute('rnaIds')
		rnap_ids = TableReader(
			os.path.join(sim_out_dir, 'RnapData')).readAttribute('rnaIds')

		synth_idx = np.array(
			[synth_ids.index(tu) for tu in tu_ids if tu in synth_ids])
		rnap_idx = np.array(
			[rnap_ids.index(tu) for tu in tu_ids if tu in rnap_ids])
		if synth_idx.size == 0:
			print('Construct TUs %s are not present in the RnaSynthProb '
				'listener for variant %d.' % (tu_ids, variant))
			return None

		basal = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'basal_prob_ppgpp_synth_prob',
			ignore_exception=True, fun=_sum_of(synth_idx))
		copies = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True, fun=_sum_of(synth_idx))
		taus = read_stacked_columns(
			cell_paths, 'Main', 'time',
			ignore_exception=True, fun=_tau_minutes)
		if basal.size == 0 or copies.size == 0 or taus.size == 0:
			return None

		share = float('nan')
		if rnap_idx.size:
			ng_init = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True, fun=_sum_of(rnap_idx))
			all_init = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True, fun=_grand_total)
			if ng_init.size and all_init.size:
				total = float(np.mean(all_init))
				if total > 0:
					share = float(np.mean(ng_init)) / total

		return dict(
			basal_prob=float(np.mean(basal)),
			n_copies=float(np.mean(copies)),
			tau=float(np.mean(taus)),
			init_share=share,
			n_cells=int(basal.shape[0]),
			)

	def _write_csv(self, plot_out_dir, plot_out_file_name, rows):
		fields = ['variant', 'tau', 'wt_coordinate', 'basal_prob', 'n_copies',
			'n_ch', 'product', 'n_over_n_ch', 'init_share', 'n_cells']
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		with open(path, 'w') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow({k: row.get(k, '') for k in fields})

	def _report(self, rows):
		"""Print the verdict, so it lands in the analysis log."""
		base, last = rows[0], rows[-1]
		print('\nNew gene dosage-normalisation diagnostic')
		print('  construct wt coordinate: %.0f' % base['wt_coordinate'])
		print('  %-8s %-9s %-12s %-10s %-10s %-10s'
			% ('variant', 'tau', 'basal_prob', 'n_ch', 'product', 'share'))
		for row in rows:
			print('  %-8d %-9.2f %-12.6g %-10.4f %-10.6g %-10.4f'
				% (row['variant'], row['tau'], row['basal_prob'],
					row['n_ch'], row['product'], row['init_share']))

		# Variant 0 is the knockout: its basal_prob is ~0 and would dominate
		# any ratio, so the trend is measured over the expressing variants.
		expressing = [r for r in rows if r['basal_prob'] > 0][1:] or rows[1:]
		if len(expressing) < 2:
			return
		first, final = expressing[0], expressing[-1]
		b_ratio = final['basal_prob'] / first['basal_prob']
		p_ratio = (final['product'] / first['product']
			if first['product'] else float('nan'))
		print('\n  across the expressing ladder:')
		print('    basal_prob  x%.3f' % b_ratio)
		print('    product     x%.3f  (flat => the normaliser sets the '
			'per-copy rate)' % p_ratio)
		if abs(p_ratio - 1) < 0.05 and b_ratio > 1.05:
			print('    VERDICT: normaliser is active on the burden axis. The '
				'per-copy channel of the dosage split is partly its doing.')
		elif abs(b_ratio - 1) < 0.05:
			print('    VERDICT: basal_prob is flat, so the ppGpp-inferred tau '
				'is not tracking. The burden-axis concern does not apply.')
		else:
			print('    VERDICT: partial. Report both ratios; neither limit '
				'holds cleanly.')
		print('    construct share of total initiation: %.3f -> %.3f'
			% (first['init_share'], final['init_share']))
		print('    (the cancellation argument assumes this is small)\n')

	def _plot(self, plot_out_dir, plot_out_file_name, rows, metadata):
		v = [r['variant'] for r in rows]
		fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))

		ax = axes[0]
		ax.plot(v, [r['basal_prob'] for r in rows], 'o-', color='C1')
		ax.set_xlabel('variant index')
		ax.set_ylabel('construct basal_prob\n(ppGpp synth prob)')
		ax.set_title('What the normaliser assigns')

		ax = axes[1]
		ax.plot(v, [r['product'] for r in rows], 'o-', color='C0')
		ax.set_xlabel('variant index')
		ax.set_ylabel('basal_prob $\\times$ $n_{CH}$')
		ax.set_title('Flat here means the normaliser\nsets the per-copy rate')

		ax = axes[2]
		ax.plot(v, [r['n_copies'] for r in rows], 'o-', color='C2',
			label='realised copies')
		ax.plot(v, [r['n_ch'] for r in rows], 's--', color='C7',
			label='Cooper$-$Helmstetter')
		ax.set_xlabel('variant index')
		ax.set_ylabel('construct copy number')
		ax.set_title('Realised vs expected copies')
		ax.legend(fontsize=8)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
