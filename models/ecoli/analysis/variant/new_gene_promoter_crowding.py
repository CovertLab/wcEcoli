"""
Test whether the new gene's per-copy transcription rate survives burden because
it is pinned at the model's physical initiation ceiling.

Why this exists
---------------
Two explanations for the construct's resilience have been measured and rejected:

  Position       Addendum 1. The dosage normaliser divides position out, so the
                 construct's assigned rate at P1 is 0.79x its rate at P4 -- the
                 same factor its copy number is higher. Not the cause.
  Being
  unregulated    Addendum 2. Among unregulated mRNAs the construct sits at the
                 99.9th percentile: it retains 96% of its per-copy rate where
                 the median unregulated gene retains 53%. Not the cause.

A third candidate was not considered in either. `transcript_initiation.py:269`
caps every promoter's initiation probability at the maximum the RNAP footprint
physically allows:

    max_p = (elongation_rate / footprint) * timestep / n_RNAPs_to_activate

Promoters above `max_p` are pinned there and the rest are rescaled up,
iteratively. The consequence is the point:

    initiations for a pinned promoter
        = max_p * n_RNAPs_to_activate
        = (elongation_rate / footprint) * timestep

`n_RNAPs_to_activate` cancels. **A pinned promoter fires at a fixed physical
rate no matter how few RNA polymerases remain**, so it is structurally immune to
the capacity collapse that takes every other gene to ~0.53.

Three things point this way. The footprint is 50 nt (`footprint_sizes.tsv`), so
at ~46 nt/s the ceiling is near 0.9 initiations per promoter per timestep, and
the construct is measured at 0.985 -> 0.948. Its rate barely moves across a
ladder where everything else halves, which is what a clamped quantity does. And
`NG001_RNA[c]` is hardcoded into a debug list three lines below the cap, then
removed from it on the next line.

What this measures
------------------
Per variant, for the construct's transcription units:

  overcrowded_frac   mean of `RnaSynthProb/tu_is_overcrowded`, i.e. the
                     fraction of timesteps the TU is pinned at the ceiling
  ceiling            mean of `max_p * summed rnaInitEvent`, the physical
                     initiation rate a pinned promoter achieves, per promoter
                     per timestep
  measured_rate      realised `rnaInitEvent / promoter_copy_number`
  ratio              measured_rate / ceiling; ~1.0 means pinned

Genome-wide context is also reported: how many TUs are pinned at all, and where
the construct ranks. If the construct is pinned and almost nothing else is, the
resilience is a property of the cap rather than of the construct's biology.

Reading the result
------------------
  overcrowded_frac high, ratio ~1      -> pinned. The 65-87% dosage share is
                                          measured in a saturated regime and
                                          must be reported as such.
  overcrowded_frac ~0                  -> not the cap. The resilience stays a
                                          known unknown and needs a targeted
                                          experiment rather than another read.

Index-matching warning: `tu_is_overcrowded` and `promoter_copy_number` are
TU-indexed against RnaSynthProb's `rnaIds`; `rnaInitEvent` against RnapData's.
Matched by ID here, never by position.

Emits a CSV alongside the PDF.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_columns)
from wholecell.io.tablereader import TableReader

IGNORE_FIRST_N_GENS = 16

# A TU counts as "pinned" if it is overcrowded in more than this fraction of
# timesteps. Only used for the genome-wide context count.
PINNED_THRESHOLD = 0.5


def _window(n_generation):
	if n_generation <= IGNORE_FIRST_N_GENS:
		return None
	return np.arange(IGNORE_FIRST_N_GENS, n_generation)


def _sem(values):
	"""Standard error of the mean across cells. NaN if fewer than two."""
	v = np.asarray(values, dtype=float).ravel()
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


def _time_mean(x):
	"""Per cell: mean over that cell's timesteps, one row per cell."""
	return x.mean(axis=0, keepdims=True)


def _time_sum(x):
	"""Per cell: total over that cell's timesteps, one row per cell."""
	return x.sum(axis=0, keepdims=True)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			print('Run has only %d generations; using all of them.'
				% self.ap.n_generation)
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)
		tu_ids = self._new_gene_tu_ids(sim_data)
		if not tu_ids:
			print('Could not locate the construct; nothing to test.')
			return
		print('Construct transcription units: %s' % tu_ids)
		classify = self._tu_classifier(sim_data, tu_ids)
		self._report_expected_ceiling(sim_data)

		rows = []
		for variant in variants:
			m = self._measure(variant, generations, tu_ids, classify)
			if m is None:
				print('No usable cells for variant %d; skipping.' % variant)
				continue
			m['variant'] = variant
			rows.append(m)
		if not rows:
			print('Nothing measurable.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._write_pinned_csv(plotOutDir, plotOutFileName, rows)
		self._report_pinned(rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, rows, metadata)

	def _report_expected_ceiling(self, sim_data):
		"""Print the ceiling implied by the constants, for a sanity check."""
		try:
			from wholecell.utils import units
			fp = sim_data.process.transcription.active_rnap_footprint_size
			print('RNAP footprint: %s' % fp)
			print('At ~46 nt/s the ceiling is near %.2f initiations per '
				'promoter per second.'
				% (46.0 / fp.asNumber(units.nt)))
		except Exception as exc:  # noqa: BLE001 - informational only
			print('Could not read the footprint size (%s); the measured '
				'ceiling below is the one that matters.' % exc)

	def _tu_classifier(self, sim_data, construct_tu_ids):
		"""
		Map every transcription unit id to a short class label, so a pinned TU
		can be named as rRNA / RNAP / r-protein rather than left as an opaque
		identifier. Uses the same rna_data flags as
		new_gene_machinery_allocation, so the two cannot disagree.
		"""
		try:
			rna_data = sim_data.process.transcription.rna_data.struct_array
			labels = {}
			construct = set(construct_tu_ids)
			for i, tu in enumerate(rna_data['id']):
				tu = str(tu)
				if tu in construct:
					labels[tu] = 'construct'
				elif rna_data['is_rRNA'][i]:
					labels[tu] = 'rRNA'
				elif rna_data['includes_RNAP'][i]:
					labels[tu] = 'RNAP subunit'
				elif rna_data['includes_ribosomal_protein'][i]:
					labels[tu] = 'r-protein'
				elif not rna_data['is_mRNA'][i]:
					labels[tu] = 'other stable'
				else:
					labels[tu] = 'mRNA'
			return labels
		except Exception as exc:  # noqa: BLE001 - naming is informative only
			print('Could not classify transcription units (%s).' % exc)
			return {}

	def _new_gene_tu_ids(self, sim_data):
		try:
			transcription = sim_data.process.transcription
			cistron_data = transcription.cistron_data.struct_array
			rna_data = transcription.rna_data.struct_array
			new_cistrons = set(
				cistron_data[cistron_data['is_new_gene']]['id'].tolist())
			if not new_cistrons:
				return []
			matrix = transcription.cistron_tu_mapping_matrix.toarray()
			cistron_ids = list(cistron_data['id'])
			rows = [cistron_ids.index(c) for c in sorted(new_cistrons)
				if c in cistron_ids]
			hit = np.asarray(matrix[np.array(rows), :].sum(axis=0)).ravel() > 0
			return [str(i) for i in rna_data['id'][hit]]
		except Exception as exc:  # noqa: BLE001
			print('Could not locate new gene TUs (%s).' % exc)
			return []

	def _measure(self, variant, generations, tu_ids, classify=None):
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
		rnap_cell = first_cell_with_table(cell_paths, 'RnapData')
		if synth_cell is None or rnap_cell is None:
			print('No cell with readable listeners for variant %d.' % variant)
			return None
		synth_ids = TableReader(os.path.join(
			synth_cell, 'simOut', 'RnaSynthProb')).readAttribute('rnaIds')
		rnap_ids = TableReader(os.path.join(
			rnap_cell, 'simOut', 'RnapData')).readAttribute('rnaIds')
		synth_idx = np.array(
			[synth_ids.index(t) for t in tu_ids if t in synth_ids])
		rnap_idx = np.array(
			[rnap_ids.index(t) for t in tu_ids if t in rnap_ids])
		if synth_idx.size == 0:
			print('Construct TUs absent from RnaSynthProb for variant %d.'
				% variant)
			return None

		# Per-TU fraction of timesteps spent pinned, for every TU.
		crowd = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'tu_is_overcrowded',
			ignore_exception=True, fun=_time_mean)
		# max_p and the RNAP pool are scalar columns, so they are read
		# unreduced and multiplied per timestep. They are perfectly
		# anti-correlated by construction, so the mean of the product is not
		# the product of the means -- do not reduce these separately.
		# The total-initiations term is rnaInitEvent summed across TUs, which
		# is n_RNAPs_to_activate. Do not reach for RnaSynthProb/total_rna_init
		# instead: it is set on the listener object but never registered in
		# tableAppend, so it reads back empty and silently NaNs the ceiling.
		# The two reads must cover the same cells. read_stacked_columns skips
		# a cell whose table is unreadable rather than failing
		# (analysis_tools.py:289 continues), so reading the two tables from
		# the full path list lets a gappy cell drop out of one read and not
		# the other. The lengths then diverge and the ceiling is lost for the
		# whole variant -- which is what NaN'd this column on Batches 1 and 2,
		# both of which carry fizzled and blocked cells. Restrict to cells
		# that have both tables before reading either.
		paired_paths = [p for p in cell_paths
			if os.path.exists(os.path.join(
				p, 'simOut', 'RnaSynthProb', 'attributes.json'))
			and os.path.exists(os.path.join(
				p, 'simOut', 'RnapData', 'attributes.json'))]
		max_p = read_stacked_columns(
			paired_paths, 'RnaSynthProb', 'max_p', ignore_exception=True)
		n_act = read_stacked_columns(
			paired_paths, 'RnapData', 'rnaInitEvent', ignore_exception=True,
			fun=lambda x: x.sum(axis=1, keepdims=True))
		copies = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True,
			fun=lambda x: np.array([[float(x[:, synth_idx].sum())]]))
		if crowd.size == 0 or copies.size == 0:
			return None
		ceiling = float('nan')
		if max_p.size and n_act.size and max_p.size == n_act.size:
			ceiling = float(np.mean(max_p.ravel() * n_act.ravel()))
		elif max_p.size != n_act.size:
			print('max_p and summed rnaInitEvent differ in length for variant '
				'%d (%d vs %d); the ceiling is not computed.'
				% (variant, max_p.size, n_act.size))

		init = np.array([[0.0]])
		if rnap_idx.size:
			init = read_stacked_columns(
				cell_paths, 'RnapData', 'rnaInitEvent',
				ignore_exception=True,
				fun=lambda x: np.array([[float(x[:, rnap_idx].sum())]]))

		per_tu = crowd.mean(axis=0)
		# WHICH transcription units are at the cap, not just how many. In rich
		# media nine are pinned at variant 0, before the construct exists, and
		# the count alone cannot say whether the rRNA operons are among them --
		# which decides whether rich's rRNA retention numbers are real.
		order = np.argsort(per_tu)[::-1]
		pinned_rows = []
		for j in order:
			if per_tu[j] <= PINNED_THRESHOLD:
				break
			tu = str(synth_ids[j])
			pinned_rows.append(dict(variant=variant, tu_id=tu,
				tu_class=(classify or {}).get(tu, 'unknown'),
				overcrowded_frac=float(per_tu[j])))
		ng_frac = float(per_tu[synth_idx].mean())
		pinned = per_tu > PINNED_THRESHOLD
		# rank 1 = the most-pinned TU in the genome
		rank = int((per_tu > per_tu[synth_idx].max()).sum() + 1)

		total_copies = float(copies.sum())
		measured = float(init.sum()) / total_copies if total_copies else 0.0
		return dict(
			overcrowded_frac=ng_frac,
			overcrowded_frac_sem=_sem(crowd[:, synth_idx].mean(axis=1)),
			ceiling=ceiling,
			measured_rate=measured,
			ratio=measured / ceiling if ceiling else float('nan'),
			n_pinned_tus=int(pinned.sum()),
			construct_rank=rank,
			n_tus=int(per_tu.size),
			n_cells=int(crowd.shape[0]),
			pinned_rows=pinned_rows,
			)

	def _write_csv(self, plot_out_dir, plot_out_file_name, rows):
		fields = ['variant', 'overcrowded_frac', 'overcrowded_frac_sem',
			'ceiling', 'measured_rate', 'ratio', 'n_pinned_tus',
			'construct_rank', 'n_tus', 'n_cells']
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		with open(path, 'w') as handle:
			w = csv.DictWriter(handle, fieldnames=fields)
			w.writeheader()
			for r in rows:
				w.writerow({k: r.get(k, '') for k in fields})

	def _write_pinned_csv(self, plot_out_dir, plot_out_file_name, rows):
		"""One row per pinned transcription unit per variant."""
		path = os.path.join(
			plot_out_dir, plot_out_file_name + '_pinned_tus.csv')
		fields = ['variant', 'tu_id', 'tu_class', 'overcrowded_frac']
		with open(path, 'w') as handle:
			w = csv.DictWriter(handle, fieldnames=fields)
			w.writeheader()
			for r in rows:
				for p in r.get('pinned_rows', []):
					w.writerow(p)

	def _report_pinned(self, rows):
		"""Name the pinned transcription units, grouped by class."""
		print('\nTranscription units at the initiation cap '
			'(overcrowded in > %.0f%% of timesteps)' % (100 * PINNED_THRESHOLD))
		for r in rows:
			pinned = r.get('pinned_rows', [])
			if not pinned:
				print('  variant %d: none' % r['variant'])
				continue
			counts = {}
			for p in pinned:
				counts[p['tu_class']] = counts.get(p['tu_class'], 0) + 1
			summary = ', '.join('%d %s' % (v, k)
				for k, v in sorted(counts.items(), key=lambda x: -x[1]))
			print('  variant %d: %d pinned  (%s)'
				% (r['variant'], len(pinned), summary))
			for p in pinned[:12]:
				print('      %-18s %-14s %.4f'
					% (p['tu_id'], p['tu_class'], p['overcrowded_frac']))
			if len(pinned) > 12:
				print('      ... and %d more' % (len(pinned) - 12))
		print('\n  If rRNA appears here at low burden and not at high burden, '
			'its early\n  per-copy rate was capped and its apparent retention '
			'is inflated.\n')

	def _report(self, rows):
		print('\nPromoter crowding for the construct')
		print('  %-8s %-14s %-10s %-12s %-8s %-10s %s'
			% ('variant', 'overcrowded', 'ceiling', 'measured', 'ratio',
				'pinned TUs', 'rank'))
		for r in rows:
			print('  %-8d %-14.3f %-10.4f %-12.4f %-8.3f %-10d %d'
				% (r['variant'], r['overcrowded_frac'], r['ceiling'],
					r['measured_rate'], r['ratio'], r['n_pinned_tus'],
					r['construct_rank']))

		expressing = [r for r in rows if r['measured_rate'] > 0]
		if not expressing:
			print('  The construct never transcribes; no verdict.')
			return
		frac = float(np.mean([r['overcrowded_frac'] for r in expressing]))
		ratio = float(np.mean([r['ratio'] for r in expressing]))
		print('\n  mean over expressing variants: overcrowded %.3f, '
			'measured/ceiling %.3f' % (frac, ratio))
		if frac > 0.5 and 0.85 < ratio < 1.15:
			print('  VERDICT: the construct is PINNED at the initiation '
				'ceiling. Its per-copy rate is set by elongation rate and '
				'footprint, not by the RNAP pool, which is why it does not '
				'fall with capacity. The 65-87% dosage share is measured in '
				'a saturated regime and must be reported as such.')
		elif frac < 0.1:
			print('  VERDICT: the construct is NOT pinned. The crowding cap '
				'does not explain its resilience, which stays a known '
				'unknown.')
		else:
			print('  VERDICT: partial. Report the fraction and the ratio; '
				'neither limit holds cleanly.')

	def _plot(self, plot_out_dir, plot_out_file_name, rows, metadata):
		v = [r['variant'] for r in rows]
		fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))

		ax = axes[0]
		ax.plot(v, [r['overcrowded_frac'] for r in rows], 'o-', color='C3')
		ax.axhline(0.5, color='C7', ls='--', lw=1)
		ax.set_ylim(-0.05, 1.05)
		ax.set_xlabel('variant index')
		ax.set_ylabel('fraction of timesteps pinned')
		ax.set_title('Is the construct at the ceiling?')

		ax = axes[1]
		ax.plot(v, [r['ceiling'] for r in rows], 's--', color='C7',
			label='ceiling')
		ax.plot(v, [r['measured_rate'] for r in rows], 'o-', color='C2',
			label='measured')
		ax.set_xlabel('variant index')
		ax.set_ylabel('initiations per copy per step')
		ax.set_title('Measured rate vs the physical cap')
		ax.legend(fontsize=8)

		ax = axes[2]
		ax.plot(v, [r['n_pinned_tus'] for r in rows], 'o-', color='C0')
		ax.set_xlabel('variant index')
		ax.set_ylabel('transcription units pinned')
		ax.set_title('How many genes hit the cap at all')

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
