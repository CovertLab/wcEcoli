"""
Ask whether the new gene's per-copy transcription rate is unusual, or just
what any unregulated mRNA does when the stringent response frees up share.

Why this exists
---------------
Under maximum burden the construct retains ~97% of its per-copy initiation rate
while the cell-wide average falls to ~51%. Three candidate explanations were
raised for that:

  1. The construct is exempt from all regulation. `adjust_new_gene_final_
     expression` asserts that new genes carry no attenuation adjustment and no
     transcription-factor delta_prob.
  2. Its expression weight is set by the variant rather than fitted.
  3. Its exp_ppgpp baseline is 1.2645x its exp_free baseline, so its weight
     rises as ppGpp rises.

Point 3 turned out to be a red herring. `set_ppgpp_expression` gives
exp_ppgpp == exp_free exactly for any gene with zero ppGpp fold change, and
`_normalize_ppgpp_expression` then divides each array by its own sum. Because
ppGpp pulls down rRNA and ribosomal proteins, those sums differ, so *every*
unregulated gene ends up at exp_ppgpp/exp_free = sum(exp_free)/sum(exp_ppgpp).
The construct is not special; it is generic.

Which produces a falsifiable prediction: if the construct's resilience is just
"unregulated gene during a stringent response", then it should sit inside the
distribution of native unregulated mRNAs. All of them receive the same ppGpp
gain, all have mRNA-like `loss` terms, and the construct sits mid-replichore
where the 1/n_CH normaliser is close to the population average.

  Construct inside the distribution  -> the x1.90 advantage over the cell
                                        average is entirely the cell average
                                        containing regulated genes. No caveat
                                        needed on the construct.
  Construct outside                  -> the residue is the variant, and that
                                        is the number to quote.

What this measures
------------------
Per transcription unit, the change in per-copy initiation rate between the
lowest-burden expressing variant and the highest:

    r = rnaInitEvent / promoter_copy_number      (both TU-indexed)
    ratio = r(max burden) / r(low burden)

Transcription units are then classified:

    ppgpp        any constituent cistron carries a non-zero ppGpp fold change
    tf           appears in transcription_regulation.delta_prob with non-zero
                 deltaV
    stable       not is_mRNA (rRNA, tRNA)
    unregulated  an mRNA that is none of the above  <- the comparison set

The construct's percentile within the unregulated-mRNA distribution is the
answer.

Index-matching warning: `rnaInitEvent` is TU-indexed against RnapData's
`rnaIds`, `promoter_copy_number` against RnaSynthProb's `rnaIds`. The two are
matched by ID here, never by position.

Memory note: both listeners have ~4,700 subcolumns and we need all of them, so
each cell is reduced to a single summed row via `fun=` before stacking.

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

# A transcription unit needs at least this many initiation events, in both
# variants, before its ratio means anything.
MIN_EVENTS = 20

# |log2 fold change| above which a gene counts as ppGpp-regulated.
PPGPP_FC_THRESHOLD = 0.1


def _window(n_generation):
	if n_generation <= IGNORE_FIRST_N_GENS:
		return None
	return np.arange(IGNORE_FIRST_N_GENS, n_generation)


def _sum_over_time(x):
	"""Per cell: total over the cell's timesteps, one row per cell."""
	return x.sum(axis=0, keepdims=True)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if len(variants) < 2:
			print('Need at least two variants. Found %d.' % len(variants))
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			print('Run has only %d generations; using all of them.'
				% self.ap.n_generation)
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		measured = {}
		for variant in variants:
			m = self._measure(variant, generations)
			if m is not None:
				measured[variant] = m
		if len(measured) < 2:
			print('Fewer than two usable variants.')
			return

		labels = self._classify(sim_data, measured[variants[-1]]['tu_ids'])
		ng_mask = labels['new_gene']
		if not ng_mask.any():
			print('Could not locate the construct; the comparison needs it.')
			return

		# The construct is knocked out in variant 0, so the low-burden
		# reference is the first variant in which it actually transcribes.
		burden = max(measured)
		low = None
		for v in sorted(measured):
			if v == burden:
				continue
			if measured[v]['init'][ng_mask].sum() > 0:
				low = v
				break
		if low is None:
			print('The construct never transcribes; nothing to compare.')
			return
		print('Comparing variant %d (low burden) against variant %d '
			'(max burden).' % (low, burden))

		lo, hi = measured[low], measured[burden]
		with np.errstate(divide='ignore', invalid='ignore'):
			r_lo = lo['init'] / lo['copies']
			r_hi = hi['init'] / hi['copies']
			ratio = r_hi / r_lo

		usable = (
			(lo['init'] >= MIN_EVENTS) & (hi['init'] >= MIN_EVENTS)
			& np.isfinite(ratio) & (ratio > 0))
		print('%d of %d transcription units clear the %d-event floor in both '
			'variants.' % (usable.sum(), usable.size, MIN_EVENTS))

		self._write_csv(plotOutDir, plotOutFileName, hi['tu_ids'], ratio,
			labels, usable, lo, hi)
		verdict = self._report(ratio, labels, usable, low, burden)
		self._plot(plotOutDir, plotOutFileName, ratio, labels, usable,
			verdict, metadata)

	def _measure(self, variant, generations):
		"""Summed initiation events and promoter copies per TU, one variant."""
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		rnap_cell = first_cell_with_table(cell_paths, 'RnapData')
		synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
		if rnap_cell is None or synth_cell is None:
			print('No cell with readable listeners for variant %d.' % variant)
			return None
		rnap_ids = TableReader(os.path.join(
			rnap_cell, 'simOut', 'RnapData')).readAttribute('rnaIds')
		synth_ids = TableReader(os.path.join(
			synth_cell, 'simOut', 'RnaSynthProb')).readAttribute('rnaIds')

		init = read_stacked_columns(
			cell_paths, 'RnapData', 'rnaInitEvent',
			ignore_exception=True, fun=_sum_over_time)
		copies = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True, fun=_sum_over_time)
		if init.size == 0 or copies.size == 0:
			return None

		# Reindex the initiation array onto the RnaSynthProb ID order, by ID.
		order = np.array([rnap_ids.index(tu) if tu in rnap_ids else -1
			for tu in synth_ids])
		missing = int((order < 0).sum())
		if missing:
			print('Warning: %d TUs present in RnaSynthProb are absent from '
				'RnapData; they are dropped.' % missing)
		keep = order >= 0

		return dict(
			tu_ids=np.array(synth_ids)[keep],
			init=init.sum(axis=0)[order[keep]],
			copies=copies.sum(axis=0)[keep],
			n_cells=init.shape[0],
			)

	def _classify(self, sim_data, tu_ids):
		"""Boolean masks over tu_ids for each regulatory category."""
		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data.struct_array
		cistron_data = transcription.cistron_data.struct_array
		idx = {str(t): i for i, t in enumerate(tu_ids)}
		n = len(tu_ids)

		def blank():
			return np.zeros(n, dtype=bool)

		# --- ppGpp-regulated: any constituent cistron with a real fold change
		ppgpp = blank()
		try:
			matrix = transcription.cistron_tu_mapping_matrix.toarray()
			cistron_ids = list(cistron_data['id'])
			rows = []
			for cid, fc in zip(transcription.ppgpp_regulated_genes,
					transcription.ppgpp_fold_changes):
				if abs(fc) > PPGPP_FC_THRESHOLD and cid in cistron_ids:
					rows.append(cistron_ids.index(cid))
			if rows:
				tu_hit = matrix[np.array(rows), :].sum(axis=0) > 0
				for j, hit in enumerate(np.asarray(tu_hit).ravel()):
					if hit and str(rna_data['id'][j]) in idx:
						ppgpp[idx[str(rna_data['id'][j])]] = True
		except Exception as exc:  # noqa: BLE001 - classification is best effort
			print('Could not classify ppGpp regulation (%s).' % exc)

		# --- TF-regulated: appears in delta_prob with a non-zero value
		tf = blank()
		try:
			dp = sim_data.process.transcription_regulation.delta_prob
			for tu_i, val in zip(dp['deltaI'], dp['deltaV']):
				if val != 0 and 0 <= tu_i < len(rna_data):
					key = str(rna_data['id'][tu_i])
					if key in idx:
						tf[idx[key]] = True
		except Exception as exc:  # noqa: BLE001
			print('Could not classify TF regulation (%s).' % exc)

		# --- stable RNA, and the construct
		stable, new_gene = blank(), blank()
		for j, rid in enumerate(rna_data['id']):
			key = str(rid)
			if key in idx and not rna_data['is_mRNA'][j]:
				stable[idx[key]] = True
		try:
			new_cistrons = set(
				cistron_data[cistron_data['is_new_gene']]['id'].tolist())
			matrix = transcription.cistron_tu_mapping_matrix.toarray()
			cistron_ids = list(cistron_data['id'])
			rows = [cistron_ids.index(c) for c in new_cistrons
				if c in cistron_ids]
			if rows:
				tu_hit = np.asarray(
					matrix[np.array(rows), :].sum(axis=0)).ravel() > 0
				for j, hit in enumerate(tu_hit):
					if hit and str(rna_data['id'][j]) in idx:
						new_gene[idx[str(rna_data['id'][j])]] = True
		except Exception as exc:  # noqa: BLE001
			print('Could not locate the new gene (%s).' % exc)

		unreg = ~ppgpp & ~tf & ~stable & ~new_gene
		return dict(ppgpp=ppgpp, tf=tf, stable=stable, new_gene=new_gene,
			unregulated=unreg)

	def _write_csv(self, plot_out_dir, plot_out_file_name, tu_ids, ratio,
			labels, usable, lo, hi):
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		fields = ['tu_id', 'ratio', 'usable', 'category', 'init_low',
			'init_burden', 'copies_low', 'copies_burden']
		with open(path, 'w') as handle:
			w = csv.DictWriter(handle, fieldnames=fields)
			w.writeheader()
			for i, tu in enumerate(tu_ids):
				if labels['new_gene'][i]:
					cat = 'new_gene'
				elif labels['stable'][i]:
					cat = 'stable'
				elif labels['ppgpp'][i]:
					cat = 'ppgpp'
				elif labels['tf'][i]:
					cat = 'tf'
				else:
					cat = 'unregulated'
				w.writerow(dict(
					tu_id=str(tu), ratio='%.6g' % ratio[i],
					usable=int(usable[i]), category=cat,
					init_low='%.6g' % lo['init'][i],
					init_burden='%.6g' % hi['init'][i],
					copies_low='%.6g' % lo['copies'][i],
					copies_burden='%.6g' % hi['copies'][i]))

	def _report(self, ratio, labels, usable, low, burden):
		"""Print the verdict and return what the plot needs."""
		ref = usable & labels['unregulated']
		ng = usable & labels['new_gene']
		print('\nPer-copy initiation rate, variant %d -> %d' % (low, burden))
		for name in ('unregulated', 'ppgpp', 'tf', 'stable'):
			m = usable & labels[name]
			if m.sum():
				v = ratio[m]
				print('  %-13s n=%-5d median %.3f   IQR %.3f - %.3f'
					% (name, m.sum(), np.median(v),
						np.percentile(v, 25), np.percentile(v, 75)))
		if not ng.any() or not ref.any():
			print('  construct or reference set empty; no verdict.')
			return None

		ng_val = float(np.mean(ratio[ng]))
		ref_v = ratio[ref]
		pct = float((ref_v < ng_val).mean() * 100)
		print('\n  construct     %.3f   -> %.1fth percentile of unregulated '
			'mRNAs' % (ng_val, pct))
		print('  unregulated median %.3f, 10th-90th %.3f - %.3f'
			% (np.median(ref_v), np.percentile(ref_v, 10),
				np.percentile(ref_v, 90)))
		if 10 <= pct <= 90:
			print('  VERDICT: the construct sits inside the unregulated '
				'distribution. Its resilience is what any unregulated mRNA '
				'does; no construct-specific caveat is needed.')
		else:
			print('  VERDICT: the construct is an outlier among unregulated '
				'mRNAs. The residue is the variant, and %.3f vs a median of '
				'%.3f is the number to quote.'
				% (ng_val, np.median(ref_v)))
		return dict(ng=ng_val, pct=pct, ref=ref_v)

	def _plot(self, plot_out_dir, plot_out_file_name, ratio, labels, usable,
			verdict, metadata):
		fig, axes = plt.subplots(1, 2, figsize=(11, 4))

		ax = axes[0]
		ref = ratio[usable & labels['unregulated']]
		if ref.size:
			ax.hist(ref, bins=40, color='C0', alpha=0.75,
				label='unregulated mRNAs')
			ax.axvline(np.median(ref), color='C0', ls='--', lw=1.5,
				label='their median')
		if verdict:
			ax.axvline(verdict['ng'], color='C2', lw=2.5,
				label='the construct')
			ax.set_title('Construct at the %.0fth percentile'
				% verdict['pct'])
		ax.set_xlabel('per-copy initiation rate, burden / low burden')
		ax.set_ylabel('transcription units')
		ax.legend(fontsize=8)

		ax = axes[1]
		names = ['unregulated', 'ppgpp', 'tf', 'stable']
		data, ticks = [], []
		for name in names:
			m = usable & labels[name]
			if m.sum():
				data.append(ratio[m])
				ticks.append('%s\nn=%d' % (name, m.sum()))
		if data:
			ax.boxplot(data, labels=ticks, showfliers=False)
		if verdict:
			ax.axhline(verdict['ng'], color='C2', lw=2,
				label='the construct')
			ax.legend(fontsize=8)
		ax.set_ylabel('per-copy rate ratio')
		ax.set_title('By regulatory category')

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
