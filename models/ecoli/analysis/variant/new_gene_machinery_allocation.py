"""
Where the shrinking transcription budget goes, and which gene classes lose
template dosage as burden slows the cell.

Why this exists
---------------
The construct-level half of the copy-number figure did not survive the drop to
expression factor 8 -- the construct decomposes like any unregulated gene, 33.2%
dosage against rRNA's 30.5%. What survived is the shared channel: a collapse in
global transcriptional capacity that every gene pays. This analysis measures
that channel directly and asks who pays for it, which the earlier scripts only
answered for rRNA.

The organising observation is chromosomal. Every gene in the transcription and
translation machinery sits near the origin:

  rRNA operons        median f = 0.121   5 of 7 inside f < 0.22
  RNAP core subunits  mean f = 0.144     rpoB/C 0.110, rpoA 0.210. The model's
                                         core is 2 x rpoA + rpoB + rpoC
                                         (complexation_reactions.tsv:106);
                                         omega is in no complex at all and
                                         sigma is not a core subunit
  ribosomal proteins  median f = 0.207   41 of 54 inside f < 0.25, 31 in one
                                         cluster at f = 0.19 - 0.23

against a genome mean of f = 0.5 by construction. Cooper-Helmstetter,
n = 2^(((1-f)C + D)/tau), gives origin-proximal genes more copies at every tau
and makes them lose more of those copies, in absolute terms, when tau rises.
So the machinery that sets transcriptional capacity is exactly the machinery
whose template dosage is most exposed to a slowdown -- a positive feedback in
principle. Whether it is a strong one is what the per-class numbers below
decide.

The medium is the second lever. Rich media is fitted at tau = 25.0 min against
minimal's 44.0 (condition_defs.tsv), so replication rounds overlap far more and
the origin-to-terminus gradient steepens from about 1.65x to about 2.4x. Run
this on a minimal batch and a rich batch and compare.

What this measures
------------------
Per variant, the pools:

  active_rnap, inactive_rnap, total_rnap       RNAP budget
  active_ribosome, inactive_ribosome, total    ribosome budget
  total_init                                   summed RnapData/rnaInitEvent

and per gene class -- rRNA operons, RNAP core subunits, ribosomal proteins,
their union, the construct, a terminus-proximal reference, and the genome as a
whole:

  n_copies       mean RnaSynthProb/promoter_copy_number over the class
  init           mean summed RnapData/rnaInitEvent over the class
  r_per_copy     init / n_copies
  rnap_engaged   standing count of polymerases on the class, from
                 RNACounts/partial_mRNA_counts (or partial_rRNA_counts for
                 rRNA) -- the direct allocation measure
  rnap_portion   rnap_engaged / active_rnap
  f_mean         mean replichore fraction of the class
  n_ch           Cooper-Helmstetter at f_mean and the measured tau

Initiation is a flux and engaged polymerases are a stock; Little's law relates
them through transit time, so they answer different questions and both are
reported. Engaged counts are the honest answer to "where did the RNAP go".

Reading the result
------------------
The feedback claim needs two things to be true, and they are separable:

  1. Machinery classes lose more template dosage than the genome average.
     Read n_copies ratios v1 -> v7 per class against the genome row. This is
     nearly guaranteed by position and is a validation, not a finding.
  2. That dosage loss actually costs them polymerase. Read rnap_portion. If
     the machinery's share of RNAP holds steady while its copy number falls,
     the loop is broken at this step and the dosage loss is being absorbed --
     which is the interesting answer and the one the minimal-media data so far
     suggests.

A rising rnap_portion for rRNA under burden would mean the cell is protecting
ribosome synthesis at the expense of everything else; a falling one would mean
the stringent response is winning. Both are publishable; state which before
looking.

Emits a CSV alongside the PDF.
"""

import csv
import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Matches every other analysis in this study. Induction is at generation 8;
# this drops the pre-induction generations and the settling that follows.
IGNORE_FIRST_N_GENS = 16

# Genes within this fraction of the replichore of the terminus form the
# reference class that the machinery is compared against.
TERMINUS_FRACTION = 0.90


def _window(n_generation):
	"""Generations to analyse: everything past the burn-in."""
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


def _tau_minutes(x):
	"""Reduce one cell's time column to its doubling time in minutes."""
	return np.array([[(x[-1, 0] - x[0, 0]) / 60.0]])


def _sum_over(idx):
	"""Per cell: sum the given subcolumns, then average over time."""
	def fn(x):
		if idx is None or len(idx) == 0:
			return np.array([[0.0]])
		return np.array([[float(np.mean(x[:, idx].sum(axis=1)))]])
	return fn


def _time_mean(x):
	"""Per cell: mean over that cell's timesteps, one row per cell, all TUs."""
	return x.mean(axis=0, keepdims=True)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		generations = _window(self.ap.n_generation)
		if generations is None:
			print('Run has only %d generations, fewer than the %d-generation '
				'burn-in; using all of them.'
				% (self.ap.n_generation, IGNORE_FIRST_N_GENS))
			generations = np.arange(self.ap.n_generation)

		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		classes = self._gene_classes(sim_data)
		if classes is None:
			return
		for name, spec in classes.items():
			print('%-20s %5d TUs   mean f = %.3f'
				% (name, len(spec['tu_ids']), spec['f_mean']))

		rows = []
		for variant in variants:
			row = self._measure(variant, generations, classes, sim_data)
			if row is None:
				print('No usable cells for variant %d; skipping.' % variant)
				continue
			row['variant'] = variant
			rows.append(row)

		if not rows:
			print('Nothing measurable.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows, classes)
		self._report(rows, classes)
		self._plot(plotOutDir, plotOutFileName, rows, classes, metadata)

	def _gene_classes(self, sim_data):
		"""
		Build the transcription-unit id list and mean replichore fraction for
		each gene class.

		rRNA, RNAP subunits and ribosomal proteins come from flags that ParCa
		already sets on rna_data, so the class definitions match the rest of
		the model rather than a hand-curated gene list.
		"""
		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data.struct_array
		cistron_data = transcription.cistron_data.struct_array
		ids = np.array([str(i) for i in rna_data['id']])
		coords = np.asarray(rna_data['replication_coordinate'], dtype=float)

		lengths = sim_data.process.replication.replichore_lengths
		frac = np.where(coords > 0, coords / lengths[0], -coords / lengths[1])

		# The construct is matched through the cistron-to-TU map rather than by
		# id prefix, the same way every other analysis in this study does it.
		new_tu = np.zeros(len(ids), dtype=bool)
		try:
			new_cistrons = sorted(set(
				cistron_data[cistron_data['is_new_gene']]['id'].tolist()))
			if new_cistrons:
				matrix = transcription.cistron_tu_mapping_matrix.toarray()
				cistron_ids = list(cistron_data['id'])
				sel = [cistron_ids.index(c) for c in new_cistrons
					if c in cistron_ids]
				new_tu = np.asarray(
					matrix[np.array(sel), :].sum(axis=0)).ravel() > 0
		except Exception as exc:  # noqa: BLE001 - construct is optional here
			print('Could not locate the construct (%s); continuing without '
				'its row.' % exc)

		is_rnap = np.asarray(rna_data['includes_RNAP'], dtype=bool)
		is_rprot = np.asarray(
			rna_data['includes_ribosomal_protein'], dtype=bool)
		is_rrna = np.asarray(rna_data['is_rRNA'], dtype=bool)

		# These two classes OVERLAP and the rows must not be summed. Both flags
		# are TU-level and E. coli co-transcribes RNAP subunits with ribosomal
		# proteins: rpsMKD-rpoA-rplQ carries rpoA beside four r-proteins, and
		# rplKAJL-rpoBC carries rpoB and rpoC beside four more. With operons on
		# (the default) four of the five RNAP-bearing transcription units are
		# also ribosomal-protein-bearing; only rpoBC is RNAP alone.
		#
		# That is a real feature rather than a nuisance -- a dosage hit to those
		# promoters lands on both machines at once, which tightens the coupling
		# this analysis exists to measure -- but it means machinery_any, the
		# union, is the only row that can be read as a total.
		overlap = int((is_rnap & is_rprot).sum())
		if overlap:
			print('NOTE: %d transcription units are in BOTH the RNAP and '
				'ribosomal-protein classes (co-transcribed operons). Do not '
				'sum those two rows; use machinery_any for a total.' % overlap)

		# The construct is an mRNA, so it sits inside any mRNA aggregate. Both
		# are reported because a rising mRNA share is uninteresting if it is
		# just the construct arriving -- mrna_not_construct is the row that
		# says whether the native mRNA sector really gained.
		is_mrna = np.asarray(rna_data['is_mRNA'], dtype=bool)

		masks = {
			'rrna': is_rrna,
			'rnap_subunits': is_rnap,
			'ribosomal_proteins': is_rprot,
			'machinery_any': is_rrna | is_rnap | is_rprot,
			'mrna': is_mrna,
			'mrna_not_construct': is_mrna & ~new_tu,
			'construct': new_tu,
			'terminus_ref': frac > TERMINUS_FRACTION,
			'genome': np.ones(len(ids), dtype=bool),
		}

		# RNACounts publishes engaged-polymerase counts for mRNA and rRNA only,
		# so the genome row's rnap_portion is mRNA + rRNA and stops short of
		# 1.0 by whatever tRNA and misc ncRNA are using. That residual is
		# recoverable as 1 - genome_rnap_portion and is worth reporting: it was
		# 17.0% of engaged RNAP at variant 0 and fell to 13.7% under burden.

		classes = {}
		for name, mask in masks.items():
			if not mask.any():
				print('Gene class %s is empty; dropping it.' % name)
				continue
			classes[name] = dict(
				tu_ids=[str(i) for i in ids[mask]],
				f_mean=float(np.mean(frac[mask])),
				n_tus=int(mask.sum()),
				)
		return classes or None

	def _measure(self, variant, generations, classes, sim_data):
		"""Pools and per-class allocation for one variant."""
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
		rnap_cell = first_cell_with_table(cell_paths, 'RnapData')
		umc_cell = first_cell_with_table(cell_paths, 'UniqueMoleculeCounts')
		if synth_cell is None or rnap_cell is None or umc_cell is None:
			print('Variant %d is missing a required listener; skipping.'
				% variant)
			return None

		synth_ids = TableReader(os.path.join(
			synth_cell, 'simOut', 'RnaSynthProb')).readAttribute('rnaIds')
		rnap_ids = TableReader(os.path.join(
			rnap_cell, 'simOut', 'RnapData')).readAttribute('rnaIds')

		taus = read_stacked_columns(cell_paths, 'Main', 'time',
			ignore_exception=True, fun=_tau_minutes)
		if taus.size == 0:
			return None
		tau = float(np.mean(taus))
		row = dict(tau=tau, tau_sem=_sem(taus), n_cells=int(taus.size))

		# ---- pools -------------------------------------------------------
		umc = TableReader(os.path.join(
			umc_cell, 'simOut', 'UniqueMoleculeCounts'))
		unique_ids = umc.readAttribute('uniqueMoleculeIds')
		counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts',
			'uniqueMoleculeCounts', ignore_exception=True)
		for label, key in (('active_rnap', 'active_RNAP'),
				('active_ribosome', 'active_ribosome')):
			if key in unique_ids and counts.size:
				row[label] = float(np.mean(counts[:, unique_ids.index(key)]))
			else:
				row[label] = float('nan')

		# Inactive ribosomes are limited by whichever subunit is scarcer, the
		# same convention new_gene_counts_and_impacts uses.
		try:
			(rnap_free, s30, s50) = read_stacked_bulk_molecules(
				cell_paths,
				([sim_data.molecule_ids.full_RNAP],
					[sim_data.molecule_ids.s30_full_complex],
					[sim_data.molecule_ids.s50_full_complex]),
				ignore_exception=True)
			row['inactive_rnap'] = float(np.mean(rnap_free))
			row['inactive_ribosome'] = float(np.mean(np.minimum(s30, s50)))
		except Exception as exc:  # noqa: BLE001 - pools are informative only
			print('Could not read free pools for variant %d (%s).'
				% (variant, exc))
			row['inactive_rnap'] = float('nan')
			row['inactive_ribosome'] = float('nan')
		row['total_rnap'] = row['active_rnap'] + row['inactive_rnap']
		row['total_ribosome'] = row['active_ribosome'] + row['inactive_ribosome']

		# Engaged polymerases, split by transcript type. rRNA has its own
		# column because rRNA transcription units are absent from mRNA_ids.
		engaged = {}
		for column, attribute in (('partial_mRNA_counts', 'mRNA_ids'),
				('partial_rRNA_counts', 'rRNA_ids')):
			try:
				reader = TableReader(os.path.join(
					first_cell_with_table(cell_paths, 'RNACounts'),
					'simOut', 'RNACounts'))
				engaged[column] = list(reader.readAttribute(attribute))
			except Exception:  # noqa: BLE001 - older runs may lack RNACounts
				engaged[column] = []

		# ---- per class ---------------------------------------------------
		# Read each wide column ONCE per variant and slice the class index sets
		# out of the in-memory array, rather than re-reading per class. The
		# per-class version did 7 classes x 3 columns = 21 passes over every
		# cell and was entirely latency-bound on Lustre -- 1h45m for one batch.
		get_n_ch = sim_data.process.replication.get_average_copy_number
		copies_all = read_stacked_columns(cell_paths, 'RnaSynthProb',
			'promoter_copy_number', ignore_exception=True, fun=_time_mean)
		init_all = read_stacked_columns(cell_paths, 'RnapData',
			'rnaInitEvent', ignore_exception=True, fun=_time_mean)
		# Total initiation is that same array summed across every TU, so it
		# does not need its own pass over the cells.
		row['total_init'] = (float(np.mean(init_all.sum(axis=1)))
			if init_all.size else float('nan'))
		# Whether each TU sat at the RNAP-footprint cap. Rich media has a much
		# larger polymerase pool, so max_p (which is inversely proportional to
		# it) is TIGHTER there and native promoters clamp at low burden --- 9
		# transcription units were pinned in rich at variant 0 with no
		# construct present. If the rRNA operons are among them their variant-1
		# per-copy rate is capped and every rRNA retention number is biased
		# toward "rRNA does better". This column settles that per class.
		crowd_all = read_stacked_columns(cell_paths, 'RnaSynthProb',
			'tu_is_overcrowded', ignore_exception=True, fun=_time_mean)

		eng_all = {}
		for column in ('partial_mRNA_counts', 'partial_rRNA_counts'):
			if engaged.get(column):
				eng_all[column] = read_stacked_columns(cell_paths, 'RNACounts',
					column, ignore_exception=True, fun=_time_mean)

		def _slice_mean(data, idx):
			"""Per-cell class MEAN from an already-read (cells x TUs) array.

			Used for fractions, where summing across a class is meaningless.
			"""
			if data is None or not data.size or idx.size == 0:
				return np.array([])
			return data[:, idx].mean(axis=1)

		def _slice(data, idx):
			"""Per-cell class total from an already-read (cells x TUs) array."""
			if data is None or not data.size or idx.size == 0:
				return np.array([])
			return data[:, idx].sum(axis=1)

		for name, spec in classes.items():
			tu_ids = spec['tu_ids']
			s_idx = np.array([synth_ids.index(t) for t in tu_ids
				if t in synth_ids], dtype=int)
			r_idx = np.array([rnap_ids.index(t) for t in tu_ids
				if t in rnap_ids], dtype=int)

			copies = _slice(copies_all, s_idx)
			init = _slice(init_all, r_idx)

			n = float(np.mean(copies)) if copies.size else float('nan')
			v = float(np.mean(init)) if init.size else float('nan')
			row['%s_n' % name] = n
			row['%s_n_sem' % name] = _sem(copies)
			row['%s_init' % name] = v
			row['%s_init_sem' % name] = _sem(init)
			row['%s_r' % name] = v / n if n else float('nan')
			row['%s_n_ch' % name] = float(
				get_n_ch(tau, np.array([spec['f_mean'] * 1.0])))

			crowd = _slice_mean(crowd_all, s_idx)
			row['%s_overcrowded' % name] = (
				float(np.mean(crowd)) if crowd.size else float('nan'))
			row['%s_overcrowded_sem' % name] = _sem(crowd)

			# Engaged polymerases. A class can span both transcript types --
			# machinery_any and genome both do -- so BOTH columns are summed.
			# Picking one by class name silently dropped the rRNA part of any
			# mixed class: machinery_any came out identical to
			# ribosomal_proteins because its 291 rRNA-engaged polymerases were
			# discarded, and the genome row reported the mRNA share (0.42)
			# while reading as though it were the whole genome.
			total_eng, found_any = 0.0, False
			for column, attribute in (('partial_mRNA_counts', 'mRNA_ids'),
					('partial_rRNA_counts', 'rRNA_ids')):
				names = engaged.get(column, [])
				if not names or column not in eng_all:
					continue
				e_idx = np.array([names.index(t) for t in tu_ids
					if t in names], dtype=int)
				if e_idx.size == 0:
					continue
				part = _slice(eng_all[column], e_idx)
				if part.size:
					total_eng += float(np.mean(part))
					found_any = True
			row['%s_rnap_engaged' % name] = (
				total_eng if found_any else float('nan'))
			row['%s_rnap_portion' % name] = (
				row['%s_rnap_engaged' % name] / row['active_rnap']
				if row['active_rnap'] and found_any else float('nan'))
			row['%s_f' % name] = spec['f_mean']

		return row

	def _write_csv(self, plot_out_dir, plot_out_file_name, rows, classes):
		fields = ['variant', 'tau', 'tau_sem', 'n_cells',
			'active_rnap', 'inactive_rnap', 'total_rnap',
			'active_ribosome', 'inactive_ribosome', 'total_ribosome',
			'total_init']
		for name in classes:
			fields += ['%s_%s' % (name, s) for s in
				('f', 'n', 'n_sem', 'n_ch', 'init', 'init_sem', 'r',
					'rnap_engaged', 'rnap_portion', 'overcrowded',
					'overcrowded_sem')]
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow({k: row.get(k, '') for k in fields})

	def _report(self, rows, classes):
		"""Print the loop-relevant comparison, so it lands in the log."""
		first, last = rows[0], rows[-1]
		print('\nMachinery allocation, variant %d -> %d'
			% (first['variant'], last['variant']))
		print('  tau %.1f -> %.1f min' % (first['tau'], last['tau']))
		print('  active RNAP     %.0f -> %.0f   (x%.3f)'
			% (first['active_rnap'], last['active_rnap'],
				last['active_rnap'] / first['active_rnap']
				if first['active_rnap'] else float('nan')))
		print('  active ribosome %.0f -> %.0f   (x%.3f)'
			% (first['active_ribosome'], last['active_ribosome'],
				last['active_ribosome'] / first['active_ribosome']
				if first['active_ribosome'] else float('nan')))
		print('  total initiation %.1f -> %.1f  (x%.3f)'
			% (first['total_init'], last['total_init'],
				last['total_init'] / first['total_init']
				if first['total_init'] else float('nan')))

		print('\n  %-20s %6s %16s %16s' % ('class', 'f', 'copies x', 'RNAP portion'))
		for name in classes:
			n0, n1 = first.get('%s_n' % name), last.get('%s_n' % name)
			p0, p1 = (first.get('%s_rnap_portion' % name),
				last.get('%s_rnap_portion' % name))
			print('  %-20s %6.3f %16s %16s'
				% (name, first.get('%s_f' % name, float('nan')),
					('x%.3f' % (n1 / n0)) if n0 else '--',
					('%.4f -> %.4f' % (p0, p1))
						if p0 and np.isfinite(p0) else '--'))
		print('\n  %-20s %14s %14s' % ('class', 'overcrowded v0', 'overcrowded vN'))
		for name in classes:
			a = first.get('%s_overcrowded' % name)
			b = last.get('%s_overcrowded' % name)
			if a is None or not np.isfinite(a):
				continue
			flag = '  <-- CLAMPED' if max(a, b) > 0.05 else ''
			print('  %-20s %14.5f %14.5f%s' % (name, a, b, flag))
		print('\n  A class clamped at low burden and free at high burden will '
			'look like it\n  RETAINS per-copy rate, because its early rate was '
			'capped. Check this\n  before reading any retention number.')
		print('\n  Dosage loss ordered by position is expected and is a '
			'validation.\n  The finding is whether RNAP portion tracks it or '
			'holds steady.\n')

	def _plot(self, plot_out_dir, plot_out_file_name, rows, classes, metadata):
		show = [n for n in
			('rrna', 'rnap_subunits', 'ribosomal_proteins',
				'mrna_not_construct', 'construct', 'terminus_ref')
			if n in classes]
		variants = [r['variant'] for r in rows]

		fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

		ax = axes[0]
		for name in show:
			base = rows[0].get('%s_n' % name)
			if not base:
				continue
			ax.plot(variants, [r['%s_n' % name] / base for r in rows],
				marker='o', label='%s (f=%.2f)' % (name, rows[0]['%s_f' % name]))
		ax.set_xlabel('variant')
		ax.set_ylabel('copy number, relative to variant %d' % variants[0])
		ax.set_title('Template dosage by gene class')
		ax.legend(fontsize=7)

		ax = axes[1]
		for name in show:
			portions = [r.get('%s_rnap_portion' % name) for r in rows]
			if not np.isfinite(portions).any():
				continue
			ax.plot(variants, portions, marker='o', label=name)
		ax.set_xlabel('variant')
		ax.set_ylabel('share of active RNAP')
		ax.set_title('Where the polymerase goes')
		ax.legend(fontsize=7)

		ax = axes[2]
		ax.plot(variants, [r['active_rnap'] for r in rows], marker='o',
			label='active RNAP')
		ax.plot(variants, [r['active_ribosome'] for r in rows], marker='s',
			label='active ribosome')
		ax.set_xlabel('variant')
		ax.set_ylabel('count')
		ax.set_title('Machinery pools')
		ax.legend(fontsize=7)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
