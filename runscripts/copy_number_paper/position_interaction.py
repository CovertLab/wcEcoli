#!/usr/bin/env python
"""
Cross-position analysis: does burden flatten the gene-dosage gradient, and by
more than Cooper-Helmstetter alone predicts?

Every registered variant analysis operates on a single output directory, so the
position contrast has nowhere to live. This script takes two or more batch
directories -- one per insertion position -- and produces the interaction table.

Only valid on batches run with PIN_WT_COORDINATE = True. Unpinned, the
normaliser divides position out of the construct's assigned rate and every
gradient below collapses to 1.0 by construction. The script checks the pin and
refuses to report if it is off.

What it measures
----------------
Per position x variant, over the generation-16-onward window, reduced per cell
so that standard errors are across cells rather than across autocorrelated
timesteps:

  tau              measured doubling time, birth to division
  n_copies         realised construct promoter copy number
  basal_prob       basal_prob_ppgpp_synth_prob for the construct. This is the
                   L1-normalised share of ONE transcription unit, already a
                   per-TU quantity. Pinned, it should be position-independent;
                   that is the fix working, and this is the column to check it on
  init_rate        realised RnapData/rnaInitEvent per TIMESTEP
  init_per_copy    init_rate / n_copies -- realised initiations per promoter per
                   timestep. The quantity to compare between positions
  init_total       rnaInitEvent summed over the whole cell cycle. DO NOT compare
                   this between batches, or ratios of it: it carries the cycle
                   length, and tau differs by up to 1.44x across a burden
                   ladder. Kept only because it is the natural per-cell total
  aborted_frac     incomplete_transcription_event / rnaInitEvent -- transcripts
                   killed by a replisome sweeping the parent domain
  protein          construct monomer counts
  n_ch             Cooper-Helmstetter at each cell's OWN measured tau and the
                   construct's REAL coordinate

Note n_ch uses the real coordinate, not the pinned one. The pin governs what
promoter strength the model assigns; realised copy number still follows the
gene's actual position, and that is what n_ch predicts.

On per-cell-cycle quantities
----------------------------
Any column accumulated over a cell cycle is proportional to that cycle's length.
Comparing such a column, or a per-copy version of it, between two positions
mixes the positional effect with the ratio of doubling times -- which on the
exp-8 ladder runs from 1.00x at variant 1 to 1.44x at variant 7 and manufactures
a trend that looks like a burden-dependent positional effect. Use init_rate and
init_per_copy for anything cross-batch.

What it answers
---------------
The dosage gradient between the extreme positions, at the bottom and the top of
the burden ladder. Predictions to beat, from Cooper-Helmstetter at the exp-8
ladder endpoints (tau 54 -> 87 min):

  P1/P6 gradient        1.594x -> 1.336x
  flattening            16.2% if tau is position-independent
  differential loss     14.5 percentage points, P1 minus P6

If tau is itself position-dependent -- P1 carries more construct, so it should
run slower than P6 at matched translation efficiency -- the flattening exceeds
16.2%. Sizing that against the ladder's own endpoints puts it at 18-21%. That
excess is the only part of this that Cooper-Helmstetter composed with the
measured burden curve does not already give you, and it is small; the honest
framing is that this measures whether the composition holds, not that the
composition is beaten.

aborted_frac is expected to be FLAT across positions. Collisions scale with
copy number and so does output, so the fraction cancels to first order. Stating
that in advance is what makes a flat result reportable rather than a null.

Usage
-----
	export PYTHONPATH="$(pwd):$PYTHONPATH"
	python runscripts/copy_number_paper/position_interaction.py \
		--label P4 --label P1 --label P6 \
		--out out/files_for_claude/copy_number_plan/results \
		out/<p4_batch> out/<p1_batch> out/<p6_batch>

All flags first, all directories last, in matching order. Interleaving them --
`--label P4 out/<p4> --label P1 out/<p1>` -- does NOT parse: sim_dirs is a greedy
nargs='+' positional, so it swallows the first directory and argparse then
rejects the second as an unrecognised argument.

Labels are optional; without them the batch directory name is used. Budget at
least three hours for two batches. Writes position_interaction.csv and prints the
summary.
"""

import argparse
import csv
import os
import pickle

import numpy as np

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (first_cell_with_table,
	read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants

# Matches new_gene_dosage_normalisation and dosage_channel_decomposition.
# Induction is at generation 8; this drops the pre-induction generations and
# the settling that follows them.
IGNORE_FIRST_N_GENS = 16

# Ladder endpoints used for the gradient. Variant 0 is the knockout and has no
# construct at all; variant 1 is full expression at zero translation
# efficiency, which transcribes but imposes no ribosome load.
MIN_BURDEN_VARIANT = 1
MAX_BURDEN_VARIANT = 7


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


def _mean_of(idx):
	"""Per cell: sum the given subcolumns, then average over time."""
	def fn(x):
		return np.array([[float(np.mean(x[:, idx].sum(axis=1)))]])
	return fn


def _total_of(idx):
	"""Per cell: sum the given subcolumns over both time and index."""
	def fn(x):
		return np.array([[float(np.sum(x[:, idx]))]])
	return fn


def _construct_tu_ids(sim_data):
	"""The construct's transcription unit ids, matched through the TU map."""
	cistron_data = sim_data.process.transcription.cistron_data.struct_array
	new_cistrons = sorted(
		set(cistron_data[cistron_data['is_new_gene']]['id'].tolist()))
	if not new_cistrons:
		return [], []

	transcription = sim_data.process.transcription
	rna_data = transcription.rna_data.struct_array
	matrix = transcription.cistron_tu_mapping_matrix.toarray()
	cistron_ids = list(cistron_data['id'])
	rows = [cistron_ids.index(c) for c in new_cistrons]
	tu_mask = matrix[rows, :].sum(axis=0) > 0
	return [str(i) for i in rna_data['id'][tu_mask]], new_cistrons


def _construct_monomer_ids(sim_data, cistron_ids):
	"""Monomer ids for the construct's cistrons."""
	monomer_data = sim_data.process.translation.monomer_data.struct_array
	mapping = dict(zip(monomer_data['cistron_id'], monomer_data['id']))
	return [str(mapping[c]) for c in cistron_ids if c in mapping]


def _coordinates(sim_data, tu_ids):
	"""Return (real, wildtype) oriC-relative coordinates for the construct."""
	rna_data = sim_data.process.transcription.rna_data
	ids = list(rna_data['id'])
	idx = [ids.index(t) for t in tu_ids if t in ids]
	real = float(np.mean([rna_data['replication_coordinate'][i] for i in idx]))
	wt = float(np.mean(
		[rna_data['wt_replication_coordinate'][i] for i in idx]))
	return real, wt


def _load_variant_sim_data(ap, variant, sim_dir, label):
	"""Load one variant's modified sim_data, falling back to the base ParCa.

	Returns (sim_data, source) where source names which pickle was read, since
	only the modified one can show PIN_WT_COORDINATE.
	"""
	try:
		path = ap.get_variant_kb(variant)
		with open(path, 'rb') as f:
			return pickle.load(f), 'simData_Modified.cPickle'
	except (AssertionError, OSError, IndexError) as exc:
		print('  %s: could not read variant %d modified sim_data (%s); falling '
			'back to the base ParCa output. THE PIN CANNOT BE VERIFIED from '
			'that file, so a pin check below is not meaningful.'
			% (label, variant, exc))

	base = os.path.join(
		sim_dir, constants.KB_DIR, constants.SERIALIZED_SIM_DATA_FILENAME)
	try:
		with open(base, 'rb') as f:
			return pickle.load(f), 'simData.cPickle (base ParCa)'
	except OSError as exc:
		print('  %s: no readable sim_data at all (%s); skipping.' % (label, exc))
		return None, None


def _measure(ap, variant, generations, tu_ids, monomer_ids):
	"""Per-cell quantities for one variant of one batch."""
	cell_paths = ap.get_cells(variant=[variant], generation=generations)
	if len(cell_paths) == 0:
		return None

	# Every read below has to cover the same cells. read_stacked_columns skips
	# a cell whose table is unreadable rather than failing, so reading Main,
	# RnaSynthProb and RnapData off the full path list would let a gappy cell
	# drop out of some reads and not others -- different row counts, and an
	# n_cells that does not describe the other columns. Restrict up front.
	cell_paths = np.array([p for p in cell_paths
		if all(os.path.exists(os.path.join(p, 'simOut', t, 'attributes.json'))
			for t in ('Main', 'RnaSynthProb', 'RnapData'))])
	if len(cell_paths) == 0:
		print('  variant %d: no cell has all three listeners; skipping.'
			% variant)
		return None

	synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
	if synth_cell is None:
		print('  variant %d: no cell with a readable RnaSynthProb; skipping.'
			% variant)
		return None
	synth_ids = TableReader(os.path.join(
		synth_cell, 'simOut', 'RnaSynthProb')).readAttribute('rnaIds')

	rnap_cell = first_cell_with_table(cell_paths, 'RnapData')
	rnap_ids = TableReader(os.path.join(
		rnap_cell, 'simOut', 'RnapData')).readAttribute('rnaIds') \
		if rnap_cell else []

	synth_idx = np.array([synth_ids.index(t) for t in tu_ids
		if t in synth_ids])
	rnap_idx = np.array([rnap_ids.index(t) for t in tu_ids if t in rnap_ids])
	if synth_idx.size == 0:
		return None

	basal = read_stacked_columns(cell_paths, 'RnaSynthProb',
		'basal_prob_ppgpp_synth_prob', ignore_exception=True,
		fun=_mean_of(synth_idx))
	copies = read_stacked_columns(cell_paths, 'RnaSynthProb',
		'promoter_copy_number', ignore_exception=True,
		fun=_mean_of(synth_idx))
	taus = read_stacked_columns(cell_paths, 'Main', 'time',
		ignore_exception=True, fun=_tau_minutes)
	if basal.size == 0 or copies.size == 0 or taus.size == 0:
		return None

	row = dict(
		variant=variant,
		n_cells=int(taus.size),
		tau=float(np.mean(taus)), tau_sem=_sem(taus),
		n_copies=float(np.mean(copies)), n_copies_sem=_sem(copies),
		basal_prob=float(np.mean(basal)), basal_prob_sem=_sem(basal),
		)

	# basal_prob is already a per-transcription-unit share, so dividing it by
	# copy number does NOT give a per-copy promoter strength -- it gives a
	# matched share divided by an unmatched copy number, which reads as a pin
	# failure when the pin is working. Kept because it is occasionally the
	# quantity wanted, but check the pin on basal_prob above.
	with np.errstate(divide='ignore', invalid='ignore'):
		per_copy = np.where(copies > 0, basal / copies, np.nan)
	row.update(
		basal_prob_over_copies=float(np.nanmean(per_copy)),
		basal_prob_over_copies_sem=_sem(per_copy),
		)

	if rnap_idx.size:
		# Per timestep, which is what compares across batches. The per-cell-cycle
		# total is kept alongside but must not be ratioed between positions.
		rates = read_stacked_columns(cell_paths, 'RnapData', 'rnaInitEvent',
			ignore_exception=True, fun=_mean_of(rnap_idx))
		inits = read_stacked_columns(cell_paths, 'RnapData', 'rnaInitEvent',
			ignore_exception=True, fun=_total_of(rnap_idx))
		row.update(init_rate=float(np.mean(rates)),
			init_rate_sem=_sem(rates),
			init_total=float(np.mean(inits)),
			init_total_sem=_sem(inits))

		# Per copy, cell by cell -- not mean(rate) / mean(copies), which is a
		# different quantity when copy number varies between cells.
		if rates.size == copies.size:
			with np.errstate(divide='ignore', invalid='ignore'):
				per_promoter = np.where(copies > 0, rates / copies, np.nan)
			row.update(
				init_per_copy=float(np.nanmean(per_promoter)),
				init_per_copy_sem=_sem(per_promoter),
				)

		# Registered in tableAppend and indexed against RnapData's own rnaIds,
		# unlike RnaSynthProb/total_rna_init which is set on the listener but
		# never registered and therefore unreadable.
		aborted = read_stacked_columns(cell_paths, 'RnapData',
			'incomplete_transcription_event', ignore_exception=True,
			fun=_total_of(rnap_idx))
		if aborted.size == inits.size and aborted.size:
			with np.errstate(divide='ignore', invalid='ignore'):
				frac = np.where(inits > 0, aborted / inits, np.nan)
			row.update(
				aborted_total=float(np.mean(aborted)),
				aborted_frac=float(np.nanmean(frac)),
				aborted_frac_sem=_sem(frac),
				)

	mono_cell = first_cell_with_table(cell_paths, 'MonomerCounts')
	if monomer_ids and mono_cell:
		try:
			all_monomers = TableReader(os.path.join(
				mono_cell, 'simOut', 'MonomerCounts')
				).readAttribute('monomerIds')
			mono_idx = np.array([all_monomers.index(m) for m in monomer_ids
				if m in all_monomers])
			if mono_idx.size:
				protein = read_stacked_columns(cell_paths, 'MonomerCounts',
					'monomerCounts', ignore_exception=True,
					fun=_mean_of(mono_idx))
				row.update(protein=float(np.mean(protein)),
					protein_sem=_sem(protein))
		except Exception as exc:  # noqa: BLE001 - protein is a nice-to-have
			print('  could not read MonomerCounts (%s)' % exc)

	return row


def analyse_batch(sim_dir, label):
	"""Every variant of one batch, plus the construct's coordinates."""
	ap = AnalysisPaths(sim_dir, variant_plot=True)
	if ap.n_generation <= IGNORE_FIRST_N_GENS:
		print('%s: only %d generations, fewer than the %d-generation burn-in.'
			% (label, ap.n_generation, IGNORE_FIRST_N_GENS))
		return None
	generations = np.arange(IGNORE_FIRST_N_GENS, ap.n_generation)

	# PIN_WT_COORDINATE is applied by the variant function, so it exists only in
	# a variant's simData_Modified.cPickle. Reading the base ParCa output --
	# kb/simData.cPickle -- can never show the pin, and the guard in report()
	# then refuses to report gradients on correctly pinned batches. Load the
	# sim_data the simulations actually ran with instead.
	variants = sorted(set(ap.get_variants()))
	pin_variant = (MIN_BURDEN_VARIANT if MIN_BURDEN_VARIANT in variants
		else variants[0])
	sim_data, source = _load_variant_sim_data(ap, pin_variant, sim_dir, label)
	if sim_data is None:
		return None

	tu_ids, cistron_ids = _construct_tu_ids(sim_data)
	if not tu_ids:
		print('%s: no new gene transcription units; skipping.' % label)
		return None
	monomer_ids = _construct_monomer_ids(sim_data, cistron_ids)
	real_coord, wt_coord = _coordinates(sim_data, tu_ids)
	print('  %s: coordinates read from %s (variant %d): real=%.0f wt=%.0f'
		% (label, source, pin_variant, real_coord, wt_coord))

	get_n_ch = sim_data.process.replication.get_average_copy_number
	rows = []
	for variant in variants:
		row = _measure(ap, variant, generations, tu_ids, monomer_ids)
		if row is None:
			continue
		# Cooper-Helmstetter at this variant's own measured tau, evaluated at
		# the construct's REAL coordinate: the pin governs assigned promoter
		# strength, not where the gene physically sits.
		row.update(
			position=label,
			real_coordinate=real_coord,
			wt_coordinate=wt_coord,
			n_ch=float(get_n_ch(row['tau'], np.array([real_coord]))),
			)
		row['n_over_n_ch'] = (row['n_copies'] / row['n_ch']
			if row['n_ch'] else float('nan'))
		rows.append(row)

	return dict(label=label, rows=rows, real_coord=real_coord,
		wt_coord=wt_coord)


def _gradient(batches, variant):
	"""Copy-number ratio between the first and last batch at one variant."""
	def copies(batch):
		for row in batch['rows']:
			if row['variant'] == variant:
				return row['n_copies']
		return None

	first, last = copies(batches[0]), copies(batches[-1])
	if not first or not last:
		return None
	return first / last


def report(batches):
	"""Print the interaction, which is the point of the script."""
	print('\n%s\nConstruct coordinates\n%s' % ('=' * 78, '=' * 78))
	print('  %-8s %14s %14s %10s' % ('position', 'real', 'wildtype', 'pinned'))
	pinned_values = set()
	for batch in batches:
		pinned = batch['real_coord'] != batch['wt_coord']
		pinned_values.add(round(batch['wt_coord']))
		print('  %-8s %14.0f %14.0f %10s'
			% (batch['label'], batch['real_coord'], batch['wt_coord'],
				'yes' if pinned else 'no (reference or unpinned)'))

	if len(pinned_values) > 1:
		print('\n  REFUSING TO REPORT GRADIENTS. The batches do not share a '
			'wildtype\n  coordinate (%s), so per-copy promoter strength '
			'differs between\n  positions and the normaliser has divided '
			'position out. Either these\n  batches predate PIN_WT_COORDINATE, '
			'or REFERENCE_WT_COORDINATE was changed\n  between them, and they '
			'cannot be compared this way.'
			% sorted(pinned_values))
		return

	print('\n%s\nPer position x variant\n%s' % ('=' * 78, '=' * 78))
	print('  %-6s %4s %6s %8s %9s %12s %10s %9s %8s'
		% ('pos', 'var', 'cells', 'tau', 'copies', 'basal_prob', 'init/copy',
			'n_ch', 'abort'))
	for batch in batches:
		for row in batch['rows']:
			print('  %-6s %4d %6d %8.1f %9.3f %12.6g %10s %9.3f %8s'
				% (row['position'], row['variant'], row['n_cells'],
					row['tau'], row['n_copies'], row['basal_prob'],
					('%.4f' % row['init_per_copy'])
						if 'init_per_copy' in row else '--',
					row['n_ch'],
					('%.4f' % row['aborted_frac'])
						if 'aborted_frac' in row else '--'))
	print('\n  basal_prob is the per-TU assigned share and should be matched '
		'between\n  positions if the pin took. init/copy is realised '
		'initiations per promoter\n  per timestep -- per timestep, so it can be '
		'compared across batches whose\n  doubling times differ.')

	if len(batches) < 2:
		return

	lo = _gradient(batches, MIN_BURDEN_VARIANT)
	hi = _gradient(batches, MAX_BURDEN_VARIANT)
	first, last = batches[0]['label'], batches[-1]['label']

	print('\n%s\nThe interaction\n%s' % ('=' * 78, '=' * 78))
	if lo is None or hi is None:
		print('  Need variants %d and %d in every batch.'
			% (MIN_BURDEN_VARIANT, MAX_BURDEN_VARIANT))
		return

	print('  %s/%s copy-number gradient' % (first, last))
	print('    at variant %d (min burden) : %.4f' % (MIN_BURDEN_VARIANT, lo))
	print('    at variant %d (max burden) : %.4f' % (MAX_BURDEN_VARIANT, hi))
	print('    flattening                : %.1f%%' % (100 * (1 - hi / lo)))
	print('\n  Cooper-Helmstetter at fixed tau predicts 16.2%% for P1/P6.')
	print('  More than that is the burden feedback: the position carrying more')
	print('  construct runs slower, so its own gradient flattens further.')


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('sim_dirs', nargs='+', help='batch output directories')
	parser.add_argument('--label', action='append', default=[],
		help='label per directory, in the same order')
	parser.add_argument('--out', default='.', help='where to write the CSV')
	args = parser.parse_args()

	labels = args.label or [
		os.path.basename(os.path.normpath(d)) for d in args.sim_dirs]
	if len(labels) != len(args.sim_dirs):
		parser.error('got %d labels for %d directories'
			% (len(labels), len(args.sim_dirs)))

	batches = []
	for sim_dir, label in zip(args.sim_dirs, labels):
		print('Reading %s (%s)' % (label, sim_dir))
		batch = analyse_batch(sim_dir, label)
		if batch and batch['rows']:
			batches.append(batch)

	if not batches:
		print('Nothing measurable.')
		return

	all_rows = [row for batch in batches for row in batch['rows']]
	fields = sorted({key for row in all_rows for key in row})
	os.makedirs(args.out, exist_ok=True)
	csv_path = os.path.join(args.out, 'position_interaction.csv')
	with open(csv_path, 'w', newline='') as f:
		writer = csv.DictWriter(f, fieldnames=fields)
		writer.writeheader()
		writer.writerows(all_rows)

	report(batches)
	print('\nWrote %s' % csv_path)


if __name__ == '__main__':
	main()
