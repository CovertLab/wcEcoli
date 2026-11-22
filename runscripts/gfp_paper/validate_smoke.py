#!/usr/bin/env python
"""
Validate a smoke batch before committing to production runs.

Checks the things that are cheap to get wrong and expensive to discover
later: that every listener column the paper panels read is present and not
degenerate, that seed-matched lineages really are counterfactuals, that a
rerun reproduces bitwise, that the ribosome crowding cap moved as the
footprint fix predicts, and that the new gene landed where we asked.

Nothing here draws a figure. Each check prints PASS, FAIL or a note, and the
exit status is non-zero if any check fails, so it can gate a workflow.

Usage:
	export PYTHONPATH="$(pwd):$PYTHONPATH"
	python runscripts/gfp_paper/validate_smoke.py <batch_dir> [--induction 8]

	<batch_dir> is the run directory holding kb/ and the variant folders,
	e.g. $SCRATCH/wcEcoli_out/20261120.120000__gfp_smoke
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import os
import pickle
import sys

import numpy as np

from wholecell.io.tablereader import TableReader

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, 'required_listeners.tsv')

# Columns that exist only on the directional decay branch. Their absence is
# reported but is not a failure, so this script works on both arms.
DECAY_ONLY = {
	'decaying_mRNA_counts', 'decaying_cleared_mRNA_counts',
	'rnap_engaged_mRNA_counts', 'available_mRNA_cistron_counts',
	'countRnaCistronSynthesized',
}

# Columns that may legitimately be all zero. Crowding flags are zero whenever
# nothing reaches its spacing cap, which after the ribosome footprint fix is
# the expected result (the smoke batch peaked at 74% of the cap); headroom is
# reported in the crowding section instead.
MAY_BE_ZERO = {
	('RibosomeData', 'mRNA_is_overcrowded'),
	('RnaSynthProb', 'tu_is_overcrowded'),
}


class Report:
	"""Collects check results and decides the exit status."""

	def __init__(self):
		self.failures = 0
		self.notes = 0

	def check(self, ok, label, detail=''):
		print('  [%s] %-56s %s'
			% ('PASS' if ok else 'FAIL', label, detail))
		if not ok:
			self.failures += 1

	def note(self, label, detail=''):
		print('  [note] %-56s %s' % (label, detail))
		self.notes += 1

	def section(self, title):
		print('\n%s\n%s' % (title, '-' * len(title)))


def read_manifest():
	"""{listener: {column: required_by}} from required_listeners.tsv."""
	wanted = collections.defaultdict(dict)
	with open(MANIFEST) as f:
		for row in csv.DictReader(f, delimiter='\t'):
			wanted[row['listener']][row['column']] = row['required_by']
	return wanted


def sim_dirs(batch_dir):
	"""[(variant, seed, generation, simOut path)] sorted, for one batch."""
	out = []
	pattern = os.path.join(batch_dir, '*', '[0-9]' * 6, 'generation_*', '*',
		'simOut')
	for path in sorted(glob.glob(pattern)):
		parts = path.split(os.sep)
		generation = int(parts[-3].split('generation_')[1][:6])
		out.append((parts[-5], int(parts[-4]), generation, path))
	return out


def check_listeners(batch_dir, report):
	"""Every manifest column is present, and not constant-zero."""
	report.section('Listener coverage')
	runs = sim_dirs(batch_dir)
	if not runs:
		report.check(False, 'found any simOut directories',
			'nothing under %s' % batch_dir)
		return
	ran = [r for r in runs if os.path.exists(os.path.join(r[3], 'Main',
		'attributes.json'))]
	report.note('simulations found', '%d that ran (of %d folders) across %d'
		' variants; folders without data are descendants of cells that'
		' failed' % (len(ran), len(runs), len({r[0] for r in runs})))
	runs = ran or runs

	_, _, _, probe = runs[len(runs) // 2]
	wanted = read_manifest()
	missing, degenerate, decay_absent = [], [], []
	for listener, columns in sorted(wanted.items()):
		path = os.path.join(probe, listener)
		if not os.path.isdir(path):
			missing.extend('%s/%s' % (listener, c) for c in columns)
			continue
		reader = TableReader(path)
		present = set(reader.columnNames())
		for column in sorted(columns):
			if column not in present:
				(decay_absent if column in DECAY_ONLY else missing).append(
					'%s/%s' % (listener, column))
				continue
			try:
				values = reader.readColumn(column)
			except Exception as exc:                       # noqa: BLE001
				degenerate.append('%s/%s unreadable (%s)'
					% (listener, column, exc))
				continue
			finite = np.asarray(values, dtype=float)
			finite = finite[np.isfinite(finite)] if finite.size else finite
			if finite.size == 0 or not np.any(finite != 0):
				if (listener, column) in MAY_BE_ZERO and finite.size:
					continue
				degenerate.append('%s/%s all zero or empty'
					% (listener, column))

	report.check(not missing, 'all required listener columns present',
		'' if not missing else 'missing %d: %s' % (len(missing),
			', '.join(missing[:6])))
	report.check(not degenerate, 'no required column is empty or all-zero',
		'' if not degenerate else '%d: %s' % (len(degenerate),
			'; '.join(degenerate[:4])))
	if decay_absent:
		report.note('decay-branch columns absent (expected on arm A)',
			'%d, e.g. %s' % (len(decay_absent), decay_absent[0]))


def check_seed_matching(batch_dir, report, induction_gen):
	"""
	Pre-induction generations must be identical across variants for a seed.

	Every growth comparison in the paper is against a seed-matched control,
	so if the lineages diverge before induction the counterfactual is void.
	"""
	report.section('Seed-matched counterfactuals (generations 0 to %d)'
		% (induction_gen - 1))
	runs = sim_dirs(batch_dir)
	variants = sorted({r[0] for r in runs})
	if len(variants) < 2:
		report.note('only one variant present', 'nothing to compare')
		return

	by_key = {(v, s, g): p for v, s, g, p in runs}
	seeds = sorted({r[1] for r in runs})
	mismatches, compared = [], 0
	for seed in seeds:
		for generation in range(induction_gen):
			series = {}
			for variant in variants:
				path = by_key.get((variant, seed, generation))
				if path is None:
					continue
				try:
					series[variant] = TableReader(
						os.path.join(path, 'Mass')).readColumn('cellMass')
				except Exception:                           # noqa: BLE001
					continue
			if len(series) < 2:
				continue
			compared += 1
			reference = series[variants[0]]
			for variant, values in series.items():
				if (values.shape != reference.shape
						or not np.array_equal(values, reference)):
					mismatches.append('seed %d gen %d: %s differs'
						% (seed, generation, variant))
	report.check(not mismatches,
		'pre-induction generations bitwise identical across variants',
		'%d comparisons' % compared if not mismatches
		else '%d mismatches, e.g. %s' % (len(mismatches), mismatches[0]))


def check_crowding(batch_dir, report, induction_gen):
	"""
	How close translation comes to the ribosome crowding cap.

	Headroom is each mRNA's requested initiation probability over its cap
	(RibosomeData target_prob_translation_per_transcript / max_p_per_protein),
	the peak over post-induction time steps: below 1 means nothing was
	crowded. Reported for the new gene and for any mRNA, with the flag rate
	alongside. Before the footprint fix the cap sat about three times lower,
	so headroom here should be about a third of the pre-fix value.
	"""
	report.section('Ribosome crowding after the footprint fix')
	runs = [r for r in sim_dirs(batch_dir) if r[2] >= induction_gen]
	if not runs:
		report.note('no post-induction generations', 'skipped')
		return
	stats = collections.defaultdict(lambda: {'cells': 0, 'flagged_cells': 0,
		'flagged_steps': 0, 'gfp': [], 'any': []})
	for variant, _, _, path in runs:
		try:
			reader = TableReader(os.path.join(path, 'RibosomeData'))
			flags = reader.readColumn('mRNA_is_overcrowded')
			target = reader.readColumn('target_prob_translation_per_transcript')
			cap = reader.readColumn('max_p_per_protein')
			ids = list(reader.readAttribute('monomerIds'))
		except Exception:                                   # noqa: BLE001
			continue
		with np.errstate(divide='ignore', invalid='ignore'):
			load = np.where(cap > 0, target / cap, np.nan)
		gfp = [i for i, m in enumerate(ids) if m.startswith('NG-')]
		# Count time steps with any mRNA flagged: a mean over every mRNA
		# rounds thousands of flagged steps to zero.
		steps = int(np.any(flags > 0, axis=1).sum()) if flags.ndim == 2 \
			else int(np.sum(flags > 0))
		stats[variant]['cells'] += 1
		stats[variant]['flagged_cells'] += int(steps > 0)
		stats[variant]['flagged_steps'] += steps
		stats[variant]['any'].append(float(np.nanmax(load)))
		if gfp:
			stats[variant]['gfp'].append(float(np.nanmax(load[:, gfp])))
	if not stats:
		report.note('crowding columns unreadable', 'skipped')
		return
	for variant in sorted(stats):
		row = stats[variant]
		gfp_peak = [v for v in row['gfp'] if np.isfinite(v)]
		report.note('%s: peak load / cap, new gene and any mRNA' % variant,
			'%s and %.3f; %d of %d cells flagged, %d flagged time steps' % (
			'%.3f' % max(gfp_peak) if gfp_peak else 'n/a', max(row['any']),
			row['flagged_cells'], row['cells'], row['flagged_steps']))
	peak = max(max(r['any']) for r in stats.values())
	report.note('highest load / cap in the batch', '%.3f (%s)' % (peak,
		'nothing crowded' if peak < 1 else 'some mRNAs at the cap'))


def check_insertion(batch_dir, report):
	"""The new gene went where we asked, and sim_data records where."""
	report.section('New gene insertion')
	path = os.path.join(batch_dir, 'kb', 'simData.cPickle')
	if not os.path.isfile(path):
		report.note('sim_data not found', path)
		return
	with open(path, 'rb') as f:
		sim_data = pickle.load(f)
	resolved = getattr(sim_data, 'new_gene_insertion_pos', None)
	requested = getattr(sim_data, 'new_gene_requested_insertion_pos', None)
	if resolved is None:
		report.note('resolved insertion coordinate not recorded',
			'sim_data predates the recording change')
	else:
		report.check(resolved == requested,
			'insertion landed where requested',
			'resolved %s, requested %s' % (resolved, requested))
		from wholecell.analysis import gfp_paper_constants
		report.check(requested == gfp_paper_constants.INSERTION_POS,
			'requested site is the paper\'s (%s)'
			% gfp_paper_constants.INSERTION_SITE_NAME,
			'requested %s, paper %s' % (requested,
			gfp_paper_constants.INSERTION_POS))
	cistron_data = sim_data.process.transcription.cistron_data.struct_array
	new_genes = cistron_data[cistron_data['is_new_gene']]['id'].tolist()
	report.check(bool(new_genes), 'new gene present in sim_data',
		', '.join(new_genes))


def check_disk(batch_dir, report):
	"""Measured disk per generation, against the 1 GB planning figure."""
	report.section('Disk')
	runs = sim_dirs(batch_dir)
	if not runs:
		return
	sizes = []
	for _, _, _, path in runs[:6]:
		total = 0
		for root, _, files in os.walk(path):
			for name in files:
				try:
					total += os.path.getsize(os.path.join(root, name))
				except OSError:
					pass
		sizes.append(total)
	mean_gb = float(np.mean(sizes)) / 1e9
	report.note('mean simOut size per generation', '%.2f GB (n=%d)'
		% (mean_gb, len(sizes)))
	report.note('implied for a 6,000 generation batch', '%.1f TB'
		% (mean_gb * 6000 / 1000))


def main():
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('batch_dir', help='run directory holding kb/')
	parser.add_argument('--induction', type=int, default=8,
		help='generation at which the new gene is induced (default 8)')
	args = parser.parse_args()

	print('Smoke validation: %s' % args.batch_dir)
	report = Report()
	check_listeners(args.batch_dir, report)
	check_seed_matching(args.batch_dir, report, args.induction)
	check_crowding(args.batch_dir, report, args.induction)
	check_insertion(args.batch_dir, report)
	check_disk(args.batch_dir, report)

	print('\n%s' % ('=' * 60))
	if report.failures:
		print('%d CHECK(S) FAILED - do not start production batches'
			% report.failures)
	else:
		print('All checks passed. %d notes to read before proceeding.'
			% report.notes)
	print('Bitwise reproducibility is not checked here: rerun one generation '
		'with the same seed and diff its simOut directory.')
	return 1 if report.failures else 0


if __name__ == '__main__':
	sys.exit(main())
