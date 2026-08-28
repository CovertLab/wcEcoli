"""
Gap G1b — the copy-number cascade within single lineages, aligned on
replication-initiation events.

Every other analysis in this study averages over cells sitting at random points
in their division cycle, which smears the one thing the copy-number loop does
discretely: replication initiation. `chromosome_replication.py:171` fires
initiation when `criticalMassPerOriC = massPerOrigin / criticalInitiationMass`
crosses 1.0, so initiation is a dateable threshold crossing and `numberOfOric` is
its digital readout. Origin-proximal gene dosage then steps, rather than drifting.

This script does two things. Per lineage it draws the cascade on one clock, in the
order the chain runs:

	criticalMassPerOriC        growth -> the trigger
	numberOfOric               the discrete event
	rRNA operon copy number    replication -> dosage
	rRNA initiation events     dosage -> transcription
	total ribosomes            transcription -> pool
	instantaneous growth rate  pool -> growth, closing the loop

Across lineages it emits one row per initiation event, so the population question
becomes a within-lineage one: conditioned on an initiation, what happened next.

Why this is not the G2 comparison over again
--------------------------------------------
G2 compared steady-state endpoints between the frozen and unfrozen ladders at
matched burden and found the dosage arm doubled (2.54 -> 4.67 points) while rRNA
share, both pools and growth did not move. So the doubled signal exists and dies
between those two links. An equilibrium comparison cannot say which link breaks,
and cannot see a transient that reconverges -- the gain of a thermostat is not
visible in a comparison of two rooms at their set points.

The test here is conditional rather than marginal. In the frozen batch the
dosage-to-transcription step is roughly twice as large, so if the loop propagates,
the frozen lineages should lose rRNA initiation and ribosomes faster in the window
that follows each initiation event, even though their endpoints match. Each
lineage contributes several events, so the contrast is within-lineage and the
statistics are far better than eight rungs of endpoint.

The closure condition, stated before looking
--------------------------------------------
The loop matters if a delayed initiation lowers rRNA synthesis enough to further
delay the next initiation. That reads out as the inter-initiation interval
lengthening PROGRESSIVELY after induction, beyond the one-off lengthening the
proteome tax alone produces. `interval_min` against `event_index` is the column
pair that answers it, and the frozen-versus-unfrozen difference in that slope is
the loop's signature. Matched slopes with a doubled dosage arm would be a
substantially harder null than G2's.

What it cannot do
-----------------
Event alignment removes cycle-phase smearing and makes the contrast
within-lineage, but it does not remove the batch-to-batch floor: two batches that
should be identical already differ by one to two minutes at the unburdened
controls. Run the control variant and read every number against it.

Output
------
	<name>.csv                  one row per initiation event, all lineages
	<name>_by_variant.csv       per-variant aggregates
	<name>.pdf / .png           the summary figure
	<name>_lineages/            one cascade plot per lineage
"""

import csv
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	first_cell_with_table, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Generation at which new gene expression is switched on. Must match
# NEW_GENE_INDUCTION_GEN in models/ecoli/sim/variants/new_gene_internal_shift.py.
INDUCTION_GEN = 8

# Generations either side of induction to read. Three before gives several
# pre-induction intervals to set the baseline against.
GENS_BEFORE = 3
GENS_AFTER = 8

# Windows around an initiation event, in minutes, over which each trace is
# compared with itself. rRNA operons sit at replichore fraction ~0.20, so a fork
# leaving the origin reaches them in about 0.20 * C = 8 min; copy number and the
# transcription that follows it therefore respond inside twenty minutes. The
# pools integrate over far longer and need a window of order a generation.
FAST_WINDOW_MIN = (2.0, 20.0)
SLOW_WINDOW_MIN = (10.0, 60.0)

# Minimum rise in numberOfOric to count as an initiation. Division halves the
# origin count, so only positive steps are events; the threshold is here to
# reject any single-timestep bookkeeping flicker.
MIN_ORIC_RISE = 1

# A lineage needs this many events on each side of induction to contribute to
# the progression test, which is a comparison of slopes.
MIN_EVENTS_PER_SIDE = 2

ROWS = [
	('crit_mass', 'mass per oriC\n/ critical', 'C7'),
	('n_oric', 'origins\nper cell', 'C3'),
	('cn_rrna', 'rRNA operon\ncopy number', 'C0'),
	('rrna_init', 'rRNA initiation\nevents', 'C2'),
	('ribosome', 'total\nribosomes', 'C4'),
	('growth_rate', 'instantaneous\ngrowth rate', 'C1'),
	]

# Which window each trace is measured over, and whether an initiation is
# expected to raise or lower it.
WINDOWS = {
	'cn_rrna': FAST_WINDOW_MIN,
	'rrna_init': FAST_WINDOW_MIN,
	'ribosome': SLOW_WINDOW_MIN,
	'growth_rate': SLOW_WINDOW_MIN,
	}


def _sem(x):
	x = np.asarray(x, dtype=float)
	x = x[np.isfinite(x)]
	if x.size < 2:
		return float('nan')
	return float(np.std(x, ddof=1) / np.sqrt(x.size))


def _window_clean(t_event, window, divisions=None):
	"""
	True if this event's measurement window is safe to read.

	The window must not straddle generation 8. An event shortly before induction
	has its POST window after it, and the change is then filed as a
	pre-induction observation. Measured, that made the pre-induction ribosome
	response slide monotonically down the burden ladder -- inside a window that
	by construction precedes any expression -- so the control the analysis leans
	on was contaminated in proportion to the effect being measured.

	Division is NOT screened here. It is a real confound: per-cell counts halve
	at division, so a window spanning one mixes the two halves and dilutes the
	step being measured. That is what produced a NEGATIVE change in
	origin-proximal copy number across an initiation on three independent
	batches -- unfrozen v4, frozen v4, rich v4, all with (C + D) / tau within
	0.03 of an integer, which is exactly when initiation coincides with
	division.

	But excluding those windows is not usable. Initiation fires at a fixed cell
	mass and therefore at a fixed cycle phase, so the exclusion is
	all-or-nothing per rung, and at a 40-minute window inside a 51-minute cycle
	it removes every event at every rung but one. The confound is removed by
	DIVISION-CORRECTING the traces instead -- see _divide_correct.

	`divisions` is accepted and ignored so callers need not know which screens
	are active.
	"""
	_ = divisions
	hi = window[1]
	return (t_event - hi < 0.0) == (t_event + hi < 0.0)


def _divide_correct(t, values, divisions):
	"""
	Undo the halving of a per-cell count at each division.

	Multiplying by 2 ** (divisions so far) turns a per-cell stock into a
	per-initial-cell one, which is continuous through division. The step at an
	initiation survives -- copy number goes 2 to 4 at initiation, and at the
	following division the count halves while the factor doubles, so the
	corrected trace holds flat rather than dropping.

	Applied only to count-like traces. A ratio such as criticalMassPerOriC and a
	rate such as the instantaneous growth rate are already continuous.
	"""
	k = np.searchsorted(np.sort(divisions), t, side='right')
	return np.asarray(values, dtype=float) * (2.0 ** k)


def _window_change(t, values, t_event, window):
	"""
	Fractional change across an initiation, same trace either side of it.

	Returns nan rather than a number whenever either side is empty, so a
	truncated lineage drops the event instead of contributing a one-sided
	comparison.
	"""
	lo, hi = window
	pre = (t >= t_event - hi) & (t <= t_event - lo)
	post = (t >= t_event + lo) & (t <= t_event + hi)
	if not (pre.any() and post.any()):
		return float('nan'), float('nan'), float('nan')
	a = float(np.nanmean(values[pre]))
	b = float(np.nanmean(values[post]))
	if not np.isfinite(a) or a == 0:
		return a, b, float('nan')
	return a, b, b / a - 1.0


def _initiation_times(t, n_oric):
	"""Times at which the origin count rises, i.e. replication initiates."""
	if t.size < 2:
		return np.array([])
	rise = np.diff(n_oric)
	idx = np.where(rise >= MIN_ORIC_RISE)[0] + 1
	return t[idx]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		import pickle
		with open(simDataFile, 'rb') as handle:
			sim_data = pickle.load(handle)

		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		lineage_dir = os.path.join(plotOutDir, plotOutFileName + '_lineages')
		if not os.path.exists(lineage_dir):
			os.makedirs(lineage_dir)

		masks = self._masks(sim_data)
		events, skipped = [], 0
		for variant in variants:
			for seed in sorted(self.ap.get_seeds(variant=variant)):
				series = self._one_lineage(variant, seed, sim_data, masks)
				if series is None:
					skipped += 1
					continue
				rows = self._events_for(variant, seed, series)
				events.extend(rows)
				self._plot_lineage(lineage_dir, variant, seed, series, rows)

		if not events:
			print('No initiation events found in any lineage.')
			return

		print('\n%d initiation events across %d lineages (%d lineages skipped).'
			% (len(events), len({(e['variant'], e['seed']) for e in events}),
			skipped))
		print('  Per-cell counts are division-corrected, so a window '
			'spanning a division is usable.')
		for key in WINDOWS:
			flag = '%s_window_clean' % key
			dropped = sum(1 for e in events if not e.get(flag, 1))
			print('  %-14s %d of %d events dropped: window straddled '
				'induction.' % (key, dropped, len(events)))
		n_dirty = sum(1 for e in events
			if not e.get('interval_clean', 1) and np.isfinite(
				e.get('interval_min', float('nan'))))
		print('  %-14s %d intervals dropped: straddled induction.'
			% ('interval', n_dirty))
		self._write_events(plotOutDir, plotOutFileName, events)
		by_variant = self._aggregate(events)
		self._write_by_variant(plotOutDir, plotOutFileName, by_variant)
		self._print_summary(by_variant)
		self._plot_summary(plotOutDir, plotOutFileName, events, by_variant)
		print('\nPer-lineage cascades in %s' % lineage_dir)

	# ---- extraction ------------------------------------------------------

	@staticmethod
	def _masks(sim_data):
		"""rRNA masks at both indexings the listeners use."""
		transcription = sim_data.process.transcription
		cistron_data = transcription.cistron_data.struct_array
		rna_data = transcription.rna_data
		return dict(
			cistron_is_rrna=dict(zip(cistron_data['id'],
				cistron_data['is_rRNA'])),
			tu_is_rrna=np.asarray(rna_data['is_rRNA'], dtype=bool),
			tu_ids=list(rna_data['id']))

	def _one_lineage(self, variant, seed, sim_data, masks):
		"""Read every trace for one lineage on a clock zeroed at induction."""
		gens = list(range(max(0, INDUCTION_GEN - GENS_BEFORE),
			INDUCTION_GEN + GENS_AFTER + 1))
		cell_paths = self.ap.get_cells(variant=[variant], seed=[seed],
			generation=gens, only_successful=True)
		if len(cell_paths) == 0:
			cell_paths = self.ap.get_cells(variant=[variant], seed=[seed],
				generation=gens)
		if len(cell_paths) == 0:
			return None

		induction = self.ap.get_cells(variant=[variant], seed=[seed],
			generation=[INDUCTION_GEN])
		if len(induction) == 0:
			return None
		try:
			t0 = TableReader(os.path.join(induction[0], 'simOut',
				'Main')).readColumn('time')[0]
		except Exception as exc:  # noqa: BLE001 - a bad lineage is skippable
			print('Variant %d seed %d: no induction time (%s).'
				% (variant, seed, exc))
			return None

		synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
		rnap_cell = first_cell_with_table(cell_paths, 'RnapData')
		if synth_cell is None or rnap_cell is None:
			return None
		cistron_ids = TableReader(os.path.join(synth_cell, 'simOut',
			'RnaSynthProb')).readAttribute('cistron_ids')
		rna_ids = TableReader(os.path.join(rnap_cell, 'simOut',
			'RnapData')).readAttribute('rnaIds')

		is_rrna_c = masks['cistron_is_rrna']
		cn_idx = np.array([i for i, c in enumerate(cistron_ids)
			if is_rrna_c.get(c, False)], dtype=int)
		tu_lookup = dict(zip(masks['tu_ids'], masks['tu_is_rrna']))
		init_idx = np.array([i for i, r in enumerate(rna_ids)
			if tu_lookup.get(r, False)], dtype=int)
		if cn_idx.size == 0 or init_idx.size == 0:
			return None

		try:
			time = read_stacked_columns(cell_paths, 'Main', 'time',
				ignore_exception=True).squeeze()
			crit = read_stacked_columns(cell_paths, 'ReplicationData',
				'criticalMassPerOriC', ignore_exception=True).squeeze()
			n_oric = read_stacked_columns(cell_paths, 'ReplicationData',
				'numberOfOric', ignore_exception=True).squeeze()
			cn_rrna = read_stacked_columns(cell_paths, 'RnaSynthProb',
				'gene_copy_number', ignore_exception=True,
				fun=lambda x: x[:, cn_idx].mean(axis=1)[:, None]).squeeze()
			rrna_init = read_stacked_columns(cell_paths, 'RnapData',
				'rnaInitEvent', ignore_exception=True,
				fun=lambda x: x[:, init_idx].sum(axis=1)[:, None]).squeeze()
			growth = read_stacked_columns(cell_paths, 'Mass',
				'instantaneous_growth_rate', ignore_exception=True).squeeze()
			ribosome = self._ribosomes(cell_paths, sim_data)
		except Exception as exc:  # noqa: BLE001 - report and drop the lineage
			print('Variant %d seed %d: read failed (%s).'
				% (variant, seed, exc))
			return None
		if ribosome is None:
			return None

		series = dict(t=(np.asarray(time) - t0) / 60.0, crit_mass=crit,
			n_oric=n_oric, cn_rrna=cn_rrna, rrna_init=rrna_init,
			ribosome=ribosome, growth_rate=growth)
		n = series['t'].size
		for key, value in series.items():
			if key.startswith('_'):
				continue
			if np.asarray(value).size != n:
				print('Variant %d seed %d: %s has %d points against %d; '
					'dropping lineage.' % (variant, seed, key,
					np.asarray(value).size, n))
				return None

		# Division boundaries, on the same clock, so a measurement window can
		# be screened for spanning one. Each cell's first timestep is a
		# division except the first in the window.
		bounds = []
		for path in sorted(cell_paths):
			try:
				tt = TableReader(os.path.join(path, 'simOut',
					'Main')).readColumn('time')
			except Exception:  # noqa: BLE001 - a missing cell is skippable
				continue
			if tt.size:
				bounds.append((tt[0] - t0) / 60.0)
		series['_divisions'] = np.asarray(sorted(bounds)[1:], dtype=float)

		# Per-cell counts are division-corrected before any window is read.
		# n_oric is left alone because it drives event DETECTION, and
		# crit_mass and growth_rate are already continuous through division.
		for key in ('cn_rrna', 'rrna_init', 'ribosome'):
			series[key] = _divide_correct(series['t'], series[key],
				series['_divisions'])
		return series

	@staticmethod
	def _ribosomes(cell_paths, sim_data):
		"""Active plus inactive ribosomes per timestep."""
		umc_cell = first_cell_with_table(cell_paths, 'UniqueMoleculeCounts')
		if umc_cell is None:
			return None
		unique_ids = TableReader(os.path.join(umc_cell, 'simOut',
			'UniqueMoleculeCounts')).readAttribute('uniqueMoleculeIds')
		counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts',
			'uniqueMoleculeCounts', ignore_exception=True)
		if counts.size == 0 or 'active_ribosome' not in unique_ids:
			return None
		active = counts[:, unique_ids.index('active_ribosome')]
		# Inactive ribosomes are limited by whichever subunit is scarcer, the
		# same convention new_gene_machinery_allocation uses.
		try:
			(s30, s50) = read_stacked_bulk_molecules(cell_paths,
				([sim_data.molecule_ids.s30_full_complex],
					[sim_data.molecule_ids.s50_full_complex]),
				ignore_exception=True)
			return active + np.minimum(np.asarray(s30).squeeze(),
				np.asarray(s50).squeeze())
		except Exception:  # noqa: BLE001 - active alone is still usable
			return active

	# ---- events ----------------------------------------------------------

	@staticmethod
	def _events_for(variant, seed, series):
		"""One row per initiation event in this lineage."""
		t = series['t']
		divisions = series['_divisions']
		times = _initiation_times(t, series['n_oric'])
		if times.size == 0:
			return []

		# event_index counts from induction: -1 is the last initiation before
		# it, +1 the first after. Zero is unused, so the two sides never share
		# an index and a slope can be fitted on each separately.
		pre = times[times < 0]
		post = times[times >= 0]
		indexed = ([(t_e, -(pre.size - i)) for i, t_e in enumerate(pre)]
			+ [(t_e, i + 1) for i, t_e in enumerate(post)])

		rows = []
		for k, (t_e, index) in enumerate(indexed):
			prev = indexed[k - 1][0] if k > 0 else float('nan')
			# An interval whose start is before induction and whose end is
			# after it is neither. Including it makes a ONE-OFF step read as a
			# positive slope, which is the null this analysis has to survive:
			# tested on synthetic constant-after-step lineages it produced
			# +1.5 min per event out of nothing. Only clean intervals are
			# fitted.
			clean = int(k > 0 and np.isfinite(prev)
				and ((prev >= 0 and t_e >= 0) or (prev < 0 and t_e < 0)))
			row = dict(variant=variant, seed=seed, event_index=index,
				t_min=float(t_e),
				interval_min=float(t_e - prev) if k > 0 else float('nan'),
				interval_clean=clean,
				post_induction=int(t_e >= 0))
			at = np.argmin(np.abs(t - t_e))
			row['crit_mass_at_event'] = float(series['crit_mass'][at])
			row['n_oric_after'] = float(series['n_oric'][at])
			row['n_oric_before'] = float(
				series['n_oric'][max(0, at - 1)])
			for key, window in WINDOWS.items():
				if not _window_clean(t_e, window, divisions):
					row['%s_pre' % key] = float('nan')
					row['%s_post' % key] = float('nan')
					row['%s_change' % key] = float('nan')
					row['%s_window_clean' % key] = 0
					continue
				a, b, change = _window_change(t, series[key], t_e, window)
				row['%s_pre' % key] = a
				row['%s_post' % key] = b
				row['%s_change' % key] = change
				row['%s_window_clean' % key] = 1
			rows.append(row)
		return rows

	@staticmethod
	def _aggregate(events):
		"""Per-variant aggregates, split either side of induction."""
		out = {}
		for variant in sorted({e['variant'] for e in events}):
			sub = [e for e in events if e['variant'] == variant]
			row = dict(variant=variant, n_events=len(sub),
				n_lineages=len({e['seed'] for e in sub}))
			for side, flag in (('pre', 0), ('post', 1)):
				iv = [e['interval_min'] for e in sub
					if e['post_induction'] == flag and e['interval_clean']
					and np.isfinite(e['interval_min'])]
				row['interval_%s' % side] = (float(np.mean(iv)) if iv
					else float('nan'))
				row['interval_%s_sem' % side] = _sem(iv)
				row['n_%s' % side] = len(iv)
			row['interval_lengthening'] = (row['interval_post']
				- row['interval_pre'])

			# The progression test: interval against event index after
			# induction, fitted per lineage then averaged, so a lineage with
			# many events cannot dominate.
			slopes = []
			for seed in {e['seed'] for e in sub}:
				pts = [(e['event_index'], e['interval_min']) for e in sub
					if e['seed'] == seed and e['post_induction'] == 1
					and e['interval_clean']
					and np.isfinite(e['interval_min'])]
				if len(pts) >= MIN_EVENTS_PER_SIDE:
					x, y = zip(*pts)
					slopes.append(float(np.polyfit(x, y, 1)[0]))
			# MEDIAN, not mean. One lineage out of sixteen at +44.9 min per
			# event carried a reported mean of +3.09 where the median was
			# +0.36; the symptom was a standard error twenty times its
			# counterpart's. The median leaves genuinely large rungs alone.
			row['interval_slope_post'] = (float(np.median(slopes)) if slopes
				else float('nan'))
			row['interval_slope_post_mean'] = (float(np.mean(slopes))
				if slopes else float('nan'))
			row['interval_slope_post_sem'] = _sem(slopes)
			row['n_lineages_slope'] = len(slopes)

			for key in WINDOWS:
				for side, flag in (('pre', 0), ('post', 1)):
					vals = [e['%s_change' % key] for e in sub
						if e['post_induction'] == flag
						and np.isfinite(e['%s_change' % key])]
					row['%s_change_%s' % (key, side)] = (float(np.median(vals))
						if vals else float('nan'))
					row['%s_change_%s_sem' % (key, side)] = _sem(vals)
			out[variant] = row
		return out

	# ---- output ----------------------------------------------------------

	@staticmethod
	def _write_events(plot_out_dir, name, events):
		path = os.path.join(plot_out_dir, name + '.csv')
		fields = ['variant', 'seed', 'event_index', 'post_induction', 't_min',
			'interval_min', 'interval_clean', 'crit_mass_at_event',
			'n_oric_before', 'n_oric_after']
		for key in WINDOWS:
			fields += ['%s_pre' % key, '%s_post' % key, '%s_change' % key,
				'%s_window_clean' % key]
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in events:
				writer.writerow(row)
		print('Wrote %s' % path)

	@staticmethod
	def _write_by_variant(plot_out_dir, name, by_variant):
		path = os.path.join(plot_out_dir, name + '_by_variant.csv')
		rows = [by_variant[v] for v in sorted(by_variant)]
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
			writer.writeheader()
			for row in rows:
				writer.writerow(row)
		print('Wrote %s' % path)

	@staticmethod
	def _print_summary(by_variant):
		print('\nInter-initiation interval, minutes')
		print('%-8s %8s %8s %14s %11s %18s' % ('variant', 'pre', 'post',
			'lengthening', 'events', 'slope per event'))
		for v in sorted(by_variant):
			r = by_variant[v]
			print('%-8d %8.1f %8.1f %14s %11s %18s' % (v,
				r['interval_pre'], r['interval_post'],
				'%+.1f' % r['interval_lengthening'],
				'%d/%d' % (r['n_pre'], r['n_post']),
				'%+.2f +/- %.2f (n=%d)' % (r['interval_slope_post'],
					r['interval_slope_post_sem'], r['n_lineages_slope'])))
		print('\n  The CLOSURE CONDITION is the slope column: a one-off '
			'lengthening is the')
		print('  proteome tax, a slope that stays positive across successive '
			'events is the loop')
		print('  feeding back. Read it against the control variant, and '
			'against the other batch.')
		print('  Intervals straddling induction are EXCLUDED from both the '
			'means and the slope;')
		print('  including them makes a one-off step read as +1.5 min per '
			'event out of nothing.')

		print('\nFractional change across an initiation event, by trace')
		print('%-8s %26s %26s' % ('variant', 'before induction',
			'after induction'))
		for v in sorted(by_variant):
			r = by_variant[v]
			for key in WINDOWS:
				print('%-8s %-14s %11s %26s' % (
					'%d' % v if key == list(WINDOWS)[0] else '', key,
					'%+.2f%%' % (100 * r['%s_change_pre' % key]),
					'%+.2f%% +/- %.2f' % (100 * r['%s_change_post' % key],
						100 * r['%s_change_post_sem' % key])))
		print('\n  An initiation RAISES origin-proximal copy number, so '
			'cn_rrna_change should be')
		print('  positive. The question is whether rrna_init and ribosome '
			'follow it, and by how')
		print('  much less after induction than before.')

	@staticmethod
	def _plot_lineage(lineage_dir, variant, seed, series, rows):
		"""The poster-style cascade for one lineage."""
		t = series['t'] / 60.0
		fig, axes = plt.subplots(len(ROWS), 1, figsize=(9, 1.35 * len(ROWS)),
			sharex=True)
		events = [r['t_min'] / 60.0 for r in rows]
		for ax, (key, label, colour) in zip(np.atleast_1d(axes), ROWS):
			ax.plot(t, series[key], lw=1.0, color=colour)
			for t_e in events:
				ax.axvline(t_e, color='0.75', lw=0.6, zorder=0)
			ax.axvline(0.0, color='C3', lw=1.4)
			if key == 'crit_mass':
				# The threshold that fires initiation.
				ax.axhline(1.0, color='0.4', lw=0.7, ls='--')
			ax.set_ylabel(label, fontsize='x-small')
			ax.tick_params(labelsize='x-small')
		np.atleast_1d(axes)[0].set_title('Replication cascade — variant %d, '
			'seed %d. Red line: induction. Grey lines: initiation events.'
			% (variant, seed), fontsize='small')
		np.atleast_1d(axes)[-1].set_xlabel(
			'hours from induction (generation %d)' % INDUCTION_GEN,
			fontsize='small')
		plt.tight_layout()
		exportFigure(plt, lineage_dir,
			'cascade_var_%02d_seed_%02d' % (variant, seed), metadata=None)
		plt.close('all')

	@staticmethod
	def _plot_summary(plot_out_dir, name, events, by_variant):
		"""Interval progression and the per-event responses, across variants."""
		variants = sorted(by_variant)
		fig, axes = plt.subplots(1, 3, figsize=(13, 4))
		cmap = plt.get_cmap('viridis')
		shade = {v: cmap(0.1 + 0.75 * i / max(1, len(variants) - 1))
			for i, v in enumerate(variants)}

		# 1. interval against event index -- the closure condition
		ax = axes[0]
		for v in variants:
			pts = {}
			for e in events:
				if (e['variant'] == v and e['interval_clean']
						and np.isfinite(e['interval_min'])):
					pts.setdefault(e['event_index'], []).append(
						e['interval_min'])
			if not pts:
				continue
			x = sorted(pts)
			ax.plot(x, [np.mean(pts[i]) for i in x], marker='o', ms=3,
				lw=1.2, color=shade[v], label='v%d' % v)
		ax.axvline(0.0, color='C3', lw=1.2)
		ax.set_xlabel('initiation event, relative to induction')
		ax.set_ylabel('inter-initiation interval (min)')
		ax.set_title('Closure condition: does it keep rising?',
			fontsize='small')
		ax.legend(fontsize='xx-small', ncol=2, frameon=False)

		# 2. the slope, per variant
		ax = axes[1]
		ax.errorbar(variants,
			[by_variant[v]['interval_slope_post'] for v in variants],
			yerr=[by_variant[v]['interval_slope_post_sem'] for v in variants],
			marker='s', lw=1.2, color='C0')
		ax.axhline(0.0, color='0.5', lw=0.8, ls='--')
		ax.set_xlabel('variant')
		ax.set_ylabel('interval slope (min per event)')
		ax.set_title('Positive = loop feeding back', fontsize='small')

		# 3. per-event response of the downstream links
		ax = axes[2]
		for key, colour in (('cn_rrna', 'C0'), ('rrna_init', 'C2'),
				('ribosome', 'C4')):
			ax.errorbar(variants,
				[100 * by_variant[v]['%s_change_post' % key]
					for v in variants],
				yerr=[100 * by_variant[v]['%s_change_post_sem' % key]
					for v in variants],
				marker='o', ms=3, lw=1.2, color=colour, label=key)
		ax.axhline(0.0, color='0.5', lw=0.8, ls='--')
		ax.set_xlabel('variant')
		ax.set_ylabel('change across an initiation (%)')
		ax.set_title('Does the step propagate?', fontsize='small')
		ax.legend(fontsize='xx-small', frameon=False)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, name, metadata=None)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
