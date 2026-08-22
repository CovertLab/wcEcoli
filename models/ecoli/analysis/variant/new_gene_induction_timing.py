"""
Gap G1 — when each link in the copy-number chain moves, relative to induction.

Every other analysis in this study is a steady-state comparison between
variants, which cannot say anything about ordering. This one uses the generation
8 induction step as a natural experiment: it aligns every seed at the moment GFP
turns on, averages the traces across seeds, and reports how long each quantity
takes to travel half of its total post-induction change.

Traces, all normalised to their own pre-induction mean so they share one axis:

	new gene monomer count            the driver
	mean rRNA operon copy number      the dosage arm's input, origin-proximal
	a terminus-proximal control gene  dosage control, falls less than rRNA
	total ribosome pool               the machinery, hit two ways
	total RNAP pool
	instantaneous growth rate

The terminus gene is a dosage control, not a null control: it loses copies too,
just fewer. At replichore fraction 0.95 Cooper-Helmstetter gives 2**(22/tau), so
over the minimal ladder's 52 -> 89 min it should fall about 11%, against about
25% for the rRNA operons at fraction 0.20. A terminus trace that rises, or that
falls as far as rRNA does, means something upstream is wrong.

What this can and cannot establish
----------------------------------
It cannot establish causation between the copy number and the pools. Both
respond to the same step, and the pools are additionally drawn down by GFP
translation directly, so whichever moves first would move first regardless of
whether gene dosage feeds back at all. A lag is expected on purely mechanical
grounds too: copy number is a near-instantaneous function of the current
replication state, whereas the pools integrate over the preceding generation.

What it can do is separate the immediate hit from the delayed one. If the pool
loss is largely complete before the copy number has moved, then the copy-number
arm cannot be carrying the early transient, and only the later, sustained part
of the pool loss is available for it to explain. That is a decomposition in
time, and it bounds the loop's contribution to the early response at roughly
zero without needing a counterfactual simulation.

It is the companion to the frozen-expectation batch, which gives the magnitude.
This gives the timing. Neither alone closes the question.

Run it on the highest-burden variant and on the control. On the control every
trace should be flat, which is the negative that makes the rest readable.
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

# Generation at which new gene expression is switched on. Must match
# NEW_GENE_INDUCTION_GEN in models/ecoli/sim/variants/new_gene_internal_shift.py.
INDUCTION_GEN = 8

# Generations either side of induction to keep. Two before is enough for a
# baseline and eight after covers the approach to the new steady state; the
# batches run 24 generations, so this is well inside the data.
GENS_BEFORE = 2
GENS_AFTER = 8

# Common time grid, in minutes relative to the start of the induction
# generation. The grid is fixed so that batches remain comparable; how much of
# it is usable is decided from the data, by COVERAGE_FRAC below.
#
# Cells sit at different points of their division cycle at any given moment, so
# binning across seeds averages over cycle phase rather than resolving it. Copy
# number oscillates strongly within a cycle, so that averaging is not cosmetic:
# a bin holding only two or three seeds samples a narrow slice of cycle phase
# instead of the mean over it, and can come out on either side of the truth.
#
# The grid is piecewise. Coverage is complete in the first hour and the
# half-times of interest all fall there, so bins are fine early; later bins only
# need to establish a plateau, so they are coarse.
BIN_MINUTES_EARLY = 2.0
EARLY_UNTIL_MIN = 60.0
BIN_MINUTES = 10.0
WINDOW_BEFORE_MIN = 150.0
WINDOW_AFTER_MIN = 900.0

# A post-induction bin is usable only if this fraction of the seeds reached it.
#
# Without this the analysis fails silently and in the worst possible direction.
# Eight generations after induction spans about 420 min in minimal media and 220
# in rich, so a fixed 900-minute window puts its tail beyond where any seed
# reaches: the tail then holds only the seeds that happened to divide slowest,
# which are the most burden-affected ones, while the output still reports all 16
# seeds. Measured on the P4 minimal batch, that produced a terminus-gene change
# of +20.9% where Cooper-Helmstetter predicts -9.6%, and an rRNA change of -9.4%
# against a predicted -24.6%. Both signs and magnitudes were set by which seeds
# survived into the tail.
#
# Note this must be counted in SEEDS, not in samples. Each seed contributes
# hundreds of timesteps to a ten-minute bin, so a threshold on the sample count
# would be met by any single seed and would mask nothing at all.
COVERAGE_FRAC = 0.8

# Fraction of the covered post-induction span used as the plateau, measured in
# time rather than in bins because the grid is piecewise.
TAIL_FRAC = 0.25

# Terminus-proximal reference, resolved by position at run time. Same convention
# as models/ecoli/analysis/multigen/copy_number_lineage_trace.py.
TERMINUS_TARGET_FRACTION = 0.95

# Fraction of the total post-induction change used as the timing statistic.
HALF = 0.5

TRACES = [
	('new_gene_monomer', 'new gene\nmonomer', 'C3'),
	('cn_rrna', 'rRNA operon\ncopy number', 'C0'),
	('cn_terminus', 'terminus gene\ncopy number', '0.55'),
	('ribosome_pool', 'total\nribosomes', 'C2'),
	('rnap_pool', 'total\nRNAP', 'C4'),
	('growth_rate', 'instantaneous\ngrowth rate', 'C1'),
	]


def _bin_edges():
	"""Return the common time grid edges and centres, in minutes.

	Coarse before induction, fine for the first EARLY_UNTIL_MIN afterwards,
	coarse again beyond it.
	"""
	edges = np.concatenate([
		np.arange(-WINDOW_BEFORE_MIN, 0.0, BIN_MINUTES),
		np.arange(0.0, EARLY_UNTIL_MIN, BIN_MINUTES_EARLY),
		np.arange(EARLY_UNTIL_MIN, WINDOW_AFTER_MIN + BIN_MINUTES,
			BIN_MINUTES)])
	return edges, 0.5 * (edges[:-1] + edges[1:])


def _plateau(centres, values, usable):
	"""
	Return the post-induction plateau value and the last covered time.

	The plateau is the mean over the last TAIL_FRAC of the *covered* span, so it
	moves with wherever the data actually ends rather than sitting at a fixed
	offset that may be past the end of every lineage.
	"""
	if not np.any(usable):
		return float('nan'), float('nan')
	t_max = float(centres[usable].max())
	tail = usable & (centres >= (1.0 - TAIL_FRAC) * t_max)
	if not np.any(tail):
		return float('nan'), t_max
	return float(np.nanmean(values[tail])), t_max


def _t_half(centres, values, baseline, final, usable):
	"""
	Return the time at which a trace first covers HALF of its total change.

	Searches only covered bins, and takes the plateau from _plateau rather than
	recomputing it, so the timing and the magnitude cannot disagree about which
	part of the window they trust.
	"""
	if not (np.isfinite(baseline) and np.isfinite(final)) or baseline == 0:
		return float('nan')
	total = final - baseline
	if abs(total / baseline) < 0.02:
		# Under a 2% move there is no step to time, only noise.
		return float('nan')
	target = baseline + HALF * total
	post, post_t = values[usable], centres[usable]
	if total < 0:
		reached = np.where(post <= target)[0]
	else:
		reached = np.where(post >= target)[0]
	if reached.size == 0:
		return float('nan')
	return float(post_t[reached[0]])


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		self._warned_unsuccessful = False
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		variants = sorted(self.ap.get_variants())
		if not variants:
			print('No variants found.')
			return

		edges, centres = _bin_edges()
		rows = []
		curves = {}
		for variant in variants:
			result = self._one_variant(variant, sim_data, edges, centres)
			if result is None:
				continue
			curves[variant] = result
			for key, _, _ in TRACES:
				rows.append(dict(variant=variant, trace=key,
					baseline=result['baseline'][key],
					final=result['final'][key],
					change=result['change'][key],
					t_half_min=result['t_half'][key],
					t_max_covered_min=result['t_max'][key],
					n_seeds=result['n_seeds']))

		if not curves:
			print('No variant had a readable induction window.')
			return

		self._print_summary(curves)
		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._plot(plotOutDir, plotOutFileName, centres, curves)

	# ---- extraction ------------------------------------------------------

	def _one_variant(self, variant, sim_data, edges, centres):
		"""Bin every trace for one variant onto the common grid."""
		seeds = sorted(self.ap.get_seeds(variant=variant))
		if not seeds:
			return None

		cistron_data = sim_data.process.transcription.cistron_data.struct_array
		monomer_data = sim_data.process.translation.monomer_data.struct_array
		replication = sim_data.process.replication

		new_cistrons = cistron_data[cistron_data['is_new_gene']]['id'].tolist()
		cistron_to_monomer = dict(zip(
			monomer_data['cistron_id'], monomer_data['id']))
		new_monomers = [cistron_to_monomer[c] for c in new_cistrons
			if c in cistron_to_monomer]

		# Accumulate sum and count per bin so seeds of unequal length combine.
		# `covered` is separate and counts SEEDS, not samples: it is incremented
		# once per seed per bin the seed reaches at all. A threshold on `counts`
		# would not work, because one seed puts hundreds of timesteps into a
		# ten-minute bin.
		totals = {key: np.zeros(centres.size) for key, _, _ in TRACES}
		counts = {key: np.zeros(centres.size) for key, _, _ in TRACES}
		covered = np.zeros(centres.size)
		n_seeds = 0

		for seed in seeds:
			series = self._one_seed(variant, seed, cistron_data, replication,
				new_monomers, sim_data)
			if series is None:
				continue
			n_seeds += 1
			t = series.pop('_time_min')
			idx = np.digitize(t, edges) - 1
			keep = (idx >= 0) & (idx < centres.size)
			idx = idx[keep]
			np.add.at(covered, np.unique(idx), 1.0)
			for key, values in series.items():
				v = np.asarray(values, dtype=float)[keep]
				good = np.isfinite(v)
				np.add.at(totals[key], idx[good], v[good])
				np.add.at(counts[key], idx[good], 1.0)

		if n_seeds == 0:
			print('Variant %d: no seed had a readable induction window.'
				% variant)
			return None

		mean = {}
		for key, _, _ in TRACES:
			with np.errstate(invalid='ignore', divide='ignore'):
				mean[key] = np.where(counts[key] > 0,
					totals[key] / np.maximum(counts[key], 1), np.nan)

		# The usable window is decided by seed coverage, not by the grid.
		usable = (centres > 0) & (covered >= COVERAGE_FRAC * n_seeds)
		pre = centres < 0
		baseline, final, change, t_half, t_max, normed = {}, {}, {}, {}, {}, {}
		for key, _, _ in TRACES:
			b = float(np.nanmean(mean[key][pre])) if np.any(
				np.isfinite(mean[key][pre])) else float('nan')
			baseline[key] = b
			f, tm = _plateau(centres, mean[key], usable)
			final[key] = f
			t_max[key] = tm
			change[key] = f / b - 1.0 if b else float('nan')
			t_half[key] = _t_half(centres, mean[key], b, f, usable)
			normed[key] = mean[key] / b if b else mean[key] * np.nan

		return dict(normed=normed, baseline=baseline, final=final,
			change=change, t_half=t_half, t_max=t_max, n_seeds=n_seeds,
			usable=usable, covered=covered)

	def _one_seed(self, variant, seed, cistron_data, replication,
			new_monomers, sim_data):
		"""Read one seed's traces on a time axis measured from induction."""
		first_gen = max(0, INDUCTION_GEN - GENS_BEFORE)
		gens = list(range(first_gen, INDUCTION_GEN + GENS_AFTER + 1))
		cell_paths = self._cells(variant, seed, gens)
		if len(cell_paths) == 0:
			return None

		# The induction instant is the start of the induction generation for
		# this seed. Without it there is nothing to align to, so skip the seed.
		induction_paths = self._cells(variant, seed, [INDUCTION_GEN])
		if len(induction_paths) == 0:
			return None
		try:
			t0 = TableReader(os.path.join(induction_paths[0], 'simOut',
				'Main')).readColumn('time')[0]
		except Exception as exc:  # noqa: BLE001 - a bad seed is skippable
			print('Variant %d seed %d: no readable induction time (%s).'
				% (variant, seed, exc))
			return None

		readable = first_cell_with_table(cell_paths, 'RnaSynthProb')
		if readable is None:
			return None
		cistron_ids = TableReader(os.path.join(readable, 'simOut',
			'RnaSynthProb')).readAttribute('cistron_ids')

		rrna_idx, term_idx = self._gene_indices(cistron_ids, cistron_data,
			replication)
		if rrna_idx.size == 0 or term_idx is None:
			return None

		time = read_stacked_columns(cell_paths, 'Main', 'time',
			ignore_exception=True).squeeze()
		if np.asarray(time).size == 0:
			return None

		cn_rrna = read_stacked_columns(cell_paths, 'RnaSynthProb',
			'gene_copy_number', ignore_exception=True,
			fun=lambda x: x[:, rrna_idx].mean(axis=1)[:, None]).squeeze()
		cn_term = read_stacked_columns(cell_paths, 'RnaSynthProb',
			'gene_copy_number', ignore_exception=True,
			fun=lambda x: x[:, [term_idx]]).squeeze()
		growth = read_stacked_columns(cell_paths, 'Mass',
			'instantaneous_growth_rate', ignore_exception=True).squeeze()

		pools = self._pools(cell_paths, sim_data)
		if pools is None:
			return None
		ribosome, rnap = pools

		monomer = self._new_gene_monomer(cell_paths, new_monomers)
		if monomer is None:
			monomer = np.full(np.asarray(time).shape, np.nan)

		series = dict(_time_min=(np.asarray(time) - t0) / 60.0,
			new_gene_monomer=monomer, cn_rrna=cn_rrna, cn_terminus=cn_term,
			ribosome_pool=ribosome, rnap_pool=rnap, growth_rate=growth)

		n = series['_time_min'].size
		for key, value in series.items():
			if np.asarray(value).size != n:
				print('Variant %d seed %d: %s has %d points against %d; '
					'skipping seed.' % (variant, seed, key,
					np.asarray(value).size, n))
				return None
		return series

	def _cells(self, variant, seed, gens):
		"""
		Return cell paths, preferring successful cells.

		Dying lineages matter here in a way they do not for a steady-state
		average: a lineage whose pools collapse would be read as a fast
		response. No sibling analysis in this study sets only_successful
		though, so fall back rather than silently return nothing if the flag
		is not populated in this batch.
		"""
		paths = self.ap.get_cells(variant=[variant], seed=[seed],
			generation=gens, only_successful=True)
		if len(paths):
			return paths
		if not self._warned_unsuccessful:
			print('No cells are flagged successful; including all cells. '
				'Check for dying lineages before reading the timings.')
			self._warned_unsuccessful = True
		return self.ap.get_cells(variant=[variant], seed=[seed],
			generation=gens)

	@staticmethod
	def _gene_indices(cistron_ids, cistron_data, replication):
		"""Return the rRNA cistron indices and one terminus-proximal index."""
		is_rrna_of = dict(zip(cistron_data['id'], cistron_data['is_rRNA']))
		coord_of = dict(zip(cistron_data['id'],
			cistron_data['replication_coordinate']))
		rrna_idx = np.array([i for i, c in enumerate(cistron_ids)
			if is_rrna_of.get(c, False)], dtype=int)

		coords = np.array([coord_of.get(c, 0) for c in cistron_ids],
			dtype=float)
		arm = np.where(coords >= 0, replication.replichore_lengths[0],
			replication.replichore_lengths[1])
		frac = np.abs(coords) / arm
		term_idx = int(np.argmin(np.abs(frac - TERMINUS_TARGET_FRACTION)))
		return rrna_idx, term_idx

	@staticmethod
	def _pools(cell_paths, sim_data):
		"""Return total ribosome and total RNAP counts per timestep."""
		umc_cell = first_cell_with_table(cell_paths, 'UniqueMoleculeCounts')
		if umc_cell is None:
			return None
		unique_ids = TableReader(os.path.join(umc_cell, 'simOut',
			'UniqueMoleculeCounts')).readAttribute('uniqueMoleculeIds')
		counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts',
			'uniqueMoleculeCounts', ignore_exception=True)
		if counts.size == 0:
			return None
		if 'active_ribosome' not in unique_ids or 'active_RNAP' not in unique_ids:
			return None
		active_ribosome = counts[:, unique_ids.index('active_ribosome')]
		active_rnap = counts[:, unique_ids.index('active_RNAP')]

		# Inactive ribosomes are limited by whichever subunit is scarcer, the
		# same convention new_gene_machinery_allocation uses.
		try:
			(rnap_free, s30, s50) = read_stacked_bulk_molecules(cell_paths,
				([sim_data.molecule_ids.full_RNAP],
					[sim_data.molecule_ids.s30_full_complex],
					[sim_data.molecule_ids.s50_full_complex]),
				ignore_exception=True)
		except Exception as exc:  # noqa: BLE001 - fall back to active only
			print('Could not read free pools (%s); using active counts.' % exc)
			return active_ribosome, active_rnap
		return (active_ribosome + np.minimum(s30, s50),
			active_rnap + np.asarray(rnap_free).squeeze())

	@staticmethod
	def _new_gene_monomer(cell_paths, new_monomers):
		"""Return the summed new gene monomer count per timestep."""
		if not new_monomers:
			return None
		try:
			read = read_stacked_bulk_molecules(cell_paths, (new_monomers,),
				ignore_exception=True)
		except Exception as exc:  # noqa: BLE001 - driver row is optional
			print('Could not read the new gene monomer (%s).' % exc)
			return None
		counts = np.asarray(read[0], dtype=float)
		return counts.sum(axis=1) if counts.ndim > 1 else counts

	# ---- output ----------------------------------------------------------

	@staticmethod
	def _print_summary(curves):
		print()
		print('G1 — covered window per variant. The grid runs to %.0f min, but '
			'only bins' % WINDOW_AFTER_MIN)
		print('reached by at least %.0f%% of seeds are used. A short covered '
			'window is not' % (COVERAGE_FRAC * 100))
		print('a failure — it is how long the lineages actually ran.')
		print()
		print('%-8s %7s %17s %10s' % ('variant', 'seeds', 'covered to (min)',
			'bins used'))
		for variant in sorted(curves):
			d = curves[variant]
			tm = [v for v in d['t_max'].values() if np.isfinite(v)]
			print('%-8d %7d %17s %10d' % (variant, d['n_seeds'],
				'%.0f' % max(tm) if tm else 'none',
				int(d['usable'].sum())))

		print()
		print('G1 — time to half of the post-induction change, in minutes.')
		print('nan means the trace moved less than 2%, i.e. there was no step '
			'to time.')
		print()
		header = '%-8s %7s' % ('variant', 'seeds')
		for key, _, _ in TRACES:
			header += ' %18s' % key
		print(header)
		for variant in sorted(curves):
			d = curves[variant]
			line = '%-8d %7d' % (variant, d['n_seeds'])
			for key, _, _ in TRACES:
				line += ' %18s' % ('%.0f (%+.1f%%)' % (d['t_half'][key],
					d['change'][key] * 100)
					if np.isfinite(d['t_half'][key])
					else 'flat (%+.1f%%)' % (d['change'][key] * 100))
			print(line)

		# The one comparison this analysis exists to make.
		print()
		for variant in sorted(curves):
			d = curves[variant]
			t_pool = d['t_half']['ribosome_pool']
			t_cn = d['t_half']['cn_rrna']
			if not (np.isfinite(t_pool) and np.isfinite(t_cn)):
				continue
			lead = t_cn - t_pool
			print('Variant %d: the ribosome pool reaches half its change at '
				'%.0f min, rRNA operon copy number at %.0f min — the pool '
				'leads by %.0f min (%.2f generations at the pre-induction '
				'doubling time).' % (variant, t_pool, t_cn, lead,
				lead / 60.0))
			if lead > 0:
				print('  The pool moves first, so the early transient cannot '
					'be dosage-mediated; only the part of the pool loss after '
					'%.0f min is available for the copy-number arm to '
					'explain.' % t_cn)
			else:
				print('  Copy number moves first, which is the ordering the '
					'loop requires. It does not by itself establish '
					'causation — see this module docstring.')

	@staticmethod
	def _write_csv(plotOutDir, plotOutFileName, rows):
		path = os.path.join(plotOutDir, plotOutFileName + '.csv')
		fields = ['variant', 'trace', 'n_seeds', 'baseline', 'final', 'change',
			't_half_min', 't_max_covered_min']
		with open(path, 'w', newline='') as f:
			writer = csv.DictWriter(f, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow(row)
		print('\nWrote %s' % path)

	@staticmethod
	def _plot(plotOutDir, plotOutFileName, centres, curves):
		variants = sorted(curves)
		fig, axes = plt.subplots(len(TRACES), 1, figsize=(8.5,
			1.5 * len(TRACES)), sharex=True)
		cmap = plt.get_cmap('viridis')
		# Where every variant has run out of covered bins.
		all_tmax = [v for d in curves.values() for v in d['t_max'].values()
			if np.isfinite(v)]
		t_uncovered = max(all_tmax) if all_tmax else float('nan')
		shades = {v: cmap(0.1 + 0.75 * i / max(1, len(variants) - 1))
			for i, v in enumerate(variants)}

		for ax, (key, label, _) in zip(np.atleast_1d(axes), TRACES):
			for variant in variants:
				d = curves[variant]
				ax.plot(centres / 60.0, d['normed'][key], lw=1.4,
					color=shades[variant], label='variant %d' % variant)
				t = d['t_half'][key]
				if np.isfinite(t):
					ax.axvline(t / 60.0, color=shades[variant], lw=0.9,
						ls=':', alpha=0.8)
			ax.axhline(1.0, color='0.7', lw=0.7, ls='--')
			ax.axvline(0.0, color='C3', lw=1.0)
			# Grey out where too few seeds reached, so a short window reads as
			# a short window rather than as a plateau.
			if np.isfinite(t_uncovered):
				ax.axvspan(t_uncovered / 60.0, centres[-1] / 60.0,
					color='0.92', zorder=0)
			ax.set_ylabel(label, fontsize='x-small')
			ax.tick_params(labelsize='x-small')

		np.atleast_1d(axes)[0].legend(fontsize='xx-small', ncol=3,
			frameon=False)
		np.atleast_1d(axes)[0].set_title('Timing of the copy-number chain '
			'across induction — dotted lines mark half of each total change',
			fontsize='small')
		np.atleast_1d(axes)[-1].set_xlabel(
			'hours from induction (generation %d)' % INDUCTION_GEN,
			fontsize='small')
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata=None)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
