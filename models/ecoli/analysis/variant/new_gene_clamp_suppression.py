"""
Measure how much transcription the initiation cap removes from the construct.

Why this exists
---------------
`new_gene_promoter_crowding` answers whether a transcription unit is pinned at
the cap. It cannot answer how hard it is pushing against it, because every
quantity it reads is post-clamp. On the P1 minimal batch that left the central
number unmeasurable: the construct is flagged overcrowded in 85% of timesteps at
variant 1, running at 96% of the physical ceiling, and from that alone the
pre-clamp demand could be 1.01x the cap or 5x it.

The distinction matters three ways:

  1. Whether the ~1.3x higher realised per-copy initiation at the
     origin-proximal position is a real positional effect or a truncation. The
     pin equalises `factor` to six significant figures, so per-copy rates should
     be matched or slightly lower at the position with more copies.
  2. Whether the construct's dosage share (48.9% at P1 against 33.2% at P4) is a
     measurement or an upper bound. Clipping the per-copy rate at the control
     rung and not at the burden rung flattens the per-copy channel and inflates
     the dosage attribution.
  3. How much of the transcription-only control's growth cost (variant 1, where
     the construct makes mRNA and no protein) is being suppressed away.

What it measures
----------------
`transcript_initiation.py` writes both sides of the clamp to the same listener,
so all of this is already on disk and needs no rerun:

  target_rna_synth_prob   pre-clamp, summed over the TU's promoters
  actual_rna_synth_prob   post-clamp, same summation
  max_p                   the per-promoter cap, scalar per timestep
  promoter_copy_number    promoters per TU
  tu_is_overcrowded       whether the clamp bound

Per variant, per target class, reduced per timestep and then averaged:

  demand_per_promoter    target / promoter_copy_number -- assigned probability
                         per promoter, the quantity the pin governs
  actual_per_promoter    actual / promoter_copy_number -- what survives
  max_p                  the cap those two are compared against
  demand_over_cap        demand_per_promoter / max_p; > 1 means this TU's own
                         assigned probability exceeds what the footprint allows
  frac_demand_over_cap   fraction of transcribing timesteps with
                         demand_per_promoter > max_p
  suppression            1 - actual/target. Positive means the clamp took
                         probability away. NEGATIVE means the TU was a net
                         beneficiary: when other promoters are capped,
                         `scale_the_rest_by` rescales everyone else upward, so an
                         uncapped TU gains
  demand_over_cap_when_clamped
                         demand_over_cap restricted to flagged timesteps

Two different things called "overcrowded"
----------------------------------------
`tu_is_overcrowded` is not the same as "this TU demanded more than the cap".
The clamp loop caps the offenders, rescales everyone else *upward* to
re-normalise, and then re-tests -- `is_overcrowded |= (probs > max_p)`. A
promoter whose own demand was comfortably under the cap can therefore be lifted
over it by other promoters being capped, and is flagged. So:

  overcrowded_frac high, frac_demand_over_cap ~0
        -> collateral. The TU is not pushing against the physical limit; it is
           being pushed into it by the rest of the transcriptome. Its assigned
           rate is intact and the flag says little about this TU.
  overcrowded_frac high, frac_demand_over_cap high
        -> self-limited. This TU genuinely demands more than the footprint
           allows, and its realised rate is set by the cap rather than by the
           pin.

`new_gene_promoter_crowding` reports only the first of those two columns, so a
high overcrowded fraction there is ambiguous between the two cases. This script
exists to resolve that ambiguity.

Reading the result
------------------
  suppression ~0                  -> the cap is decorative. Realised rates are
                                     assigned rates and every cross-position
                                     comparison stands as measured.
  suppression small and positive  -> grazing. The clamp binds but removes
                                     little; correct the dosage share by the
                                     suppression and move on.
  suppression large               -> the realised per-copy rate is set by the
                                     footprint, not by the pin. Cross-position
                                     per-copy comparisons are censored and the
                                     dosage share is an upper bound only.
  suppression negative            -> this TU is gaining from other promoters
                                     being capped, which is worth knowing before
                                     attributing its resilience to its biology.

Ratios are formed per timestep and then averaged, never as a ratio of means:
copy number changes within a cell cycle, so mean(target)/mean(copies) is a
different quantity from mean(target/copies).

Index-matching note: every column read here is TU-indexed against
RnaSynthProb's own `rnaIds` except `max_p`, which is scalar. Nothing is matched
by position across tables, and because all five columns live in one table they
cover the same cells by construction -- the length check below is belt and
braces.

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

# Probabilities below this are treated as "not transcribing" rather than divided
# through, which would otherwise produce enormous meaningless ratios.
MIN_PROB = 1e-12


def _window(n_generation):
	if n_generation <= IGNORE_FIRST_N_GENS:
		return None
	return np.arange(IGNORE_FIRST_N_GENS, n_generation)


def _sem(values):
	"""Return the standard error of the mean. NaN if fewer than two values."""
	v = np.asarray(values, dtype=float).ravel()
	v = v[np.isfinite(v)]
	if v.size < 2:
		return float('nan')
	return float(np.std(v, ddof=1) / np.sqrt(v.size))


def _summed_over(idx):
	"""Per cell: sum the given subcolumns, keeping one row per timestep.

	Time is deliberately preserved. Reading these columns unreduced would stack
	(timesteps x 3266) per cell and run to gigabytes; collapsing time instead
	would destroy the per-timestep ratios this analysis is built on.
	"""
	def fn(x):
		return np.asarray(x)[:, idx].sum(axis=1, keepdims=True)
	return fn


def _nanmean(x):
	v = np.asarray(x, dtype=float).ravel()
	v = v[np.isfinite(v)]
	return float(np.mean(v)) if v.size else float('nan')


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

		rows = []
		for variant in variants:
			measured = self._measure(variant, generations, sim_data)
			if measured is None:
				print('No usable cells for variant %d; skipping.' % variant)
				continue
			rows.extend(measured)

		if not rows:
			print('Nothing measurable.')
			return

		self._write_csv(plotOutDir, plotOutFileName, rows)
		self._report(rows)
		self._plot(plotOutDir, plotOutFileName, rows, metadata)

	# ------------------------------------------------------------ targets ----

	def _target_indices(self, sim_data, synth_ids):
		"""Return {class name: TU indices} for the construct and rRNA."""
		transcription = sim_data.process.transcription
		rna_data = transcription.rna_data
		cistron_data = transcription.cistron_data.struct_array

		new_cistrons = sorted(
			set(cistron_data[cistron_data['is_new_gene']]['id'].tolist()))
		construct_ids = []
		if new_cistrons:
			matrix = transcription.cistron_tu_mapping_matrix.toarray()
			cistron_ids = list(cistron_data['id'])
			cistron_rows = [cistron_ids.index(c) for c in new_cistrons]
			tu_mask = matrix[cistron_rows, :].sum(axis=0) > 0
			construct_ids = [str(i) for i in
				rna_data.struct_array['id'][tu_mask]]

		rrna_ids = [str(i) for i in
			rna_data.struct_array['id'][rna_data['is_rRNA']]]

		targets = {}
		for name, ids in (('construct', construct_ids), ('rrna', rrna_ids)):
			idx = np.array([synth_ids.index(t) for t in ids if t in synth_ids])
			if idx.size:
				targets[name] = idx
			else:
				print('  %s: no transcription units present in RnaSynthProb.'
					% name)
		return targets

	# ------------------------------------------------------------ measure ----

	def _measure(self, variant, generations, sim_data):
		cell_paths = self.ap.get_cells(
			variant=[variant], generation=generations)
		if len(cell_paths) == 0:
			return None

		# Every column below comes out of RnaSynthProb, so restricting to cells
		# with a readable RnaSynthProb up front guarantees the five reads cover
		# the same cells. read_stacked_columns skips an unreadable cell rather
		# than failing, so without this a gappy cell could drop out of one read
		# and not another and the per-timestep ratios would misalign.
		cell_paths = np.array([p for p in cell_paths
			if os.path.exists(os.path.join(
				p, 'simOut', 'RnaSynthProb', 'attributes.json'))])
		if len(cell_paths) == 0:
			print('  variant %d: no cell with a readable RnaSynthProb.'
				% variant)
			return None

		synth_cell = first_cell_with_table(cell_paths, 'RnaSynthProb')
		if synth_cell is None:
			return None
		synth_ids = TableReader(os.path.join(
			synth_cell, 'simOut', 'RnaSynthProb')).readAttribute('rnaIds')

		targets = self._target_indices(sim_data, synth_ids)
		if not targets:
			return None

		max_p = read_stacked_columns(cell_paths, 'RnaSynthProb', 'max_p',
			ignore_exception=True)
		if max_p.size == 0:
			print('  variant %d: max_p is empty.' % variant)
			return None
		max_p = max_p.ravel()

		rows = []
		for name, idx in sorted(targets.items()):
			row = self._measure_one(variant, name, idx, cell_paths, max_p)
			if row is not None:
				rows.append(row)
		return rows

	def _measure_one(self, variant, name, idx, cell_paths, max_p):
		reader = lambda column: read_stacked_columns(cell_paths,
			'RnaSynthProb', column, ignore_exception=True,
			fun=_summed_over(idx)).ravel()
		target = reader('target_rna_synth_prob')
		actual = reader('actual_rna_synth_prob')
		copies = reader('promoter_copy_number')
		crowded = reader('tu_is_overcrowded')

		lengths = {'target': target.size, 'actual': actual.size,
			'copies': copies.size, 'overcrowded': crowded.size,
			'max_p': max_p.size}
		if len(set(lengths.values())) != 1 or target.size == 0:
			print('  variant %d, %s: column lengths disagree (%s); skipping.'
				% (variant, name, lengths))
			return None

		with np.errstate(divide='ignore', invalid='ignore'):
			demand = np.where(copies > 0, target / copies, np.nan)
			realised = np.where(copies > 0, actual / copies, np.nan)
			suppression = np.where(target > MIN_PROB, 1.0 - actual / target,
				np.nan)
			over_cap = np.where(max_p > 0, demand / max_p, np.nan)

		# The construct is silent at the knockout variant, where every ratio is
		# an undefined 0/0. Report it as absent rather than as zero suppression.
		expressing = np.isfinite(demand) & (demand > MIN_PROB)
		clamped = expressing & (crowded > 0.5)

		return dict(
			variant=variant,
			target_class=name,
			n_cells=int(len(cell_paths)),
			n_timesteps=int(target.size),
			frac_expressing=float(np.mean(expressing)),
			overcrowded_frac=float(np.mean(crowded > 0.5)),
			# How often this TU's OWN demand exceeds the cap, as distinct from
			# how often it ends up flagged after everyone else is rescaled.
			frac_demand_over_cap=(float(np.mean(over_cap[expressing] > 1.0))
				if expressing.any() else float('nan')),
			max_p=float(np.mean(max_p)),
			demand_per_promoter=_nanmean(demand[expressing]),
			demand_per_promoter_sem=_sem(demand[expressing]),
			actual_per_promoter=_nanmean(realised[expressing]),
			actual_per_promoter_sem=_sem(realised[expressing]),
			demand_over_cap=_nanmean(over_cap[expressing]),
			demand_over_cap_when_clamped=_nanmean(over_cap[clamped]),
			suppression=_nanmean(suppression[expressing]),
			suppression_sem=_sem(suppression[expressing]),
			suppression_when_clamped=_nanmean(suppression[clamped]),
			)

	# -------------------------------------------------------------- output ----

	FIELDS = ['variant', 'target_class', 'n_cells', 'n_timesteps',
		'frac_expressing', 'overcrowded_frac', 'frac_demand_over_cap', 'max_p',
		'demand_per_promoter', 'demand_per_promoter_sem',
		'actual_per_promoter', 'actual_per_promoter_sem', 'demand_over_cap',
		'demand_over_cap_when_clamped', 'suppression', 'suppression_sem',
		'suppression_when_clamped']

	def _write_csv(self, plot_out_dir, plot_out_file_name, rows):
		path = os.path.join(plot_out_dir, plot_out_file_name + '.csv')
		with open(path, 'w') as handle:
			w = csv.DictWriter(handle, fieldnames=self.FIELDS)
			w.writeheader()
			for r in rows:
				w.writerow({k: r.get(k, '') for k in self.FIELDS})
		print('\nWrote %s' % path)

	def _report(self, rows):
		for name in sorted({r['target_class'] for r in rows}):
			subset = [r for r in rows if r['target_class'] == name]
			print('\nClamp suppression: %s' % name)
			print('  %-8s %-10s %-11s %-11s %-11s %-11s %-11s %s'
				% ('variant', 'flagged', 'own>cap', 'demand/pr', 'actual/pr',
					'max_p', 'demand/cap', 'suppression'))
			for r in subset:
				print('  %-8d %-10.4f %-11.4f %-11.5f %-11.5f %-11.5f '
					'%-11.4f %+.5f'
					% (r['variant'], r['overcrowded_frac'],
						r['frac_demand_over_cap'], r['demand_per_promoter'],
						r['actual_per_promoter'], r['max_p'],
						r['demand_over_cap'], r['suppression']))

			expressing = [r for r in subset if r['frac_expressing'] > 0.01]
			if not expressing:
				print('  Never transcribes; no verdict.')
				continue
			supp = max(r['suppression'] for r in expressing
				if np.isfinite(r['suppression']))
			flagged = max(r['overcrowded_frac'] for r in expressing)
			own = max(r['frac_demand_over_cap'] for r in expressing
				if np.isfinite(r['frac_demand_over_cap']))
			print('  peak flagged %.4f, peak own-demand-over-cap %.4f, '
				'peak suppression %+.5f' % (flagged, own, supp))

			if flagged > 0.05 and own < 0.05:
				print('  NOTE: %s is flagged often but rarely demands more '
					'than the cap itself. It is being lifted over the cap by '
					'other promoters being capped, not pushing against the '
					'footprint. Do not read its flagged fraction as '
					'saturation.' % name)
			if supp < -0.005:
				print('  NOTE: suppression is NEGATIVE, so %s is a net '
					'beneficiary of the clamp -- it gains probability when '
					'other promoters are capped and rescaling lifts it.'
					% name)

			if abs(supp) < 0.01:
				print('  VERDICT: the cap is decorative for %s. Realised '
					'rates are assigned rates, and cross-position per-copy '
					'comparisons stand as measured.' % name)
			elif supp < 0.10:
				print('  VERDICT: %s is GRAZING the cap. It binds but removes '
					'little. Correct any dosage share by the suppression '
					'above and treat per-copy comparisons as good to a few '
					'percent.' % name)
			else:
				print('  VERDICT: %s is CENSORED. Its realised per-copy rate '
					'is set by the RNAP footprint rather than by the assigned '
					'probability, so per-copy comparisons between positions '
					'are truncated and any dosage share is an upper bound '
					'only.' % name)

	def _plot(self, plot_out_dir, plot_out_file_name, rows, metadata):
		classes = sorted({r['target_class'] for r in rows})
		fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
		colours = {'construct': '#c2410c', 'rrna': '#4a3aa7'}

		for name in classes:
			subset = [r for r in rows if r['target_class'] == name]
			v = [r['variant'] for r in subset]
			colour = colours.get(name, '#767469')

			axes[0].plot(v, [r['demand_per_promoter'] for r in subset], 'o-',
				color=colour, label='%s demand' % name)
			axes[0].plot(v, [r['actual_per_promoter'] for r in subset], 's--',
				color=colour, alpha=0.55, label='%s realised' % name)
			axes[1].plot(v, [r['suppression'] for r in subset], 'o-',
				color=colour, label=name)
			axes[2].plot(v, [r['overcrowded_frac'] for r in subset], 'o-',
				color=colour, label='%s flagged' % name)
			axes[2].plot(v, [r['frac_demand_over_cap'] for r in subset],
				's--', color=colour, alpha=0.55,
				label='%s own demand > cap' % name)

		max_p = [r['max_p'] for r in rows if r['target_class'] == classes[0]]
		v0 = [r['variant'] for r in rows if r['target_class'] == classes[0]]
		axes[0].plot(v0, max_p, ':', color='#767469', label='max_p (cap)')
		axes[0].set_ylabel('probability per promoter')
		axes[0].set_title('Assigned vs surviving')

		axes[1].axhline(0.0, color='#767469', lw=0.8)
		axes[1].set_ylabel('1 - actual / target')
		axes[1].set_title('Fraction the clamp removes')

		axes[2].set_ylabel('fraction of timesteps')
		axes[2].set_title('Flagged vs genuinely over the cap')

		for ax in axes:
			ax.set_xlabel('variant')
			ax.legend(fontsize=7, frameon=False)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_file_name, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
